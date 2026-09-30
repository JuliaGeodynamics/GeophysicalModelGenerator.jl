# Downloads topography from the GMT data server without going through GMT itself.
#
# The server publishes its grids in two forms: the coarse resolutions (06m and coarser) as
# single NetCDF files, and the finer ones as tiles of JPEG2000. Both are plain HTTP
# downloads, so the only thing GMT was really needed for was decoding -- which `OpenJpeg_jll`
# does. That matters because GMT pulls in GDAL_jll, whose HDF5 bounds conflict with the
# PETSc-based stack that LaMEM uses, and because GMT's own share directory is not reliably
# relocatable.

using Downloads, OpenJpeg_jll, Scratch

"""
How long to wait for one tile before giving up on that mirror and trying the next. The
transfers are a few megabytes, so a minute is generous; without a limit a throttled or
stalled connection blocks indefinitely.
"""
const TILE_TIMEOUT = 60.0

"""
The GMT data server mirrors. The first that answers is used; the others are tried in turn if
a download fails, so one server being offline does not fail the import.
"""
const GMT_SERVERS = [
    "https://oceania.generic-mapping-tools.org/server/earth",
    "https://data.generic-mapping-tools.org/server/earth",
    "https://au.generic-mapping-tools.org/server/earth",
    "https://sa.generic-mapping-tools.org/server/earth",
]

"""
How long to wait for the server's index file, which is a hundred kilobytes: if that does not
arrive within a few seconds the mirror is not going to deliver the tiles either.
"""
const INDEX_TIMEOUT = 15.0

"""
Once every mirror has failed to deliver a file, how long (seconds) to leave the server alone
before asking it again. Without this, every later call in the session -- and in every other
test worker, hence the marker file rather than a variable -- waits out the full retry budget
too, and a run that should fail in a minute takes hours.
"""
const SERVER_RETRY_AFTER = 600.0

"""
    tile_size_deg(res)

The width and height, in degrees, of one tile at resolution `res`, or `nothing` if that
resolution is published as a single global file rather than as tiles.
"""
function tile_size_deg(dataset::AbstractString, res::AbstractString, reg::AbstractString)
    e = server_entry(dataset, res, reg)
    isnothing(e) && return nothing                      # not tiled: one file for the globe
    endswith(e.name, ".grd") && return nothing
    return e.tile_deg
end

"""
    server_entry(dataset, res, reg)

The line that `gmt_data_server.txt` holds for this set, as a named tuple. Everything that
varies from set to set -- how large a tile is, how the integers in it scale to metres, which
coarser set fills the gaps -- is published there, so it is read rather than assumed.
"""
function server_entry(dataset::AbstractString, res::AbstractString, reg::AbstractString)
    index = server_index_file()
    isnothing(index) && return nothing
    want = "$(dataset)_$(res)_$(reg)"
    for line in eachline(index)
        startswith(line, "#") && continue
        f = split(line, '\t')
        length(f) < 11 && continue
        name = rstrip(f[2], ['/'])
        replace(name, ".grd" => "") == want || continue
        return (name = name,
                scale = something(tryparse(Float64, f[5]), 1.0),
                offset = something(tryparse(Float64, f[6]), 0.0),
                tile_deg = something(tryparse(Int, f[8]), 0),
                filler = strip(f[11]))
    end
    return nothing
end

"""
The copy of the server's index that ships with the package, for when the server cannot be
asked for the current one. It changes rarely (new sets are added now and then; the existing
lines stay), so a stale copy is far better than none: without it a tiled resolution would be
looked for as a single global file, which does not exist.
"""
const BUNDLED_SERVER_INDEX = joinpath(@__DIR__, "..", "assets", "gmt_data_server.txt")

# whether this session has already tried to download the index
const SERVER_INDEX_TRIED = Ref(false)

"""
    server_index_file()

The server's own index of its datasets: downloaded once into the cache and read from there
after that. If it cannot be downloaded, the copy bundled with the package is used instead.
The download is attempted once per session, not once per call: the index is consulted
several times per import, and with the server down each attempt would cost a full round of
the mirrors.
"""
function server_index_file()
    index = joinpath(topo_cache_dir(), "gmt_data_server.txt")
    (isfile(index) && filesize(index) > 0) && return index
    if !SERVER_INDEX_TRIED[] && server_down_for() == 0
        SERVER_INDEX_TRIED[] = true
        for server in GMT_SERVERS
            root = replace(server, "/server/earth" => "")
            try
                Downloads.download("$root/gmt_data_server.txt", index; timeout = INDEX_TIMEOUT)
                return index
            catch
                rm(index; force = true)               # no half-written index
            end
        end
    end
    return BUNDLED_SERVER_INDEX
end

"""
    server_down_for()

How many seconds are left of leaving the server alone after every mirror failed to deliver a
file, or `0.0` if it may be asked. The marker is a file in the cache rather than a variable
so that separate processes -- the parallel test workers -- share it.
"""
function server_down_for()
    marker = joinpath(topo_cache_dir(), "server_unreachable_until")
    isfile(marker) || return 0.0
    until = something(tryparse(Float64, read(marker, String)), 0.0)
    remaining = until - time()
    remaining > 0 && return remaining
    rm(marker; force = true)
    return 0.0
end

function mark_server_down()
    return write(joinpath(topo_cache_dir(), "server_unreachable_until"), string(time() + SERVER_RETRY_AFTER))
end

"""
    reset_topo_server()

Forget that the GMT data server was unreachable, so that the next `import_topo` asks it
again right away rather than after `SERVER_RETRY_AFTER` seconds.
"""
reset_topo_server() = rm(joinpath(topo_cache_dir(), "server_unreachable_until"); force = true)

"""
    tile_size_px(res, deg, reg)

How many pixels a tile has along each side: the tile width in degrees divided by the spacing
that `res` names (`m` = arcminutes, `s` = arcseconds), plus one for gridline registration.
"""
function tile_size_px(res::AbstractString, deg::Integer, reg::AbstractString)
    unit = last(res)
    n = parse(Int, res[1:(end - 1)])
    spacing_deg = unit == 'm' ? n / 60 : n / 3600
    npix = round(Int, deg / spacing_deg)
    # Pixel registration ("p") samples cell centres, so a tile of `deg` degrees holds exactly
    # `npix` of them. Gridline registration ("g") samples the cell corners instead and so
    # repeats the shared edge: one more row and column. Getting this wrong shifts everything.
    return reg == "g" ? npix + 1 : npix
end

"""
    tile_filename(lat, lon, dataset, res, reg)

The name of the tile whose south-west corner is at (`lat`, `lon`).

The 3 arcsecond set is the odd one out: it is SRTM repackaged, and is named after that rather
than after the dataset it is served under.
"""
function tile_filename(lat::Integer, lon::Integer, dataset::AbstractString,
                       res::AbstractString, reg::AbstractString)
    ns = lat < 0 ? "S" : "N"
    ew = lon < 0 ? "W" : "E"
    corner = string(ns, lpad(abs(lat), 2, '0'), ew, lpad(abs(lon), 3, '0'))
    res in ("03s", "01s") && return "$corner.SRTMGL$(res == "03s" ? 3 : 1).jp2"
    return "$corner.$(dataset)_$(res)_$(reg).jp2"
end

"""
    tile_directory(dataset, res, reg)

The directory on the server that holds the tiles of this set.
"""
function tile_directory(dataset::AbstractString, res::AbstractString, reg::AbstractString)
    res in ("03s", "01s") && return "earth_relief_$(res)_g"
    return "$(dataset)_$(res)_$(reg)"
end

"""
    dataset_scale(dataset, res, reg)

The scale factor and offset that turn the integers in a tile into metres. The server
publishes these in `gmt_data_server.txt`, one line per set: `earth_relief_15s_p` is stored in
half metres, for instance, while the SRTM sets are stored directly in metres. Reading them
rather than assuming keeps this correct if a set is ever re-encoded.

Falls back to `(1, 0)` if the index cannot be fetched, which is right for every set that is
not scaled.
"""
function dataset_scale(dataset::AbstractString, res::AbstractString, reg::AbstractString)
    e = server_entry(dataset, res, reg)
    return isnothing(e) ? (1.0, 0.0) : (e.scale, e.offset)
end

"""
    filler_dataset(dataset, res, reg)

Some sets cover only part of the globe and are published with a coarser set to fill the rest:
the 3 arcsecond relief is SRTM, which is land only, and the server names
`earth_relief_15s_p` as what belongs in the sea. GMT resamples that filler into the gaps, and
so do we -- otherwise the ocean comes out flat at zero.

Returns `nothing` when a set needs no filler.
"""
function filler_dataset(dataset::AbstractString, res::AbstractString, reg::AbstractString)
    e = server_entry(dataset, res, reg)
    isnothing(e) && return nothing
    (isempty(e.filler) || e.filler == "-") && return nothing
    m = match(r"^(.*)_(\d+[dms])_([gp])$", e.filler)
    return isnothing(m) ? nothing : (String(m[1]), String(m[2]), String(m[3]))
end

"""
    native_registration(dataset, res, reg)

Which registration this set is actually published in. Most come in both, and then the
caller's `reg` is used; the SRTM-derived sets are gridline only and 15s is pixel only, so
asking for the other one would 404 on every tile.
"""
function native_registration(dataset::AbstractString, res::AbstractString, reg::AbstractString)
    res in ("03s", "01s") && return "g"
    index = server_index_file()
    have = String[]
    for line in eachline(index)
        startswith(line, "#") && continue
        f = split(line, '\t')
        length(f) < 4 && continue
        name = replace(rstrip(f[2], ['/']), ".grd" => "")
        startswith(name, "$(dataset)_$(res)_") && push!(have, String(strip(f[4])))
    end
    isempty(have) && return reg
    return reg in have ? reg : first(have)
end

"""
    subset_indices(coords, lo, hi, px)

The samples that make up `lo..hi`, following what GMT returns. With gridline registration the
samples are the cell corners and GMT widens the request out to the nodes that enclose it;
with pixel registration they are the cell centres and a cell belongs to the region when its
centre does.
"""
function subset_indices(coords::AbstractVector, lo::Real, hi::Real, px::Real, reg::AbstractString)
    if reg == "g"
        # nodes sit on the cell corners, and GMT widens the request out to the nodes that
        # enclose it -- but no further, so the tolerance stops just short of a full step
        tol = px - 1.0e-9
        return findall(c -> (lo - tol) < c < (hi + tol), coords)
    else
        # nodes are cell centres, and a cell belongs to the region when its centre does
        tol = 1.0e-9
        return findall(c -> (lo - tol) <= c <= (hi + tol), coords)
    end
end

"""
    topo_cache_dir()

Where downloaded tiles are kept, so that repeated imports of the same region do not download
them again. Uses a scratch space, which is removed when the package is.
"""
topo_cache_dir() = @get_scratch!("topo_tiles")

"""
    read_topo_netcdf(file)

Read a global `.grd` grid, returning `(lon, lat, z)`. Filled in by the `NCDatasets`
extension; the tiled resolutions do not go through here.
"""
function read_topo_netcdf(file::AbstractString)
    return error("""
        reading the global grid $(basename(file)) needs NCDatasets. Either

            using NCDatasets

        or ask for one of the tiled resolutions (05m and finer), which need nothing extra.
        """)
end

"""
    download_tile(url, dest)

Fetch one file, trying each mirror in turn. Returns `true` if it was downloaded, `false` if
the file is not on the server -- which is not an error: the 3 arcsecond set covers land only,
so a tile that is entirely ocean simply does not exist.
"""
function download_tile(relative_url::AbstractString, dest::AbstractString;
                      maxattempts::Integer = 5, timeout::Real = TILE_TIMEOUT)
    isfile(dest) && filesize(dest) > 0 && return true       # cached

    wait = server_down_for()
    wait > 0 && error("the GMT data server could not be reached a moment ago, so $(basename(dest)) was not requested; it will be tried again in $(round(Int, wait)) s. To retry now, call GeophysicalModelGenerator.reset_topo_server().")

    tmp = dest * ".part"                                    # never leave a half file behind
    for attempt in 1:max(1, maxattempts)
        for server in GMT_SERVERS
            try
                # `Downloads.download` waits forever by default. The data server throttles
                # when several jobs pull tiles at once -- which is exactly what CI does --
                # and without a timeout a throttled connection hangs the whole run rather
                # than failing over to the next mirror.
                Downloads.download("$server/$relative_url", tmp; timeout = timeout)
                mv(tmp, dest; force = true)
                return true
            catch err
                rm(tmp, force = true)
                # a 404 means there is no such tile, and no mirror will have it either
                if err isa Downloads.RequestError && err.response.status == 404
                    return false
                end
            end
        end
        # back off a little before going round the mirrors again, so a server that is busy
        # is given a chance rather than hammered
        attempt < maxattempts && sleep(min(2.0^attempt, 10.0))
    end
    # Every mirror failed every time: the server is down, or throttling us. That is not the
    # same as a tile that does not exist (which is filled with sea level), so it is an error
    # -- and it is remembered, so the next call does not wait all of this out again.
    mark_server_down()
    return error("could not download $relative_url from any GMT data server mirror in $maxattempts attempts")
end

"""
    decode_jp2(file, n)

Decode a JPEG2000 tile into an `n x n` matrix of `Int16` elevations, using the decoder that
comes with `OpenJpeg_jll`. The file stores its rows from north to south.
"""
function decode_jp2(file::AbstractString, n::Integer)
    raw = file * ".raw"
    if !(isfile(raw) && filesize(raw) == n * n * 2)
        run(pipeline(`$(opj_decompress()) -i $file -o $raw`, stdout = devnull, stderr = devnull))
    end
    A = Array{Int16}(undef, n, n)
    read!(raw, A)
    return permutedims(A)                # to (row, column) = (lat, lon)
end

"""
    grid_from_tiles(limits, dataset, res, reg)

Assemble the tiles that cover `limits` and cut out the requested region. Tiles that are not
on the server are taken to be at sea level, which is what the land-only 3 arcsecond set needs.
"""
function grid_from_tiles(limits, dataset::AbstractString, res::AbstractString, reg::AbstractString;
                         maxattempts::Integer = 5)
    lonmin, lonmax, latmin, latmax = limits
    # Not every set is published in both registrations: the SRTM-derived ones are gridline
    # only, 15s is pixel only. GMT serves whichever the set actually has rather than
    # resampling, so follow the set and only honour `reg` where there is a choice.
    tile_reg = native_registration(dataset, res, reg)
    deg = tile_size_deg(dataset, res, tile_reg)
    n = tile_size_px(res, deg, tile_reg)
    px = deg / (tile_reg == "g" ? n - 1 : n)         # degrees per pixel

    # The tiles are named after their south-west corner. Longitude tiles start at the
    # antimeridian and latitude tiles at the south pole, so for a tile larger than 30 degrees
    # the latitude anchors are -90, -90+deg, ... rather than multiples of `deg`.
    snap(v, origin) = origin + Int(floor((v - origin) / deg)) * deg
    lat0 = snap(latmin, -90)
    lat1 = snap(latmax - 1.0e-9, -90)
    lon0 = snap(lonmin, -180)
    lon1 = snap(lonmax - 1.0e-9, -180)

    lats = lat0:deg:lat1
    lons = lon0:deg:lon1

    # Gridline-registered tiles repeat the edge they share with the next tile along, so a
    # tile contributes `step` new rows/columns and only the last one of a row contributes
    # its final edge as well.
    step = tile_reg == "g" ? n - 1 : n
    Z = zeros(Int16, length(lats) * step + (tile_reg == "g"),
                     length(lons) * step + (tile_reg == "g"))
    dir = topo_cache_dir()
    sub = tile_directory(dataset, res, tile_reg)

    missing_tiles = Tuple{Int, Int}[]
    for (j, la) in enumerate(lats), (i, lo) in enumerate(lons)
        name = tile_filename(la, lo, dataset, res, tile_reg)
        file = joinpath(dir, name)
        if !download_tile("$dataset/$sub/$name", file; maxattempts = maxattempts)
            push!(missing_tiles, (la, lo))      # filled from the coarser set below
            continue
        end
        A = decode_jp2(file, n)
        Z[((j - 1) * step + 1):((j - 1) * step + n),
          ((i - 1) * step + 1):((i - 1) * step + n)] = A[end:-1:1, :]
    end

    # the coordinates of each sample: cell centres for pixel registration, cell corners
    # (so starting exactly on the tile edge) for gridline registration
    off = tile_reg == "g" ? 0.0 : 0.5
    lon = [lon0 + (k - 1 + off) * px for k in 1:size(Z, 2)]
    lat = [lat0 + (k - 1 + off) * px for k in 1:size(Z, 1)]
    ix = subset_indices(lon, lonmin, lonmax, px, tile_reg)
    iy = subset_indices(lat, latmin, latmax, px, tile_reg)

    isempty(ix) && error("no data between longitudes $lonmin and $lonmax")
    isempty(iy) && error("no data between latitudes $latmin and $latmax")

    # the tiles hold integers in units the server declares, not metres
    scale, offset = dataset_scale(dataset, res, tile_reg)
    out = Float64.(Z[iy, ix]) .* scale .+ offset

    # where a tile was missing there is no data at all, only the zeros we started with.
    # Fill those from the coarser set that the server names for the purpose, so the sea has
    # its depth rather than being flat.
    filler = filler_dataset(dataset, res, tile_reg)
    if !isnothing(filler)
        # SRTM records the sea as a flat zero rather than leaving it out, and has no tile at
        # all where there is no land, so both cases have to be filled -- every zero, not only
        # the tiles that were missing. GMT does the same, which is why its 03s has bathymetry.
        fill_gaps!(out, lon[ix], lat[iy], filler)
    end

    return lon[ix], lat[iy], out
end

"""
    fill_gaps!(Z, lon, lat, filler)

Replace the cells of `Z` that carry no data -- the flat zeros that SRTM writes for the sea,
and whatever a missing tile left behind -- with values interpolated from the coarser `filler`
set that the server names for the purpose.
"""
function fill_gaps!(Z, lon, lat, filler)
    any(iszero, Z) || return Z
    fdataset, fres, freg = filler
    flon, flat, FZ = grid_from_tiles([minimum(lon), maximum(lon), minimum(lat), maximum(lat)],
                                     fdataset, fres, freg)
    for j in axes(Z, 1), i in axes(Z, 2)
        iszero(Z[j, i]) || continue
        Z[j, i] = bilinear(flon, flat, FZ, lon[i], lat[j])
    end
    return Z
end

"""
    bilinear(xs, ys, A, x, y)

Sample `A` at (`x`, `y`) by interpolating between the four surrounding points, which is what
GMT does when it resamples a coarser set into a finer grid. Falls back to the nearest edge
value outside the grid.
"""
function bilinear(xs::AbstractVector, ys::AbstractVector, A::AbstractMatrix, x::Real, y::Real)
    i = searchsortedlast(xs, x)
    j = searchsortedlast(ys, y)
    i = clamp(i, 1, length(xs) - 1)
    j = clamp(j, 1, length(ys) - 1)

    tx = (x - xs[i]) / (xs[i + 1] - xs[i])
    ty = (y - ys[j]) / (ys[j + 1] - ys[j])
    tx = clamp(tx, 0.0, 1.0)
    ty = clamp(ty, 0.0, 1.0)

    return (1 - tx) * (1 - ty) * A[j, i] + tx * (1 - ty) * A[j, i + 1] +
           (1 - tx) * ty * A[j + 1, i] + tx * ty * A[j + 1, i + 1]
end

"""
    grid_from_single_file(limits, dataset, res, reg)

The coarse resolutions are published as one NetCDF file for the whole globe, which is small
enough to download and cut down here. Reading it needs `NCDatasets`, which is a weak
dependency -- the tiled resolutions, which are the interesting ones, need nothing extra.
"""
function grid_from_single_file(limits, dataset::AbstractString, res::AbstractString, reg::AbstractString;
                               maxattempts::Integer = 5)
    lonmin, lonmax, latmin, latmax = limits
    name = "$(dataset)_$(res)_$(reg).grd"
    file = joinpath(topo_cache_dir(), name)
    download_tile("$dataset/$name", file; maxattempts = maxattempts) ||
        error("could not download $name from the GMT data server")

    lon, lat, Z = read_topo_netcdf(file)
    px = length(lon) > 1 ? abs(lon[2] - lon[1]) : 1.0
    ix = subset_indices(lon, lonmin, lonmax, px, reg)
    iy = subset_indices(lat, latmin, latmax, px, reg)

    isempty(ix) && error("no data between longitudes $lonmin and $lonmax")
    isempty(iy) && error("no data between latitudes $latmin and $latmax")

    return lon[ix], lat[iy], Float64.(Z[iy, ix])
end

"""
    import_topo(limits; dataset="earth_relief", res="01m", reg="p")

Download topography (and bathymetry) for the region `limits = [lonmin, lonmax, latmin, latmax]`
from the GMT data server, and return it as `GeoData` whose `Topography` field is in km.

This does not require `GMT`: the grids are fetched over HTTP and decoded here.

# Datasets
- `"earth_relief"` (the default) blends land topography with ocean bathymetry
- `"earth_gebco"` is GEBCO bathymetry
- `"earth_synbath"` is GEBCO with synthetic bathymetry where the sea floor is unsurveyed

# Resolutions
`"01d"`, `"30m"`, `"20m"`, `"15m"`, `"10m"`, `"06m"` come as one file for the globe;
`"05m"` through `"01m"`, `"30s"`, `"15s"` and `"03s"` are tiled.

Note that `"03s"` is SRTM, which covers **land only** -- the sea is flat zero there, and the
finest resolution that carries real bathymetry is `"15s"`.

# Example
```julia
julia> Topo = import_topo([-18.7, -17.1, 28.0, 29.2], res="03s")
julia> Topo = import_topo(lon=[4,20], lat=[37,49])
```

Tiles are cached, so importing the same region again does not download them a second time.
"""
function import_topo(limits; dataset::AbstractString = "earth_relief",
                     res::AbstractString = "01m", reg::AbstractString = "g",
                     file::Union{Nothing, AbstractString} = nothing,
                     maxattempts::Integer = 5)

    # accept the `file="@earth_relief_01m"` spelling that this function used to take
    if !isnothing(file)
        dataset, res, reg = parse_topo_file(file)
    end

    # longitudes west of Greenwich may be given either way round
    limits = collect(float.(limits))
    if limits[1] > limits[2]
        limits[1], limits[2] = limits[2], limits[1]
    end
    if limits[3] > limits[4]
        limits[3], limits[4] = limits[4], limits[3]
    end

    lon, lat, Z = isnothing(tile_size_deg(dataset, res, native_registration(dataset, res, reg))) ?
        grid_from_single_file(limits, dataset, res, reg; maxattempts = maxattempts) :
        grid_from_tiles(limits, dataset, res, reg; maxattempts = maxattempts)

    Lon, Lat, Depth = lonlatdepth_grid(lon, lat, 0)
    @views Depth[:, :, 1] = 1.0e-3 * Z'            # the field is in km
    return GeoData(Lon, Lat, Depth, (Topography = Depth * km,))
end

import_topo(; lat = [37, 49], lon = [4, 20], kwargs...) =
    import_topo([lon[1], lon[2], lat[1], lat[2]]; kwargs...)

"""
    parse_topo_file(file)

Take apart the `"@earth_relief_01m"` spelling into the dataset, resolution and registration
that name it on the server.
"""
function parse_topo_file(file::AbstractString)
    s = lstrip(file, '@')
    # the tutorials write "@earth_relief_01m.grd"; the suffix names the format GMT used to
    # download into and says nothing about which set is wanted
    endswith(s, ".grd") && (s = s[1:(end - 4)])
    endswith(s, ".nc") && (s = s[1:(end - 3)])
    m = match(r"^(.*)_(\d+[dms])(?:_([gp]))?$", s)
    isnothing(m) && error("cannot make sense of the topography file name \"$file\"")
    # a name without a trailing _g/_p leaves the registration open, and GMT then serves
            # whichever the set has, preferring gridline -- the same default as `reg`
            return String(m[1]), String(m[2]), isnothing(m[3]) ? "g" : String(m[3])
end
