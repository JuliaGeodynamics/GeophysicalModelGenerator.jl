# Import of InSAR ground-motion data from the European Ground Motion Service (EGMS),
# https://egms.land.copernicus.eu

# --- ETRS89-LAEA Europe (EPSG:3035) -------------------------------------------------------
# EGMS gives its coordinates as easting/northing in this projection. It is a Lambert
# azimuthal equal-area projection on the GRS80 ellipsoid; the formulas are those of EPSG
# guidance note 7-2 (method 9820), written out here so that no projection library is needed.

const LAEA_a = 6378137.0                    # GRS80 semi-major axis [m]
const LAEA_e2 = (2 - 1 / 298.257222101) / 298.257222101
const LAEA_e = sqrt(LAEA_e2)
const LAEA_lat0 = deg2rad(52.0)
const LAEA_lon0 = deg2rad(10.0)
const LAEA_FE = 4321000.0                   # false easting [m]
const LAEA_FN = 3210000.0                   # false northing [m]

laea_q(sinφ) = (1 - LAEA_e2) * (sinφ / (1 - LAEA_e2 * sinφ^2) - log((1 - LAEA_e * sinφ) / (1 + LAEA_e * sinφ)) / (2 * LAEA_e))

const LAEA_qP = laea_q(1.0)
const LAEA_β0 = asin(laea_q(sin(LAEA_lat0)) / LAEA_qP)
const LAEA_Rq = LAEA_a * sqrt(LAEA_qP / 2)
const LAEA_D = LAEA_a * cos(LAEA_lat0) / (sqrt(1 - LAEA_e2 * sin(LAEA_lat0)^2) * LAEA_Rq * cos(LAEA_β0))

"""
    lon, lat = laea_to_lonlat(easting, northing)

Converts ETRS89-LAEA Europe (EPSG:3035) coordinates in meters, as used by EGMS, to
longitude and latitude in degrees.
"""
function laea_to_lonlat(easting::Real, northing::Real)
    x, y = easting - LAEA_FE, northing - LAEA_FN
    ρ = sqrt((x / LAEA_D)^2 + (LAEA_D * y)^2)
    ρ == 0 && return rad2deg(LAEA_lon0), rad2deg(LAEA_lat0)
    C = 2 * asin(ρ / (2 * LAEA_Rq))
    sinβ0, cosβ0 = sincos(LAEA_β0)
    β = asin(cos(C) * sinβ0 + LAEA_D * y * sin(C) * cosβ0 / ρ)
    λ = LAEA_lon0 + atan(x * sin(C), LAEA_D * ρ * cosβ0 * cos(C) - LAEA_D^2 * y * sinβ0 * sin(C))
    e2, e4, e6 = LAEA_e2, LAEA_e2^2, LAEA_e2^3
    φ = β + (e2 / 3 + 31 * e4 / 180 + 517 * e6 / 5040) * sin(2β) +
        (23 * e4 / 360 + 251 * e6 / 3780) * sin(4β) + (761 * e6 / 45360) * sin(6β)
    return rad2deg(λ), rad2deg(φ)
end

"""
    easting, northing = lonlat_to_laea(lon, lat)

Converts longitude and latitude in degrees to ETRS89-LAEA Europe (EPSG:3035) coordinates in
meters. The inverse of [`laea_to_lonlat`](@ref).
"""
function lonlat_to_laea(lon::Real, lat::Real)
    β = asin(laea_q(sind(lat)) / LAEA_qP)
    dλ = deg2rad(lon) - LAEA_lon0
    sinβ0, cosβ0 = sincos(LAEA_β0)
    B = LAEA_Rq * sqrt(2 / (1 + sinβ0 * sin(β) + cosβ0 * cos(β) * cos(dλ)))
    easting = LAEA_FE + B * LAEA_D * cos(β) * sin(dλ)
    northing = LAEA_FN + B / LAEA_D * (cosβ0 * sin(β) - sinβ0 * cos(β) * cos(dλ))
    return easting, northing
end

# --- Reading EGMS files -------------------------------------------------------------------

using Dates, NearestNeighbors, LinearAlgebra

export import_EGMS, EGMS_velocity, EGMS_points, EGMS_on_grid, EGMSData

"""
    EGMSData

InSAR measurement points of the European Ground Motion Service, as read by
[`import_EGMS`](@ref).

- `lon`, `lat`     : position of every point [degrees]
- `height`         : height of every point above the geoid [m]
- `dates`          : the acquisition dates of the time series
- `displacement`   : displacement [mm], of which only the change with time is meaningful; `displacement[i, p]`
                     is the value at `dates[i]` for point `p`
- `fields`         : the other numerical columns of the file, such as `mean_velocity`
                     [mm/yr] (the velocity EGMS computed over the full period)
"""
struct EGMSData
    lon::Vector{Float64}
    lat::Vector{Float64}
    height::Vector{Float64}
    dates::Vector{Date}
    displacement::Matrix{Float32}
    fields::NamedTuple
end

Base.length(d::EGMSData) = length(d.lon)

function Base.show(io::IO, d::EGMSData)
    println(io, "EGMSData ")
    println(io, "  points    : $(length(d))")
    if length(d) > 0
        println(io, "  lon       ϵ [ $(minimum(d.lon)) : $(maximum(d.lon))]")
        println(io, "  lat       ϵ [ $(minimum(d.lat)) : $(maximum(d.lat))]")
    end
    println(io, "  dates     : $(length(d.dates)) from $(first(d.dates)) to $(last(d.dates))")
    return println(io, "  fields    : $(keys(d.fields))")
end

"""
    data = import_EGMS(file; limits = nothing)
    data = import_EGMS(files::Vector; limits = nothing)
    data = import_EGMS(directory; limits = nothing)

Reads one or more `.csv` files of the [European Ground Motion Service](https://egms.land.copernicus.eu)
(EGMS), which provides InSAR-derived ground motion for Europe, into an [`EGMSData`](@ref)
structure. Unzip the files that you downloaded from the EGMS Explorer first. If a directory
is given, all `EGMS_*.csv` files in it are read.

All three product levels can be read: *Basic* (L2a) and *Calibrated* (L2b), which give the
motion along the line of sight of the satellite, and *Ortho* (L3), which gives the vertical
(`_U_` in the file name) or east-west (`_E_`) component on a 100 m grid.

Optionally, only the points within `limits = [lon_min, lon_max, lat_min, lat_max]` are kept.

Several files can only be combined if they have the same acquisition dates, which is the
case for the Ortho tiles of one release and one component.

# Example
```julia
julia> using Dates
julia> limits = [14.7, 15.3, 37.5, 38.0]            # Etna
julia> data = import_EGMS("EGMS_L3_E47N16_100km_U_2020_2025_1.csv", limits = limits)
julia> vel  = EGMS_velocity(data, from = Date(2021, 1, 1), to = Date(2021, 12, 31))  # mm/yr
julia> Topo = import_topo(limits, res = "03s")
julia> Surf = EGMS_on_grid(data, Topo, from = Date(2021, 1, 1), to = Date(2021, 12, 31), smooth = true)
julia> write_paraview(Surf, "Etna_EGMS")
```
See also [`EGMS_velocity`](@ref), [`EGMS_points`](@ref) and [`EGMS_on_grid`](@ref).
"""
function import_EGMS(files::AbstractVector{<:AbstractString}; limits = nothing)
    isempty(files) && error("No EGMS files given")
    data = [read_EGMS_file(f, limits) for f in files]
    d = first(data)
    length(data) == 1 && return d
    for (f, o) in zip(files, data)
        o.dates == d.dates || error("$f has different acquisition dates than $(first(files)); only files of the same release, level and component (and for L2 the same track) can be combined")
        keys(o.fields) == keys(d.fields) || error("$f has different columns than $(first(files))")
    end
    fields = NamedTuple{keys(d.fields)}(Tuple(reduce(vcat, [o.fields[k] for o in data]) for k in keys(d.fields)))
    return EGMSData(
        reduce(vcat, [o.lon for o in data]), reduce(vcat, [o.lat for o in data]), reduce(vcat, [o.height for o in data]),
        d.dates, reduce(hcat, [o.displacement for o in data]), fields
    )
end

function import_EGMS(path::AbstractString; limits = nothing)
    isdir(path) || return import_EGMS([path]; limits = limits)
    files = filter(f -> startswith(basename(f), "EGMS_") && endswith(lowercase(f), ".csv"), readdir(path, join = true))
    isempty(files) && error("No EGMS_*.csv files found in $path; unzip the downloaded files first")
    return import_EGMS(files; limits = limits)
end

# Reads a single file. The columns are found by their name in the header, so that the
# different levels (which have different columns) are all handled by the same code.
function read_EGMS_file(file::AbstractString, limits)
    isfile(file) || error("File $file not found")
    endswith(lowercase(file), ".zip") && error("$file is a zip archive; unzip it first and pass the .csv file")

    return open(file) do io
        header = readline(io)
        delim = count(==(';'), header) > count(==(','), header) ? ';' : ','
        names = lowercase.(strip.(split(header, delim)))

        idates = findall(n -> occursin(r"^\d{8}$", n), names)
        isempty(idates) && error("$file does not look like an EGMS file: no date columns (such as 20190105) in the header")
        dates = [Date(n, dateformat"yyyymmdd") for n in names[idates]]

        col(name) = findfirst(==(name), names)
        ilon, ilat, ieast, inorth = col("longitude"), col("latitude"), col("easting"), col("northing")
        has_lonlat = !isnothing(ilon) && !isnothing(ilat)
        has_lonlat || (!isnothing(ieast) && !isnothing(inorth)) ||
            error("$file has neither longitude/latitude nor easting/northing columns")
        iheight = something(col("height_ortho"), col("height"), Some(nothing))     # name in new and in old releases

        # all remaining columns are kept if they hold numbers (the point id, for example, does not)
        skip = Set(filter(!isnothing, [ilon, ilat, ieast, inorth, iheight]))
        iother = [i for i in eachindex(names) if !(i in skip) && !(i in idates)]
        first_row = position(io)
        row = split(readline(io), delim)
        seek(io, first_row)
        iother = [i for i in iother if i <= length(row) && !isnothing(tryparse(Float64, row[i]))]

        lon, lat, height = Float64[], Float64[], Float64[]
        other = [Float64[] for _ in iother]
        disp = Float32[]
        num(s) = something(tryparse(Float64, s), NaN)

        for line in eachline(io)
            isempty(strip(line)) && continue
            row = split(line, delim)
            length(row) == length(names) || error("$file: a row has $(length(row)) columns, but the header has $(length(names))")
            lo, la = has_lonlat ? (num(row[ilon]), num(row[ilat])) : laea_to_lonlat(num(row[ieast]), num(row[inorth]))
            if !isnothing(limits)
                (limits[1] <= lo <= limits[2] && limits[3] <= la <= limits[4]) || continue
            end
            push!(lon, lo); push!(lat, la)
            push!(height, isnothing(iheight) ? 0.0 : num(row[iheight]))
            for (k, i) in enumerate(iother)
                push!(other[k], num(row[i]))
            end
            for i in idates
                push!(disp, something(tryparse(Float32, row[i]), NaN32))
            end
        end

        fields = NamedTuple{Tuple(Symbol.(names[iother]))}(Tuple(other))
        EGMSData(lon, lat, height, dates, reshape(disp, length(dates), length(lon)), fields)
    end
end

# --- Velocities ---------------------------------------------------------------------------

"""
    vel = EGMS_velocity(data::EGMSData; from = first(data.dates), to = last(data.dates), seasonal = nothing)

Computes the average velocity [mm/yr] of every point between the dates `from` and `to`, by
fitting a straight line through the displacement time series within that period. Points
with too few measurements in the period get `NaN`.

Many points move up and down with the seasons, which would bias the velocity, in particular
over short periods. With `seasonal = true` a yearly signal is therefore fitted together with
the straight line, as EGMS does for its own `mean_velocity` (`data.fields.mean_velocity`),
which covers the full period of the file. That needs at least a year of data, so by default
it is done if the period is at least one year long.

# Example
```julia
julia> using Dates
julia> vel = EGMS_velocity(data, from = Date(2021, 1, 1), to = Date(2022, 12, 31))
```
"""
function EGMS_velocity(data::EGMSData; from = first(data.dates), to = last(data.dates), seasonal = nothing)
    ind = findall(d -> Date(from) <= d <= Date(to), data.dates)
    length(ind) >= 2 || error("Fewer than two acquisition dates between $from and $to; the data covers $(first(data.dates)) to $(last(data.dates))")
    t = [Dates.value(d - data.dates[ind[1]]) / 365.25 for d in data.dates[ind]]     # years
    seasonal = something(seasonal, t[end] >= 1)
    seasonal && t[end] < 1 && error("Fitting a seasonal signal needs at least one year of data; use `seasonal = false` or a longer period")

    # displacement = offset + velocity*t (+ yearly sine and cosine); the velocity is the 2nd unknown
    A = seasonal ? [ones(length(t)) t sinpi.(2 .* t) cospi.(2 .* t)] : [ones(length(t)) t]
    w = pinv(A)[2, :]                       # weights that give the velocity if no value is missing

    vel = fill(NaN, length(data))
    for p in eachindex(vel)
        d = view(data.displacement, ind, p)
        if !any(isnan, d)
            vel[p] = sum(w[k] * d[k] for k in eachindex(w))
        else
            ok = findall(!isnan, d)
            length(ok) > size(A, 2) && (vel[p] = (A[ok, :] \ Float64.(d[ok]))[2])
        end
    end
    return vel
end

"""
    Points = EGMS_points(data::EGMSData; from = first(data.dates), to = last(data.dates), seasonal = nothing)

Returns the measurement points as a `GeoData` set of points, located at their own height,
with the `velocity` [mm/yr] between `from` and `to` (see [`EGMS_velocity`](@ref)) and the
other columns of the file as fields. Save it with `write_paraview(Points, "name", PointsData = true)`.
"""
function EGMS_points(data::EGMSData; kwargs...)
    fields = merge((velocity = EGMS_velocity(data; kwargs...),), data.fields)
    return GeoData(data.lon, data.lat, data.height ./ 1000, fields)
end

# --- Interpolation on a grid --------------------------------------------------------------

"""
    Grid_vel = EGMS_on_grid(data::EGMSData, Grid::GeoData; from = first(data.dates), to = last(data.dates), seasonal = nothing, radius = nothing, smooth = false, minpoints = 10, east = nothing)
    Grid_vel = EGMS_on_grid(data::EGMSData, lon::AbstractRange, lat::AbstractRange; kwargs...)

Puts the velocity [mm/yr] between the dates `from` and `to` (see [`EGMS_velocity`](@ref)) on
a surface `Grid`, for example a topography that was loaded with `import_topo`, or on a
regular grid at sea level that is spanned by the ranges `lon` and `lat`.

Every grid point gets the average velocity of the measurement points that are within
`radius` [km] of it, and `NaN` if there are none: InSAR has no data over water, vegetation
or snow, and those gaps are kept rather than filled. By default `radius` is half the
diagonal of a grid cell, but at least 0.1 km (the spacing of the EGMS Ortho product).

With `smooth = true` the points are weighted with a Gaussian of their distance (of width
`radius/3`) instead of equally, which gives a smooth field; `radius` is then three times
larger by default. Being a weighted average, the result never exceeds the range of the data.

Where there is a lot of data, such a `radius` would average over more than needed and blur
sharp steps in the velocity, such as across a fault. The radius is therefore reduced locally
to the smallest one that still holds `minpoints` measurement points, but not below the size
of a grid cell. So `radius` is only used in full where data is sparse. Set `minpoints = Inf`
to always use the full `radius`.

`Grid_vel` is `Grid` with the additional fields `velocity`, `npoints` (the number of
measurement points that were averaged) and the vector `Velocity = (Ve, Vn, Vz)`, which can
be shown as arrows in Paraview with the `Glyph` filter. For this, `data` should be the
vertical component (the `_U_` file of the Ortho level), which becomes `Vz`. Pass the
east-west component (the `_E_` file) as `east` to fill `Ve`, which is also added as the field
`velocity_east`; without it `Ve` is zero. `Vn` is always zero, as EGMS has no north-south
component.
"""
function EGMS_on_grid(data::EGMSData, Grid::GeoData; radius = nothing, smooth = false, minpoints = 10, east = nothing, kwargs...)
    size(Grid.lon.val, 3) == 1 || error("Grid should be a surface (a single level in the third direction)")
    vel = EGMS_velocity(data; kwargs...)
    lon, lat = Grid.lon.val, Grid.lat.val

    # local flat coordinates in km, which is accurate enough to find neighbours
    kmdeg = 111.195
    c = cosd(sum(lat) / length(lat))
    dx = size(lon, 1) > 1 ? kmdeg * c * sum(abs, diff(lon, dims = 1)) / length(diff(lon, dims = 1)) : 0.0
    dy = size(lat, 2) > 1 ? kmdeg * sum(abs, diff(lat, dims = 2)) / length(diff(lat, dims = 2)) : 0.0
    rcell = max(hypot(dx, dy) / 2, 0.1)
    radius = something(radius, rcell * (smooth ? 3 : 1))

    keep = findall(!isnan, vel)
    velocity = fill(NaN, size(lon))
    npoints = zeros(Int64, size(lon))
    if !isempty(keep)
        P = permutedims([kmdeg * c .* data.lon[keep] kmdeg .* data.lat[keep]])
        tree = KDTree(P)
        for i in eachindex(lon)
            x, y = kmdeg * c * lon[i], kmdeg * lat[i]
            r = radius
            idx = inrange(tree, [x, y], r)
            isempty(idx) && continue
            if length(idx) > minpoints && r > rcell     # plenty of data: shrink the radius, to keep sharp features
                r = min(r, max(rcell, maximum(knn(tree, [x, y], minpoints)[2])))
                idx = inrange(tree, [x, y], r)
            end
            npoints[i] = length(idx)
            w = smooth ? [exp(-((P[1, j] - x)^2 + (P[2, j] - y)^2) / (r / 3)^2) for j in idx] : ones(length(idx))
            velocity[i] = sum(w[k] * vel[keep[j]] for (k, j) in enumerate(idx)) / sum(w)
        end
    end
    fields = (velocity = velocity, npoints = npoints)
    Ve = zeros(size(lon))
    if !isnothing(east)
        Ve = EGMS_on_grid(east, Grid; radius = radius, smooth = smooth, minpoints = minpoints, kwargs...).fields.velocity
        fields = merge(fields, (velocity_east = Ve,))
    end
    # copies, as writing to Paraview rotates the vector components in place
    fields = merge(fields, (Velocity = (copy(Ve), zeros(size(lon)), copy(velocity)),))
    return GeoData(Grid.lon.val, Grid.lat.val, Grid.depth.val, merge(Grid.fields, fields), Grid.atts)
end

function EGMS_on_grid(data::EGMSData, lon::AbstractRange, lat::AbstractRange; kwargs...)
    Lon, Lat, Depth = lonlatdepth_grid(lon, lat, 0.0)
    return EGMS_on_grid(data, GeoData(Lon, Lat, Depth, (elevation = zeros(size(Lon)),)); kwargs...)
end
