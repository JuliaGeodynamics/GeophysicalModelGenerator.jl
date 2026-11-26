# Gravity anomalies of GeoData (on the reference ellipsoid) and CartData (flat Earth) grids,
# at the grid columns or at arbitrary stations, with adaptive subdivision of cells close to the stations

using NearestNeighbors

const WGS84 = (a = 6378137.0, f = 1 / 298.257223563)

"""
    dg, dgdx, dgdy = voxel_grav(Grid::GeoData, RHO; refMod = "AVG", zRef = nothing, stations = nothing, rhoTol = 1e-9, gradients = false, cutoff = Inf,
                               subdivide = 2.5, ellipsoid = (a = 6378137.0, f = 1 / 298.257223563), backend = nothing, float_type = Float64)

Gravity anomaly [mGal] of the density `RHO` [kg/m³], given at the points of the regular geographic grid `Grid` (as created by `lonlatdepth_grid`), or of a curvilinear grid whose depth levels follow surfaces such as the topography (see [`deform_grid`](@ref)).

The computation is done on the reference ellipsoid (WGS84 by default; `ellipsoid = (a = 6371e3, f = 0)` is a sphere). Depths and heights are heights above the ellipsoid in km; heights relative to sea level are used as such, i.e. the geoid undulation is ignored. Every cell is a point mass at its position on the ellipsoid (subdivided close to the stations), and `dg` is the component of its attraction along the ellipsoid normal at the station (downwards positive), i.e. the gravity disturbance caused by `RHO - refMod`.

- `refMod`:    reference density for every depth level, or `"NE"`, `"SE"`, `"SW"`, `"NW"`, `"AVG"` as in [`voxel_grav`](@ref).
- `zRef`:      height of the stations above every grid column [km]: a number, or an array of size `size(Grid.lon)[1:2]`, e.g. the topography. By default 0 (sea level), or the top of the grid for a curvilinear grid.
- `stations`:  `GeoData` with arbitrary station positions (e.g. measurement points; `depth` is their height in km), instead of the grid columns. `zRef` is then ignored.
- `gradients`: also compute the horizontal gradients of `dg` towards east (`dgdx`) and north (`dgdy`) [mGal/m]; otherwise they are `nothing`.
- `cutoff`:    ignore cells whose horizontal distance from a station is larger than `cutoff` [km]; e.g. 166.7 km is the standard outer radius of terrain corrections. Default: all cells.
- `subdivide`: cells closer to a station than `subdivide` times their size are split recursively into smaller point masses (`0` switches this off).
- `rhoTol`, `backend`, `float_type`: as in [`voxel_grav`](@ref).

On a regular grid, only the part of the model below the highest station is used. The cells of a curvilinear grid are hexahedra with vertical edges, whose inclined top and bottom are taken into account by the full tensor of their second moments. `dg`, `dgdx` and `dgdy` have size `size(Grid.lon)[1:2]`, or the size of `stations.lon`. [`voxel_grav_adjoint`](@ref) gives the adjoint of `dg`.
"""
function voxel_grav(
        Grid::GeoData, RHO::AbstractArray{<:Real, 3}; refMod = "AVG", zRef = nothing, stations::Union{Nothing, GeoData} = nothing, rhoTol = 1.0e-9,
        gradients = false, cutoff = Inf, subdivide = 2.5, ellipsoid = WGS84, backend = nothing, float_type = Float64
    )
    # the fields of GeoData are not concretely typed: unwrap them here and continue in a type-stable function
    st = isnothing(stations) ? nothing : (stations.lon.val, stations.lat.val, stations.depth.val)
    M = isnothing(stations) ? Matrix{promote_type(Float64, eltype(RHO))} : Array{promote_type(Float64, eltype(RHO)), 3}
    r = voxel_grav_ellipsoid(Grid.lon.val, Grid.lat.val, Grid.depth.val, RHO, refMod, zRef, st, rhoTol, gradients, cutoff, subdivide, ellipsoid, backend, float_type)
    return gradients ? r::Tuple{M, M, M} : r::Tuple{M, Nothing, Nothing}
end

function voxel_grav_ellipsoid(Lon, Lat, Depth, RHO, refMod, zRef, st, rhoTol, gradients, cutoff, subdivide, ellipsoid, backend, float_type)
    check_backend(backend)
    size(RHO) == size(Lon) || error("RHO has size $(size(RHO)) but Grid has size $(size(Lon)).")
    geo, win, shape = ellipsoid_setup(Lon, Lat, Depth, zRef, st, ellipsoid, 1.0e3 * cutoff)
    RefMod = reference_model(refMod, RHO)
    DRHO = cell_average(RHO[:, :, geo.ind] .- reshape(RefMod[geo.ind], 1, 1, :))
    return points_forward(geo, win, DRHO, rhoTol, 1.0e3 * cutoff, gradients, subdivide, backend, float_type, shape)
end

"""
    sens = voxel_grav_adjoint(Grid::GeoData, R; zRef = nothing, stations = nothing, cutoff = Inf, subdivide = 2.5, ellipsoid = (a = 6378137.0, f = 1 / 298.257223563),
                              backend = nothing, float_type = Float64)

Adjoint of `voxel_grav(Grid::GeoData, RHO; ...)[1]`: `sens[i, j, k] = ∂(sum(R .* dg)) / ∂RHO[i, j, k]` [mGal per kg/m³], for a fixed reference-model vector. `R` has the size of `dg`. See [`voxel_grav_adjoint`](@ref) for the Cartesian version.
"""
function voxel_grav_adjoint(
        Grid::GeoData, R::AbstractArray; zRef = nothing, stations::Union{Nothing, GeoData} = nothing, cutoff = Inf, subdivide = 2.5,
        ellipsoid = WGS84, backend = nothing, float_type = Float64
    )
    st = isnothing(stations) ? nothing : (stations.lon.val, stations.lat.val, stations.depth.val)
    return voxel_grav_adjoint_ellipsoid(Grid.lon.val, Grid.lat.val, Grid.depth.val, R, zRef, st, cutoff, subdivide, ellipsoid, backend, float_type)::Array{Float64, 3}
end

function voxel_grav_adjoint_ellipsoid(Lon, Lat, Depth, R, zRef, st, cutoff, subdivide, ellipsoid, backend, float_type)
    check_backend(backend)
    geo, win, shape = ellipsoid_setup(Lon, Lat, Depth, zRef, st, ellipsoid, 1.0e3 * cutoff)
    size(R) == shape || error("R has size $(size(R)) but there are $(shape) stations.")
    S = points_adjoint(geo, win, vec(R), 1.0e3 * cutoff, subdivide, backend, float_type)
    return cells_to_nodes(S, size(Lon), geo.ind)
end

"""
    dg, dgdx, dgdy = voxel_grav(Grid::CartData, RHO; stations = nothing, gradients = false, kwargs...)

Flat-Earth gravity anomaly [mGal] and, with `gradients = true`, its horizontal gradients [mGal/m] of the density `RHO` given at the points of the Cartesian grid `Grid` (coordinates in km). Without `stations`, this is `voxel_grav(X, Y, Z, RHO; lengthUnit = "km", gradients, kwargs...)` (`zRef` and `cutoff` in km) without file output. `stations` (`CartData`, in km) gives arbitrary station positions instead of the grid columns; `refMod`, `rhoTol`, `cutoff`, `subdivide`, `backend` and `float_type` work as for `GeoData`.

`Grid` may also be a curvilinear grid whose levels follow surfaces such as the topography (see [`deform_grid`](@ref) and [`curvilinear_grid`](@ref)). Its cells are hexahedra with vertical edges, and the stations are by default on top of the grid columns (`zRef` gives other heights); all cells are used.
"""
function voxel_grav(Grid::CartData, RHO::AbstractArray{<:Real, 3}; stations::Union{Nothing, CartData} = nothing, gradients = false, kwargs...)
    T = promote_type(Float64, eltype(RHO))
    if isnothing(stations) && is_curvilinear(Grid.z.val)
        st, kw = column_stations(Grid, kwargs)
        r = voxel_grav_flat_points(Grid.x.val, Grid.y.val, Grid.z.val, RHO, st, gradients; kw...)
        return gradients ? r::Tuple{Matrix{T}, Matrix{T}, Matrix{T}} : r::Tuple{Matrix{T}, Nothing, Nothing}
    end
    if isnothing(stations)
        kw = merge((; lengthUnit = "km", outName = nothing, printing = false), values(kwargs))
        r = voxel_grav(Grid.x.val, Grid.y.val, Grid.z.val, RHO; gradients, kw...)
        return gradients ? r::Tuple{Matrix{T}, Matrix{T}, Matrix{T}} : r::Tuple{Matrix{T}, Nothing, Nothing}
    end
    r = voxel_grav_flat_points(Grid.x.val, Grid.y.val, Grid.z.val, RHO, (stations.x.val, stations.y.val, stations.z.val), gradients; kwargs...)
    return gradients ? r::Tuple{Array{T, 3}, Array{T, 3}, Array{T, 3}} : r::Tuple{Array{T, 3}, Nothing, Nothing}
end

# stations at the columns of a curvilinear CartData grid, at the heights zRef (by default the top of the grid), and the
# other keyword arguments
function column_stations(Grid, kwargs)
    X, Y, Z = Grid.x.val, Grid.y.val, Grid.z.val
    zo = column_heights(get(kwargs, :zRef, nothing), Z, size(X)[1:2])
    return (X[:, :, 1], Y[:, :, 1], zo), Base.structdiff(values(kwargs), NamedTuple{(:zRef,)})
end

function voxel_grav_flat_points(X, Y, Z, RHO, st, gradients; refMod = "AVG", rhoTol = 1.0e-9, cutoff = Inf, subdivide = 2.5, backend = nothing, float_type = Float64)
    check_backend(backend)
    geo, win, shape = flat_setup(X, Y, Z, st, 1.0e3 * cutoff)
    RefMod = reference_model(refMod, RHO)
    DRHO = cell_average(RHO[:, :, geo.ind] .- reshape(RefMod[geo.ind], 1, 1, :))
    return points_forward(geo, win, DRHO, rhoTol, 1.0e3 * cutoff, gradients, subdivide, backend, float_type, shape)
end

"""
    sens = voxel_grav_adjoint(Grid::CartData, R; stations = nothing, kwargs...)

Adjoint of `voxel_grav(Grid::CartData, RHO; ...)[1]` (flat Earth, coordinates in km).
"""
function voxel_grav_adjoint(Grid::CartData, R::AbstractArray; stations::Union{Nothing, CartData} = nothing, kwargs...)
    if isnothing(stations) && is_curvilinear(Grid.z.val)
        st, kw = column_stations(Grid, kwargs)
        return voxel_grav_adjoint_flat_points(Grid.x.val, Grid.y.val, Grid.z.val, R, st; kw...)::Array{Float64, 3}
    end
    isnothing(stations) && return voxel_grav_adjoint(Grid.x.val, Grid.y.val, Grid.z.val, R; lengthUnit = "km", kwargs...)::Array{Float64, 3}
    return voxel_grav_adjoint_flat_points(Grid.x.val, Grid.y.val, Grid.z.val, R, (stations.x.val, stations.y.val, stations.z.val); kwargs...)::Array{Float64, 3}
end

function voxel_grav_adjoint_flat_points(X, Y, Z, R, st; cutoff = Inf, subdivide = 2.5, backend = nothing, float_type = Float64)
    check_backend(backend)
    geo, win, shape = flat_setup(X, Y, Z, st, 1.0e3 * cutoff)
    size(R) == shape || error("R has size $(size(R)) but there are $(shape) stations.")
    S = points_adjoint(geo, win, vec(R), 1.0e3 * cutoff, subdivide, backend, float_type)
    return cells_to_nodes(S, size(X), geo.ind)
end

# ----------------------------------------------------------------------------------------------------------------------
# geometry

# Earth-centred Cartesian coordinates [m] of longitude and geodetic latitude [°] and height above the ellipsoid h [m]
function geodetic_to_ecef(lon, lat, h, ell)
    e2 = ell.f * (2 - ell.f)
    sφ, cφ = sincosd(lat)
    sλ, cλ = sincosd(lon)
    N = ell.a / sqrt(1 - e2 * sφ^2)
    return (N + h) * cφ * cλ, (N + h) * cφ * sλ, (N * (1 - e2) + h) * sφ
end

# coordinate vectors of a rectilinear grid (first, second, third dimension)
function grid_vectors(A, B, C, msg)
    a, b, c = A[:, 1, 1], B[1, :, 1], C[1, 1, :]
    mono(v) = all(>(0), diff(v)) || all(<(0), diff(v))
    if !(all(A .== a) && all(B .== b') && all(C .== reshape(c, 1, 1, :)) && mono(a) && mono(b) && mono(c))
        error(msg)
    end
    return a, b, c
end

# coordinate vectors of the first two dimensions, and the heights as a vector (rectilinear grid) or as a 3D array
# (curvilinear grid whose levels follow surfaces, e.g. from deform_grid, with monotonic columns)
function grid_axes(A, B, C, msg)
    is_curvilinear(C) || return grid_vectors(A, B, C, msg)
    a, b = A[:, 1, 1], B[1, :, 1]
    mono(v) = all(>(0), diff(v)) || all(<(0), diff(v))
    dz = diff(C, dims = 3)
    all(A .== a) && all(B .== b') && mono(a) && mono(b) && (all(>=(0), dz) || all(<=(0), dz)) || error(msg)
    return a, b, C
end

# station heights above the grid columns: zRef (a number or an array), or by default sea level (rectilinear grid)
# or the top of the grid (curvilinear grid, whose top usually is the topography)
column_heights(::Nothing, z::AbstractVector, n) = zeros(n)
column_heights(::Nothing, z::AbstractArray{<:Any, 3}, n) = dropdims(maximum(z, dims = 3), dims = 3)
function column_heights(zRef, z, n)
    zo = zRef isa Number ? fill(float(zRef), n) : float.(zRef)
    size(zo) == n || error("zRef has size $(size(zo)) but the grid has $n points in the horizontal plane.")
    return zo
end

# Stations (position s, upward normal n, east e and north q unit vectors; heights sh) and cells (centre p, volume vol,
# second moments m1, m2, m3 along east, north, up and, for a curvilinear grid, the mixed moments mo = (m12, m13, m23),
# otherwise nothing; node coordinates bx, by, bz of the cell corners, as lon/lat/height for the ellipsoid (ell isa
# NamedTuple) or x/y/z for the flat Earth (ell = nothing); all lengths in m). On a rectilinear grid (bz a vector), only
# the cells below the highest station are used; on a curvilinear grid (bz a 3D array), all cells. Positions are relative
# to the first station o, which keeps Float32 computations accurate.
function points_geometry(bx, by, bz_all, slon, slat, sh, ell)
    if bz_all isa AbstractVector
        ind = findall(<=(maximum(sh)), bz_all)
        length(ind) >= 2 || error("The model must have at least two depth levels below the highest station.")
        bz = bz_all[ind]
    else
        ind, bz = collect(axes(bz_all, 3)), bz_all
    end
    ns = length(sh)
    if isnothing(ell)
        o = (slon[1], slat[1], sh[1])
        sx, sy, sz = slon .- o[1], slat .- o[2], sh .- o[3]
        nx, ny, nz = zeros(ns), zeros(ns), ones(ns)
        ex, ey = ones(ns), zeros(ns)
        qx, qy, qz = zeros(ns), ones(ns), zeros(ns)
    else
        s = geodetic_to_ecef.(slon, slat, sh, Ref(ell))
        o = s[1]
        sx, sy, sz = getindex.(s, 1) .- o[1], getindex.(s, 2) .- o[2], getindex.(s, 3) .- o[3]
        nx, ny, nz = cosd.(slat) .* cosd.(slon), cosd.(slat) .* sind.(slon), sind.(slat)
        ex, ey = -sind.(slon), cosd.(slon)
        qx, qy, qz = -sind.(slat) .* cosd.(slon), -sind.(slat) .* sind.(slon), cosd.(slat)
    end
    n = (length(bx) - 1, length(by) - 1, size(bz, ndims(bz)) - 1)
    px, py, pz, vol = zeros(n), zeros(n), zeros(n), zeros(n)
    m1, m2, m3 = zeros(n), zeros(n), zeros(n)          # second moments (L²/12 for a box) along east, north, up (x, y, z)
    mo = bz isa AbstractVector ? nothing : (zeros(n), zeros(n), zeros(n))
    for J in CartesianIndices(vol)
        p, v, _, _, m = cell_props(ell, o, cell_box(bx, by, bz, J))
        px[J], py[J], pz[J] = p
        vol[J] = v
        m1[J], m2[J], m3[J] = m[1], m[2], m[3]
        isnothing(mo) || ((mo[1][J], mo[2][J], mo[3][J]) = (m[4], m[5], m[6]))
    end
    return (; ind, ell, o, sx, sy, sz, nx, ny, nz, ex, ey, qx, qy, qz, sh, px, py, pz, vol, m1, m2, m3, mo, bx, by, bz)
end

# a cell of a rectilinear grid as a box (u1, u2, v1, v2, w1, w2), and of a curvilinear grid as a hexahedron
# (u1, u2, v1, v2, zb, zt): its horizontal bounds and the heights of the corners (u1, v1), (u2, v1), (u1, v2), (u2, v2)
# of its bottom (zb) and top (zt)
cell_box(bx, by, bz::AbstractVector, J) = (bx[J[1]], bx[J[1] + 1], by[J[2]], by[J[2] + 1], bz[J[3]], bz[J[3] + 1])
function cell_box(bx, by, bz::AbstractArray{<:Any, 3}, J)
    i, j, k = Tuple(J)
    corners(l) = (bz[i, j, l], bz[i + 1, j, l], bz[i, j + 1, l], bz[i + 1, j + 1, l])
    return (bx[i], bx[i + 1], by[j], by[j + 1], corners(k), corners(k + 1))
end

# centre (relative to o), volume, size, extent and second moments (m11, m22, m33, m12, m13, m23; per unit volume, about
# the centre, along east, north, up) of a cell (x, y, z or lon, lat, height). A box has an analytic solution; a
# hexahedron with vertical edges and inclined top or bottom is integrated exactly with 2×2×2 Gauss points.
function cell_props(ell, o, box::NTuple{6, Real})
    p, v, size, ext = box_props(ell, o, box)
    return p, v, size, ext, (ext[1]^2 / 12, ext[2]^2 / 12, ext[3]^2 / 12, 0.0, 0.0, 0.0)
end

function cell_props(ell, o, (u1, u2, v1, v2, zb, zt))
    allequal(zb) && allequal(zt) && return cell_props(ell, o, (u1, u2, v1, v2, zb[1], zt[1]))
    thick = abs.(zt .- zb)
    if isnothing(ell)
        a, b = u2 - u1, v2 - v1
        vol, c, m = hexahedron_moments(a, b, zb, zt)
        p = ((u1 + u2) / 2 + c[1] - o[1], (v1 + v2) / 2 + c[2] - o[2], c[3] - o[3])
    else
        # in the local east/north/up frame at the centre of the cell
        λ, φ, h = (u1 + u2) / 2, (v1 + v2) / 2, (sum(zb) + sum(zt)) / 8
        Rλ, Rφ = ellipsoid_radii(ell, φ, h)
        a, b = Rλ * deg2rad(u2 - u1), Rφ * deg2rad(v2 - v1)
        vol, c, m = hexahedron_moments(a, b, zb, zt)
        p = geodetic_to_ecef(λ + rad2deg(c[1] / Rλ), φ + rad2deg(c[2] / Rφ), c[3], ell) .- o
    end
    return p, vol, max(abs(a), abs(b), maximum(thick)), (abs(a), abs(b), sum(thick) / 4), m
end

function box_props(::Nothing, o, (u1, u2, v1, v2, w1, w2))
    p = ((u1 + u2) / 2 - o[1], (v1 + v2) / 2 - o[2], (w1 + w2) / 2 - o[3])
    du, dv, dw = abs(u2 - u1), abs(v2 - v1), abs(w2 - w1)
    return p, du * dv * dw, max(du, dv, dw), (du, dv, dw)
end

function box_props(ell::NamedTuple, o, (λ1, λ2, φ1, φ2, h1, h2))
    λ, φ, h = (λ1 + λ2) / 2, (φ1 + φ2) / 2, (h1 + h2) / 2
    p = geodetic_to_ecef(λ, φ, h, ell) .- o
    Rλ, Rφ = ellipsoid_radii(ell, φ, h)
    a, b, c = Rλ * abs(deg2rad(λ2 - λ1)), Rφ * abs(deg2rad(φ2 - φ1)), abs(h2 - h1)
    return p, a * b * c, max(a, b, c), (a, b, c)
end

# radii of the circles of latitude ((N + h) cos φ) and of the meridian (M + h) at latitude φ and height h
function ellipsoid_radii(ell, φ, h)
    e2 = ell.f * (2 - ell.f)
    w = 1 - e2 * sind(φ)^2
    M, N = ell.a * (1 - e2) / w^1.5, ell.a / sqrt(w)
    return (N + h) * cosd(φ), M + h
end

const GAUSS2 = (0.5 - 0.5 / sqrt(3), 0.5 + 0.5 / sqrt(3))     # 2-point Gauss-Legendre nodes on [0, 1]

# value at (ξ, η) of the unit square, interpolated bilinearly from the corner values z = (z11, z21, z12, z22)
bilinear(z, ξ, η) = (1 - η) * ((1 - ξ) * z[1] + ξ * z[2]) + η * ((1 - ξ) * z[3] + ξ * z[4])

# Volume, centroid (x, y relative to the centre of the footprint, z absolute) and second moments about the centroid of
# a hexahedron with vertical edges, a footprint of a × b and the corner heights zb (bottom) and zt (top). The integrands
# have at most degree 3 in each reference coordinate, which 2 Gauss points per direction integrate exactly.
function hexahedron_moments(a, b, zb, zt)
    z0 = (sum(zb) + sum(zt)) / 8                         # reference height, for accuracy
    V = sx = sy = sz = sxx = syy = szz = sxy = sxz = syz = 0.0
    for ξ in GAUSS2, η in GAUSS2
        bot, top = bilinear(zb, ξ, η) - z0, bilinear(zt, ξ, η) - z0
        w = abs(a * b * (top - bot)) / 8
        x, y = (ξ - 0.5) * a, (η - 0.5) * b
        for ζ in GAUSS2
            z = bot + ζ * (top - bot)
            V += w
            sx, sy, sz = sx + w * x, sy + w * y, sz + w * z
            sxx, syy, szz = sxx + w * x * x, syy + w * y * y, szz + w * z * z
            sxy, sxz, syz = sxy + w * x * y, sxz + w * x * z, syz + w * y * z
        end
    end
    V > 0 || return 0.0, (0.0, 0.0, z0), (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    cx, cy, cz = sx / V, sy / V, sz / V
    m = (sxx / V - cx^2, syy / V - cy^2, szz / V - cz^2, sxy / V - cx * cy, sxz / V - cx * cz, syz / V - cy * cz)
    return V, (cx, cy, z0 + cz), m
end

# the part (ξa..ξb, ηa..ηb, ζa..ζb) of a cell, in its reference coordinates, which again is a cell with vertical edges
function subcell((u1, u2, v1, v2, zb, zt), (ξa, ξb), (ηa, ηb), (ζa, ζb))
    corners = ((ξa, ηa), (ξb, ηa), (ξa, ηb), (ξb, ηb))
    b, t = map(c -> bilinear(zb, c...), corners), map(c -> bilinear(zt, c...), corners)
    return (u1 + ξa * (u2 - u1), u1 + ξb * (u2 - u1), v1 + ηa * (v2 - v1), v1 + ηb * (v2 - v1), b .+ ζa .* (t .- b), b .+ ζb .* (t .- b))
end
# a box is halved as on a regular grid
function subcell(box::NTuple{6, Real}, (ξa, ξb), (ηa, ηb), (ζa, ζb))
    (u1, u2, v1, v2, w1, w2) = box
    return (u1 + ξa * (u2 - u1), u1 + ξb * (u2 - u1), v1 + ηa * (v2 - v1), v1 + ηb * (v2 - v1), w1 + ζa * (w2 - w1), w1 + ζb * (w2 - w1))
end

# field of a cell (point mass plus quadrupole), with the mixed moments mo for a curvilinear grid
@inline quad_field(::Nothing, J, c, t1, t2, h, m1, m2, m3) = box_field(c, t1, t2, h, m1, m2, m3)
@inline quad_field(mo, J, c, t1, t2, h, m1, m2, m3) = @inbounds box_field(c, t1, t2, h, m1, m2, m3, mo[1][J], mo[2][J], mo[3][J])

# geometry, station windows and output shape for a GeoData grid, with stations at the grid columns (st = nothing) or at st
function ellipsoid_setup(Lon, Lat, Depth, zRef, st, ell, cutoff)
    lon, lat, z = grid_axes(Lon, Lat, Depth, "Grid must be a regular lon/lat/depth grid (see lonlatdepth_grid), with lon varying along the first, lat along the second and depth along the third dimension, or a curvilinear grid whose levels follow surfaces (see deform_grid).")
    if isnothing(st)
        zo = column_heights(zRef, z, (length(lon), length(lat)))
        slon, slat = vec(lon .+ 0 .* lat'), vec(0 .* lon .+ lat')
        geo = points_geometry(lon, lat, 1.0e3 .* z, slon, slat, 1.0e3 .* vec(zo), ell)
        return geo, grid_windows(geo, lon, lat, cutoff), (length(lon), length(lat))
    end
    geo = points_geometry(lon, lat, 1.0e3 .* z, vec(st[1]), vec(st[2]), 1.0e3 .* vec(st[3]), ell)
    return geo, list_windows(geo, cutoff), size(st[1])
end

# geometry, station windows and output shape for a CartData grid with stations at st = (x, y, z) [km]
function flat_setup(X, Y, Z, st, cutoff)
    x, y, z = grid_axes(X, Y, Z, "Grid must be a regular x/y/z grid, with x varying along the first, y along the second and z along the third dimension, or a curvilinear grid whose levels follow surfaces (see deform_grid).")
    geo = points_geometry(1.0e3 .* x, 1.0e3 .* y, 1.0e3 .* z, 1.0e3 .* vec(st[1]), 1.0e3 .* vec(st[2]), 1.0e3 .* vec(st[3]), nothing)
    return geo, list_windows(geo, cutoff), size(st[1])
end

# ----------------------------------------------------------------------------------------------------------------------
# which stations a cell column can reach within the cutoff: lon/lat index ranges for stations on the grid columns,
# index lists (from a KD-tree) for arbitrary stations; the exact horizontal distance is checked in the kernels

struct GridWindows{A, B}
    rφ::A
    rλ::B
    nλ::Int
end
struct ListWindows{A}
    lists::A
end

station_blocks(w::GridWindows, J) = (w.rλ[J[1], J[2]] .+ (iφ - 1) * w.nλ for iφ in w.rφ[J[2]])
station_blocks(w::ListWindows, J) = (w.lists[J[1], J[2]],)

# stations on the lon/lat grid: the angle γ between station and cell satisfies |Δφ| ≤ γ and |Δλ| cos(φmax) ≤ γ
function grid_windows(geo, lon, lat, cutoff)
    lonc, latc = (lon[1:(end - 1)] .+ lon[2:end]) ./ 2, (lat[1:(end - 1)] .+ lat[2:end]) ./ 2
    nλ, nφ = length(lon), length(lat)
    Rmin = geo.ell.a * (1 - geo.ell.f)^2 + min(0.0, minimum(geo.bz))    # lower bound of the distance to the Earth's centre
    if !isfinite(cutoff) || 1.1 * cutoff >= Rmin
        return GridWindows(fill(1:nφ, length(latc)), fill(1:nλ, length(lonc), length(latc)), nλ)
    end
    γ = rad2deg(asin(1.1 * cutoff / Rmin))                   # degrees, with a 10% margin
    range_of(i) = isempty(i) ? (1:0) : (first(i):last(i))
    rφ = [range_of(findall(φ -> abs(φ - φc) <= γ, lat)) for φc in latc]
    rλ = Matrix{UnitRange{Int}}(undef, length(lonc), length(latc))
    for j in eachindex(latc)
        φmax = min(90.0, maximum(abs, lat[rφ[j]]; init = abs(latc[j])) + γ)
        Δλ = cosd(φmax) > 0 ? γ / cosd(φmax) : 360.0
        for i in eachindex(lonc)
            rλ[i, j] = Δλ >= 180 ? (1:nλ) : range_of(findall(λ -> abs(λ - lonc[i]) <= Δλ, lon))
        end
    end
    return GridWindows(rφ, rλ, nλ)
end

# arbitrary stations: all of them, or those within reach of the column centre (KD-tree; 2D on the flat Earth)
function list_windows(geo, cutoff)
    n = (length(geo.bx) - 1, length(geo.by) - 1)
    isfinite(cutoff) || return ListWindows(fill(1:length(geo.sx), n))
    if isnothing(geo.ell)
        tree = KDTree(permutedims(hcat(geo.sx, geo.sy)))
        lists = [inrange(tree, [geo.px[i, j, 1], geo.py[i, j, 1]], cutoff) for i in 1:n[1], j in 1:n[2]]
    else
        tree = KDTree(permutedims(hcat(geo.sx, geo.sy, geo.sz)))
        reach = 1.1 * cutoff + maximum(geo.sh) - minimum(geo.bz)
        lists = [inrange(tree, [geo.px[i, j, 1], geo.py[i, j, 1], geo.pz[i, j, 1]], reach) for i in 1:n[1], j in 1:n[2]]
    end
    return ListWindows(sort!.(lists))
end

# ----------------------------------------------------------------------------------------------------------------------
# kernels

# anomaly (component along the station normal) and optionally its gradients towards east and north
function points_forward(geo, win, DRHO, rhoTol, cutoff, gradients, subdivide, backend, float_type, shape)
    dg, gx, gy = isnothing(backend) ? point_masses_forward(geo, win, DRHO, rhoTol, cutoff, Val(gradients)) :
        grav_points_forward_ka(backend, float_type, geo, DRHO, rhoTol, cutoff, gradients)
    subdivide > 0 && subdivision_correction!(dg, gx, gy, geo, DRHO, rhoTol, cutoff, subdivide)
    reshape_to(A) = isnothing(A) ? A : reshape(A, shape)
    return reshape_to(dg), reshape_to(gx), reshape_to(gy)
end

function points_adjoint(geo, win, R, cutoff, subdivide, backend, float_type)
    S = isnothing(backend) ? point_masses_adjoint(geo, win, R, cutoff) : reshape(grav_points_adjoint_ka(backend, float_type, geo, R, cutoff), size(geo.vol))
    subdivide > 0 && subdivision_correction_adjoint!(S, geo, R, cutoff, subdivide)
    return S
end

# every cell is a point mass plus quadrupole (box_field), in the east/north/up frame of the station
function point_masses_forward(geo, win, DRHO, rhoTol, cutoff, ::Val{grad}) where {grad}
    (; sx, sy, sz, nx, ny, nz, ex, ey, qx, qy, qz, px, py, pz, vol, m1, m2, m3, mo) = geo
    c2 = cutoff^2
    dg = zeros(eltype(DRHO), length(sx))
    gx = grad ? zero(dg) : nothing
    gy = grad ? zero(dg) : nothing
    for J in CartesianIndices(DRHO)
        abs(DRHO[J]) > rhoTol || continue
        c = G_grav * mGal * vol[J] * DRHO[J]
        for b in station_blocks(win, J)
            @inbounds @simd for i in b
                dx, dy, dz = sx[i] - px[J], sy[i] - py[J], sz[i] - pz[J]
                h = dx * nx[i] + dy * ny[i] + dz * nz[i]        # height of the station above the cell
                t1 = dx * ex[i] + dy * ey[i]                    # east
                t2 = dx * qx[i] + dy * qy[i] + dz * qz[i]       # north
                @fastmath gz, ge, gn = quad_field(mo, J, ifelse(t1 * t1 + t2 * t2 <= c2, c, zero(c)), t1, t2, h, m1[J], m2[J], m3[J])
                dg[i] += gz
                if grad
                    gx[i] += ge
                    gy[i] += gn
                end
            end
        end
    end
    return dg, gx, gy
end

# transpose of point_masses_forward (anomaly only): sensitivity of every cell
function point_masses_adjoint(geo, win, R, cutoff)
    (; sx, sy, sz, nx, ny, nz, ex, ey, qx, qy, qz, px, py, pz, vol, m1, m2, m3, mo) = geo
    c2 = cutoff^2
    S = zeros(size(vol))
    for J in CartesianIndices(S)
        s = 0.0
        for b in station_blocks(win, J)
            @inbounds @simd for i in b
                dx, dy, dz = sx[i] - px[J], sy[i] - py[J], sz[i] - pz[J]
                h = dx * nx[i] + dy * ny[i] + dz * nz[i]
                t1 = dx * ex[i] + dy * ey[i]
                t2 = dx * qx[i] + dy * qy[i] + dz * qz[i]
                @fastmath s += quad_field(mo, J, ifelse(t1 * t1 + t2 * t2 <= c2, R[i], zero(R[i])), t1, t2, h, m1[J], m2[J], m3[J])[1]
            end
        end
        S[J] = G_grav * mGal * vol[J] * s
    end
    return S
end

# ----------------------------------------------------------------------------------------------------------------------
# adaptive subdivision: cells closer to a station than `ratio` times their size are split recursively (in halves along
# every direction that is longer than half the largest one) until each part is far enough or MAX_SUBDIVISION levels are
# reached, similar to Tesseroids (Uieda et al., 2016). Every part is a point mass with its quadrupole term (box_field);
# this is as accurate as 2×2×2 Gauss-Legendre nodes per part, but several times cheaper.

const MAX_SUBDIVISION = 10

function subdivided_field(ell, o, cell, s, n, e, q, ratio, level)
    p, vol, size, ext, m = cell_props(ell, o, cell)
    dx, dy, dz = s[1] - p[1], s[2] - p[2], s[3] - p[3]
    if level >= MAX_SUBDIVISION || dx * dx + dy * dy + dz * dz >= (ratio * size)^2
        # point mass + quadrupole in the frame of the station, as in the main loop (zero for a station at the centre)
        dx * dx + dy * dy + dz * dz > 0 || return (0.0, 0.0, 0.0)
        t1, t2, h = dx * e[1] + dy * e[2] + dz * e[3], dx * q[1] + dy * q[2] + dz * q[3], dx * n[1] + dy * n[2] + dz * n[3]
        m[4] == m[5] == m[6] == 0 && return box_field(vol, t1, t2, h, m[1], m[2], m[3])
        return box_field(vol, t1, t2, h, m...)
    end
    halves(split) = split ? ((0.0, 0.5), (0.5, 1.0)) : ((0.0, 1.0),)
    f = (0.0, 0.0, 0.0)
    for ξ in halves(2ext[1] > size), η in halves(2ext[2] > size), ζ in halves(2ext[3] > size)
        f = f .+ subdivided_field(ell, o, subcell(cell, ξ, η, ζ), s, n, e, q, ratio, level + 1)
    end
    return f
end

# call f(i, J, correction) for every station i that is closer than ratio * size to cell J (and within the cutoff)
function foreach_near_pair(f, geo, tree, cells, cutoff, ratio)
    (; ell, o, sx, sy, sz, nx, ny, nz, ex, ey, qx, qy, qz, sh, px, py, pz, vol, m1, m2, m3, mo, bx, by, bz) = geo
    lowest_station = minimum(sh)
    for J in cells
        cell = cell_box(bx, by, bz, J)
        size = cell_props(ell, o, cell)[3]
        lowest_station - max(maximum(cell[5]), maximum(cell[6])) >= ratio * size && continue      # all stations far above the cell
        for i in inrange(tree, [px[J], py[J], pz[J]], ratio * size)
            s, n = (sx[i], sy[i], sz[i]), (nx[i], ny[i], nz[i])
            e, q = (ex[i], ey[i], 0.0), (qx[i], qy[i], qz[i])
            dx, dy, dz = s[1] - px[J], s[2] - py[J], s[3] - pz[J]
            h = dx * n[1] + dy * n[2] + dz * n[3]
            dx * dx + dy * dy + dz * dz - h * h <= cutoff^2 || continue
            single = quad_field(mo, J, vol[J], dx * e[1] + dy * e[2], dx * q[1] + dy * q[2] + dz * q[3], h, m1[J], m2[J], m3[J])   # as in the main loop
            f(i, J, subdivided_field(ell, o, cell, s, n, e, q, ratio, 0) .- single)
        end
    end
    return nothing
end

# the cells in chunks, for the tasks of a threaded loop
cell_chunks(cells) = collect(Iterators.partition(cells, max(1, cld(length(cells), 4 * Threads.nthreads()))))

# The subdivision correction depends only on the geometry, so it is computed once as a sparse operator: for every
# station i and nearby cell k, the change of (dg, dgdx, dgdy) per unit density. The operators of the last few
# geometries are cached, so that repeated calls (inversions, sampling) do not compute them again.
struct SubdivisionOperator
    station::Vector{Int}
    cell::Vector{Int}
    gz::Vector{Float64}
    gx::Vector{Float64}
    gy::Vector{Float64}
end

const SUBDIVISION_CACHE = Pair{UInt, SubdivisionOperator}[]
const SUBDIVISION_CACHE_SIZE = 4
const SUBDIVISION_LOCK = ReentrantLock()

function subdivision_operator(geo, cutoff, ratio)
    key = hash((geo.sx, geo.sy, geo.sz, geo.nx, geo.ny, geo.nz, geo.ex, geo.ey, geo.qx, geo.qy, geo.qz, geo.sh, geo.bx, geo.by, geo.bz, geo.o, geo.ell, cutoff, ratio))
    return lock(SUBDIVISION_LOCK) do
        k = findfirst(p -> p.first == key, SUBDIVISION_CACHE)
        isnothing(k) || return SUBDIVISION_CACHE[k].second
        op = compute_subdivision_operator(geo, cutoff, ratio)
        length(SUBDIVISION_CACHE) >= SUBDIVISION_CACHE_SIZE && popfirst!(SUBDIVISION_CACHE)
        push!(SUBDIVISION_CACHE, key => op)
        return op
    end
end

function compute_subdivision_operator(geo, cutoff, ratio)
    tree = KDTree(permutedims(hcat(geo.sx, geo.sy, geo.sz)))
    L = LinearIndices(geo.vol)
    tasks = map(cell_chunks(vec(CartesianIndices(geo.vol)))) do chunk
        Threads.@spawn begin
            st, ce, z, x, y = Int[], Int[], Float64[], Float64[], Float64[]
            foreach_near_pair(geo, tree, chunk, cutoff, ratio) do i, J, corr
                push!(st, i)
                push!(ce, L[J])
                push!(z, G_grav * mGal * corr[1])
                push!(x, G_grav * mGal * corr[2])
                push!(y, G_grav * mGal * corr[3])
            end
            (st, ce, z, x, y)
        end
    end
    parts = fetch.(tasks)
    return SubdivisionOperator((reduce(vcat, getindex.(parts, k); init = (k <= 2 ? Int[] : Float64[])) for k in 1:5)...)
end

function subdivision_correction!(dg, gx, gy, geo, DRHO, rhoTol, cutoff, ratio)
    op = subdivision_operator(geo, cutoff, ratio)
    @inbounds for k in eachindex(op.station)
        ρ = DRHO[op.cell[k]]
        abs(ρ) > rhoTol || continue
        i = op.station[k]
        dg[i] += op.gz[k] * ρ
        if !isnothing(gx)
            gx[i] += op.gx[k] * ρ
            gy[i] += op.gy[k] * ρ
        end
    end
    return nothing
end

function subdivision_correction_adjoint!(S, geo, R, cutoff, ratio)
    op = subdivision_operator(geo, cutoff, ratio)
    @inbounds for k in eachindex(op.station)
        S[op.cell[k]] += op.gz[k] * R[op.station[k]]
    end
    return nothing
end

# versions for KernelAbstractions backends, defined in ext/KernelAbstractions_utils.jl
function grav_points_forward_ka end
function grav_points_adjoint_ka end
