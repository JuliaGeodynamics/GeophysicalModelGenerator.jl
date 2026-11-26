# This is voxel_gravity.jl
#
# This function allows the computation of bouguer anomalies over a rectilinear grid
#
# Author: Arne Spang 06/2021

using Printf
using Statistics

export voxel_grav, voxel_grav_adjoint

"""
    voxel_grav(X::Array{Float64, 3}, Y::Array{Float64, 3}, Z::Array{Float64, 3}, RHO::Array{<:Real, 3};
    refMod="AVG", lengthUnit="m", rhoTol=1e-9, Topo=[], outName="Bouguer", printing=true, zRef=0,
    backend=nothing, float_type=Float64, gradients=true, cutoff=Inf, subdivide=2.5)

Computes Bouguer anomalies and gradients

Required arguments:
===
- `X`,`Y`,`Z`:       3D matrices with the coordinates of the grid (X should vary in the first dimension, Y in the second, Z (vertical) in the third). The spacing may vary along each direction.
- `RHO`:         3D matrix with the density at each grid point [kg/m^3]. Any number type works, so that `voxel_grav` can be differentiated with ForwardDiff (with `outName = nothing`).


Optional arguments:
===
- `refMod`:      1D vector with the reference density for each depth. Alternatively, the strings "NE", "SE", "SW", "NW", "AVG" can be used. In that case, one of the corners of `RHO` is used as reference model.In case of "AVG" the reference model is the average of each depth slice.
- `lengthUnit`:  The unit of the coordinates and topography file. Either "m" or "km"
- `rhoTol`:      density differences smaller than `rhoTol` will be ignored [kg/m^3]
- `Topo`:        2D matrix with the topography of the surface, of size `size(X)[1:2]` (only relevant for the paraview output)
- `outName`:     name of the paraview output (do not include file type); `nothing` writes no file (e.g. inside an inversion loop)
- `printing`:    activate printing of additional information [`true` or `false`]
- `zRef`:        height of the stations at which the anomaly is computed (same unit as `Z`): a number, or an array of size `size(X)[1:2]` with the height of every surface point (e.g. the topography, see the free-air and Bouguer examples in the documentation). Only the part of the model below `maximum(zRef)` is taken into account.
- `backend`:     `nothing` (default) computes on one CPU core. A KernelAbstractions.jl backend computes in parallel, e.g. `CPU()` (all threads), `MetalBackend()` (Apple GPU, needs `float_type = Float32`) or `CUDABackend()`. Requires `using KernelAbstractions` (or a GPU package that loads it).
- `float_type`:  number type of the computation with a `backend`; `Float32` is faster on GPUs but less accurate (relative error ≈ 1e-5).
- `gradients`:   compute the horizontal gradients of the anomaly; `false` saves time and returns `nothing` for them.
- `cutoff`:      ignore cells whose horizontal distance from a station is larger than `cutoff` (same unit as the coordinates); e.g. 166.7 km is the standard outer radius of terrain corrections. Default: all cells.
- `subdivide`:   cells closer to a station than `subdivide` times their size are split recursively into smaller point masses, which keeps the result accurate close to the stations; `0` switches this off.

Returns the anomaly `dg` [mGal] and its horizontal gradients `gradX = ∂dg/∂x`, `gradY = ∂dg/∂y` [mGal/m] at the surface grid points (`x`, `y` of the grid, `z = zRef`), computed analytically from the point masses.
`dg` is linear in the density; [`voxel_grav_adjoint`](@ref) gives its adjoint (sensitivity).
"""
function voxel_grav(
        X::Array{Float64, 3}, Y::Array{Float64, 3}, Z::Array{Float64, 3}, RHO::Array{<:Real, 3};
        refMod = "AVG", lengthUnit = "m", rhoTol = 1.0e-9, Topo = [], outName = "Bouguer", printing = true, zRef = 0,
        backend = nothing, float_type = Float64, gradients = true, cutoff = Inf, subdivide = 2.5
    )
    check_backend(backend)
    sz = size(X)
    cut = lowercase(lengthUnit) == "km" ? 1.0e3 * cutoff : float(cutoff)

    ## check input
    X, Y, Z, RHO, RefMod, rhoTol, Topo, orient = checkInput(X, Y, Z, RHO, refMod, lengthUnit, rhoTol, Topo, outName, printing)

    geo = grav_geometry(X, Y, Z, station_heights(zRef, sz, lengthUnit, orient))

    # density anomaly below zRef, averaged from the nodes to the cell centres
    RHO = RHO[:, :, geo.ind] .- reshape(RefMod[geo.ind], 1, 1, :)
    DRHO = cell_average(RHO)

    ############# compute bouguer anomaly #############
    if printing
        @printf "Computing Bouguer anomaly:"
        @time dg, gradX, gradY = isnothing(backend) ? computeBoug(geo, DRHO, rhoTol, cut, Val(gradients)) : grav_forward_ka(backend, float_type, geo, DRHO, rhoTol, cut, gradients)
        @printf "\n"
    else
        dg, gradX, gradY = isnothing(backend) ? computeBoug(geo, DRHO, rhoTol, cut, Val(gradients)) : grav_forward_ka(backend, float_type, geo, DRHO, rhoTol, cut, gradients)
    end
    subdivide > 0 && subdivision_correction!(dg, gradX, gradY, grid_points_geometry(X, Y, Z, geo), DRHO, rhoTol, cut, subdivide)

    ############## write info and output ##############
    if printing
        frac = count(x -> abs(x) > rhoTol, DRHO) / length(DRHO)
        @printf "%.3f %% of the domain contained anomalous densities. If this is more than expected, adjust rhoTol for faster computation.\n" 100 * frac
    end
    if !isnothing(outName)
        coords = cat(X[:, :, 1], Y[:, :, 1], Topo, dims = 4)
        coords = permutedims(coords, [4, 1, 2, 3])
        vtkfile = vtk_grid(outName, coords)
        vtkfile["Bouguer Anomaly [mGal]"] = dg
        if gradients
            vtkfile["Boug_Gradient_X [mGal/m]"] = gradX
            vtkfile["Boug_Gradient_Y [mGal/m]"] = gradY
        end
        outfiles = vtk_save(vtkfile)
        printing && @printf "Wrote output to: %s\n\n" outfiles
    end
    ###################################################

    if orient == 1
        return dg, gradX, gradY
    else
        P(A) = isnothing(A) ? A : permutedims(A, (2, 1))
        return P(dg), P(gradX), P(gradY)
    end
end


"""
    sens = voxel_grav_adjoint(X, Y, Z, R; lengthUnit = "m", zRef = 0, backend = nothing, float_type = Float64, cutoff = Inf, subdivide = 2.5)

Adjoint of [`voxel_grav`](@ref). Returns `sens`, of the same size as `X`, with

    sens[i, j, k] = ∂(sum(R .* dg)) / ∂RHO[i, j, k]      [mGal per kg/m³]

where `dg` is the anomaly computed by `voxel_grav(X, Y, Z, RHO; refMod, zRef)` for a fixed reference-model vector `refMod`
and `R` is a weight for every surface point (same size as `dg`). Since `dg` is linear in `RHO`, `sens` does not depend on `RHO`.

- `R = ones(size(dg))`: sensitivity of the summed anomaly to the density of every grid point.
- `R` = 1 at one station and 0 elsewhere: change of the anomaly at that station, `δg ≈ sum(sens .* δRHO)`.
- `R = dg - dg_obs`: gradient of the misfit `sum((dg - dg_obs) .^ 2) / 2` with respect to `RHO`, for gradient-based inversion.

Grid points above `maximum(zRef)` have zero sensitivity. The cost is that of one call to `voxel_grav` with `rhoTol = 0`; `backend`, `float_type`, `cutoff` and `subdivide` work as there.
"""
function voxel_grav_adjoint(X, Y, Z, R; lengthUnit = "m", zRef = 0, backend = nothing, float_type = Float64, cutoff = Inf, subdivide = 2.5)
    check_backend(backend)
    sz = size(X)
    cut = lowercase(lengthUnit) == "km" ? 1.0e3 * cutoff : float(cutoff)
    size(R) == sz[1:2] || error("R has size $(size(R)) but the grid has $(sz[1:2]) points in the horizontal plane.")
    X, Y, Z, _, _, _, _, orient = checkInput(X, Y, Z, zero(X), zeros(size(X, 3)), lengthUnit, 0.0, [], nothing, false)
    orient == 2 && (R = permutedims(R))

    geo = grav_geometry(X, Y, Z, station_heights(zRef, sz, lengthUnit, orient))
    S = isnothing(backend) ? computeBougAdjoint(geo, R, cut) : grav_adjoint_ka(backend, float_type, geo, R, cut)
    subdivide > 0 && subdivision_correction_adjoint!(S, grid_points_geometry(X, Y, Z, geo), R, cut, subdivide)

    sens = cells_to_nodes(S, size(X), geo.ind)
    return orient == 1 ? sens : permutedims(sens, (2, 1, 3))
end


function checkInput(X, Y, Z, RHO, refMod, lengthUnit, rhoTol, Topo, outName, printing)
    # orientation
    if X[1, 1, 1] ≠ X[2, 1, 1] && X[1, 1, 1] == X[1, 2, 1] && X[1, 1, 1] == X[1, 1, 2]
        orientation = 1
    elseif X[1, 1, 1] == X[2, 1, 1] && X[1, 1, 1] ≠ X[1, 2, 1] && X[1, 1, 1] == X[1, 1, 2]
        orientation = 2
        X = permutedims(X, [2, 1, 3])
        Y = permutedims(Y, [2, 1, 3])
        Z = permutedims(Z, [2, 1, 3])
        RHO = permutedims(RHO, [2, 1, 3])
        isempty(Topo) || (Topo = permutedims(Topo))
    else
        error("Coordinate orientation looks wrong!")
    end

    # dimensions
    if !(size(X) == size(Y) && size(X) == size(Z) && size(X) == size(RHO))
        error("X, Y, Z, RHO must be 3D matrices of the same size.")
    end

    # rectilinear grid: X varies only along the first, Y along the second and Z along the third dimension
    x, y, z = X[:, 1, 1], Y[1, :, 1], Z[1, 1, :]
    mono(v) = all(>(0), diff(v)) || all(<(0), diff(v))
    if !(all(X .== x) && all(Y .== y') && all(Z .== reshape(z, 1, 1, :)) && mono(x) && mono(y) && mono(z))
        error("X, Y, Z must form a rectilinear grid: X may only vary along the first, Y along the second and Z along the third dimension, each monotonically.")
    end

    # Topo
    if !isempty(Topo)
        size(Topo) == size(X)[1:2] || error("Topo has size $(size(Topo)) but the grid has $(size(X)[1:2]) points in the horizontal plane.")
    else
        Topo = zeros(size(X, 1), size(X, 2))
    end

    # tolerance
    if printing
        @printf "Using density tolerance: %.3e\n" rhoTol
    end

    # length unit
    if lengthUnit == "m" || lengthUnit == "M"
        if printing
            @printf "Assuming coordinates to be in meters.\n\n"
        end
    elseif lengthUnit == "km" || lengthUnit == "KM"
        if printing
            @printf "Converting coordinates to meters.\n\n"
        end
        X = X .* 1000
        Y = Y .* 1000
        Z = Z .* 1000
        Topo = Topo .* 1000
    else
        error("lengthUnit should be \"m\" or \"km\".")
    end

    RefMod = reference_model(refMod, RHO)

    # outName
    if !(outName isa Union{String, Nothing})
        error("outName must be a string or nothing.")
    end

    return X, Y, Z, RHO, RefMod, rhoTol, Topo, orientation

end

# reference density for every depth level
function reference_model(refMod, RHO)
    nz = size(RHO, 3)
    if !(typeof(refMod) == String)
        if length(refMod) == nz
            RefMod = refMod
        else
            error("RefMod must have the same length as the third dimension of X, Y, Z, RHO.")
        end
    else
        if refMod == "NE"
            RefMod = RHO[end, end, :]
        elseif refMod == "SE"
            RefMod = RHO[end, 1, :]
        elseif refMod == "SW"
            RefMod = RHO[1, 1, :]
        elseif refMod == "NW"
            RefMod = RHO[1, end, :]
        elseif refMod == "AVG"
            RefMod = vec(mean(RHO, dims = (1, 2)))
        else
            error("RefMod should be NE, SE, SW, NW, AVG or a vector with one value for each depth.")
        end
    end
    return RefMod
end

const G_grav = 6.67408e-11   # gravitational constant [m³/(kg s²)]
const mGal = 1.0e5           # m/s² to mGal

# station heights [m] for every surface point, in the orientation used internally
function station_heights(zRef, sz, lengthUnit, orient)
    zo = zRef isa Number ? fill(float(zRef), sz[1:2]) : float.(zRef)
    size(zo) == sz[1:2] || error("zRef has size $(size(zo)) but the grid has $(sz[1:2]) points in the horizontal plane.")
    orient == 2 && (zo = permutedims(zo))
    return lowercase(lengthUnit) == "km" ? 1000 .* zo : zo
end

# surface points and station heights zo, and centres and volumes of the cells below the highest station
# (grid spacing may vary along x, y and z)
function grav_geometry(X, Y, Z, zo)
    ind = findall(<=(maximum(zo)), Z[1, 1, :])
    x, y, z = X[:, 1, 1], Y[1, :, 1], Z[1, 1, ind]
    mid(v) = (v[1:(end - 1)] .+ v[2:end]) ./ 2
    h(v) = abs.(diff(v))
    vol = h(x) .* h(y)' .* reshape(h(z), 1, 1, :)
    # second moments of the cells (L²/12) along x, y, z
    mx, my, mz = h(x) .^ 2 ./ 12, h(y) .^ 2 ./ 12, h(z) .^ 2 ./ 12
    return (; ind, x, y, xc = mid(x), yc = mid(y), zc = mid(z), vol, zo, mx, my, mz)
end

"""
    gz, ge, gn = box_field(c, t1, t2, h, m1, m2, m3)

Attraction of a homogeneous box with mass factor `c` (`G * mGal * density * volume`) and second moments `m1, m2, m3`
(`L²/12` of its edges along east, north and up), at a station at offset `(t1, t2, h)` (east, north, up; station - box
centre): the point mass plus the quadrupole term of the box, whose error decreases with the fourth power of
size/distance. Returns the downward component and its derivatives towards east and north.
"""
@inline function box_field(c, t1, t2, h, m1, m2, m3)
    # integer constants only, so that the number type of the arguments is kept (Metal has no Float64)
    d2 = t1 * t1 + t2 * t2 + h * h
    q = 1 / d2
    i3 = sqrt(q) * q
    S = m1 + m2 + m3
    T = m1 * t1 * t1 + m2 * t2 * t2 + m3 * h * h
    gz = c * h * i3 * (1 + (15 * T * q - (6 * m3 + 3 * S)) * q / 2)
    i5 = i3 * q
    b = i5 * (-3 + (30 * m3 + 15 * S - 105 * T * q) * q / 2)
    return gz, c * h * t1 * (b + 15 * m1 * i5 * q), c * h * t2 * (b + 15 * m2 * i5 * q)
end

"""
    gz, ge, gn = box_field(c, t1, t2, h, m11, m22, m33, m12, m13, m23)

As `box_field(c, t1, t2, h, m11, m22, m33)`, for a cell with the full tensor of second moments (per unit volume, about
its centroid), e.g. a hexahedron whose top and bottom are inclined. With `m12 = m13 = m23 = 0` it equals the box version.
"""
@inline function box_field(c, t1, t2, h, m11, m22, m33, m12, m13, m23)
    # quadrupole potential c (1/r + (3 r'Mr - r² tr M) / (2 r⁵)); gz = -∂/∂h and its derivatives towards t1, t2
    d2 = t1 * t1 + t2 * t2 + h * h
    q = 1 / d2
    i3 = sqrt(q) * q
    i5 = i3 * q
    S = m11 + m22 + m33
    P1 = m11 * t1 + m12 * t2 + m13 * h
    P2 = m12 * t1 + m22 * t2 + m23 * h
    P3 = m13 * t1 + m23 * t2 + m33 * h
    T = t1 * P1 + t2 * P2 + h * P3
    B = (15 * T * q - 3 * S) * q / 2
    gz = c * (h * i3 * (1 + B) - 3 * P3 * i5)
    K = i5 * (15 * P3 * q - h * (3 + 3 * B) + h * (3 * S - 30 * T * q) * q)
    return gz, c * (t1 * K + i5 * (15 * h * P1 * q - 3 * m13)), c * (t2 * K + i5 * (15 * h * P2 * q - 3 * m23))
end

# the grid of voxel_grav (coordinates in m) as stations and cells for the general point-mass functions
function grid_points_geometry(X, Y, Z, geo)
    x, y = X[:, 1, 1], Y[1, :, 1]
    return points_geometry(x, y, Z[1, 1, :], vec(x .+ 0 .* y'), vec(0 .* x .+ y'), vec(geo.zo), nothing)
end

# mean of the 8 corner nodes of every cell
function cell_average(A)
    n = size(A) .- 1
    C = zeros(eltype(A), n)
    for k in 0:1, j in 0:1, i in 0:1
        C .+= @view A[(1:n[1]) .+ i, (1:n[2]) .+ j, (1:n[3]) .+ k]
    end
    return C ./ 8
end

# transpose of cell_average, for the cells between the depth levels ind of a grid of size sz
function cells_to_nodes(S, sz, ind)
    sens = zeros(sz)
    n = size(S)
    for k in 0:1, j in 0:1, i in 0:1
        @views sens[(1:n[1]) .+ i, (1:n[2]) .+ j, ind[(1:n[3]) .+ k]] .+= S ./ 8
    end
    return sens
end

# stations within the horizontal distance cutoff of every cell centre, along one axis (a range, as the coordinates are monotonic)
function station_window(x, xc, cutoff)
    return map(xc) do c
        i = findall(xi -> abs(xi - c) <= cutoff, x)
        isempty(i) ? (1:0) : (first(i):last(i))
    end
end

# anomaly at the stations, and optionally its horizontal gradients: every cell acts as a point mass at its centre
function computeBoug(geo, DRHO, rhoTol, cutoff, ::Val{grad}) where {grad}
    (; x, y, xc, yc, zc, vol, zo, mx, my, mz) = geo
    dx, dy = x .- xc', y .- yc'            # station - cell
    rx, ry = station_window(x, xc, cutoff), station_window(y, yc, cutoff)
    c2 = cutoff^2
    dg = zeros(eltype(DRHO), length(x), length(y))
    gx = grad ? zero(dg) : nothing
    gy = grad ? zero(dg) : nothing
    for jZ in eachindex(zc), jY in eachindex(yc), jX in eachindex(xc)
        abs(DRHO[jX, jY, jZ]) > rhoTol || continue
        c = G_grav * mGal * vol[jX, jY, jZ] * DRHO[jX, jY, jZ]
        m1, m2, m3 = mx[jX], my[jY], mz[jZ]
        for iY in ry[jY]   # iX innermost: column-major
            dyi = dy[iY, jY]
            @inbounds @simd for iX in rx[jX]
                dxi = dx[iX, jX]
                h = zo[iX, iY] - zc[jZ]        # height of the station above the cell
                @fastmath gz, ge, gn = box_field(ifelse(dxi * dxi + dyi * dyi <= c2, c, zero(c)), dxi, dyi, h, m1, m2, m3)
                dg[iX, iY] += gz
                if grad
                    gx[iX, iY] += ge
                    gy[iX, iY] += gn
                end
            end
        end
    end
    return dg, gx, gy
end

# transpose of computeBoug (anomaly only): sensitivity of every cell
function computeBougAdjoint(geo, R, cutoff)
    (; x, y, xc, yc, zc, vol, zo, mx, my, mz) = geo
    dx, dy = x .- xc', y .- yc'
    rx, ry = station_window(x, xc, cutoff), station_window(y, yc, cutoff)
    c2 = cutoff^2
    S = zeros(length(xc), length(yc), length(zc))
    for jZ in eachindex(zc), jY in eachindex(yc), jX in eachindex(xc)
        s = 0.0
        m1, m2, m3 = mx[jX], my[jY], mz[jZ]
        for iY in ry[jY]
            dyi = dy[iY, jY]
            @inbounds @simd for iX in rx[jX]
                dxi = dx[iX, jX]
                h = zo[iX, iY] - zc[jZ]
                @fastmath s += box_field(ifelse(dxi * dxi + dyi * dyi <= c2, R[iX, iY], zero(R[iX, iY])), dxi, dyi, h, m1, m2, m3)[1]
            end
        end
        S[jX, jY, jZ] = G_grav * mGal * vol[jX, jY, jZ] * s
    end
    return S
end

# parallel versions of computeBoug and computeBougAdjoint for KernelAbstractions backends (CPU threads, GPUs),
# defined in ext/KernelAbstractions_utils.jl
function grav_forward_ka end
function grav_adjoint_ka end

function check_backend(backend)
    isnothing(backend) || hasmethod(grav_forward_ka, Tuple{typeof(backend), Any, Any, Any, Any, Any, Any}) ||
        error("backend = $backend requires KernelAbstractions.jl: add `using KernelAbstractions`.")
    return nothing
end
