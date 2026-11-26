# Curvilinear grids whose horizontal grid levels follow surfaces such as the topography or internal interfaces
export deform_grid, curvilinear_grid, sbp_nodes

"""
    Grid_new = deform_grid(Grid::CartData, surfaces; levels = nothing, min_thickness = 0.0)

Deforms the 3D grid `Grid` vertically, such that its horizontal grid levels follow `surfaces`, for example the
topography and internal interfaces. This creates a curvilinear (structured) grid, as required by many seismic wave
propagation codes. The `x` and `y` coordinates and all fields of `Grid` are kept; only `z` changes.

Each surface can be:
- a `CartData` surface (of size `(nx, ny, 1)`), which may have a different resolution than `Grid` but must cover it;
- a number, for a flat surface;
- a function `z(x, y)`, or `z(x)` for a 2D grid (of size `(nx, 1, nz)`).

Surface `i` is placed on the grid level `levels[i]` (the third index of `Grid`). By default, each surface is placed on
the level whose mean height is closest to the mean height of the surface. The bottom and top levels of `Grid` stay
where they are, unless a surface is placed on them. In every column, the points between two such fixed levels are
stretched linearly, so the relative spacing of `Grid` is kept within each layer.

Surfaces must not cross. Where a surface lies above the surface on the next higher level (for example, a basin that
pinches out), it is moved down such that each cell in between has a thickness of at least `min_thickness` (in km).
The heights of the result are strictly monotonic in every column; otherwise (e.g. surfaces that touch with
`min_thickness = 0`) an error is thrown.

# Example
A grid that follows the topography at the top and becomes flat at a depth of 10 km:
```julia
julia> x, y = range(-20, 20, 81), range(-20, 20, 81);
julia> X, Y, _ = xyz_grid(x, y, 0);
julia> H = 2 .* exp.(-(X .^ 2 .+ Y .^ 2) ./ 50);
julia> Topo = CartData(X, Y, H, (Topography = H,));
julia> X, Y, Z = xyz_grid(-20:0.5:20, -20:0.5:20, -30:0.5:0);
julia> Grid = CartData(X, Y, Z, (Z = Z,));
julia> Grid_topo = deform_grid(Grid, (Topo, -10.0))
CartData
    size    : (81, 81, 61)
    x       ϵ [ -20.0 : 20.0]
    y       ϵ [ -20.0 : 20.0]
    z       ϵ [ -30.0 : 2.0]
    fields  : (:Z,)
```
"""
function deform_grid(Grid::CartData, surfaces; levels = nothing, min_thickness = 0.0)
    X, Y, Z = coordinate_grids(Grid)
    surfaces = surfaces isa Tuple || surfaces isa AbstractVector ? surfaces : (surfaces,)
    nz = size(Z, 3)
    nz > 1 || throw(ArgumentError("the grid must be 3D, with more than one level in z"))
    all(>(0), diff(Z; dims = 3)) || throw(ArgumentError("z must increase with the third index of the grid"))
    if isnothing(levels)
        zlevel = vec(mean(Z, dims = (1, 2)))
        levels = [argmin(abs.(zlevel .- mean(surface_heights(S, X, Y)))) for S in surfaces]
    end
    levels = collect(levels)
    length(levels) == length(surfaces) || throw(ArgumentError("levels must give a grid level for each surface"))
    allunique(levels) || throw(ArgumentError("each surface must be placed on a different grid level, got levels $levels; increase the vertical resolution of the grid or pass levels"))
    all(k -> 1 <= k <= nz, levels) || throw(ArgumentError("levels must be between 1 and $nz"))

    # fixed levels from top to bottom, each with its heights
    order = sortperm(levels, rev = true)
    anchors = [levels[order]...]
    heights = [surface_heights(surfaces[i], X, Y) for i in order]
    top = anchors[1] == nz ? nothing : Z[:, :, nz]
    bottom = anchors[end] == 1 ? nothing : Z[:, :, 1]

    Znew = similar(Z, Float64)
    for I in CartesianIndices(view(Z, :, :, 1))
        h = [s[I] for s in heights]
        # a surface that lies above the one above it is moved down
        for i in 2:length(h)
            h[i] = min(h[i], h[i - 1] - min_thickness * (anchors[i - 1] - anchors[i]))
        end
        x, y = X[I, 1], Y[I, 1]
        isnothing(top) || h[1] < top[I] || throw(ArgumentError("a surface lies above the top of the grid at (x, y) = ($x, $y); place it on the top level or extend the grid"))
        isnothing(bottom) || h[end] > bottom[I] || throw(ArgumentError("a surface lies below the bottom of the grid at (x, y) = ($x, $y); place a surface on the bottom level or extend the grid"))
        a = anchors
        isnothing(top) || (a = [nz; a]; h = [top[I]; h])
        isnothing(bottom) || (a = [a; 1]; h = [h; bottom[I]])
        stretch_column!(view(Znew, I, :), view(Z, I, :), a, h)
        all(>(0), diff(view(Znew, I, :))) ||
            throw(ArgumentError("the grid levels touch or cross at (x, y) = ($x, $y), where surfaces touch or cross; pass min_thickness > 0 to keep them apart"))
    end
    return CartData(X, Y, Znew, Grid.fields, Grid.atts)
end

"""
    Grid = curvilinear_grid(x, y, surfaces, nz; spacing = :uniform, min_thickness = 0.0)
    Grid = curvilinear_grid(x, y, surfaces; dz, spacing = :uniform, min_thickness = 0.0)

Creates a curvilinear grid with the columns at `x` and `y` (vectors or ranges, in km; a single `y`, e.g. `0`, gives a
2D grid of size `(nx, 1, nz)`), whose layers are bounded by `surfaces`, given from the top to the bottom. One grid
level lies exactly on every surface. A surface can be a `CartData` surface (with any resolution that covers `x` and
`y`), a number for a flat surface, or a function `z(x, y)` (or `z(x)` for a 2D grid); the bottom is usually a number.
The heights increase with the third index.

Layer `i`, between `surfaces[i]` and `surfaces[i+1]`, has `nz[i]` cells in every column. Alternatively, `dz` [km] gives
the target vertical spacing where a layer is thickest, and the number of cells of each layer follows from it. The
points of a layer are distributed between its bounding surfaces in every column:
- `spacing = :uniform`: with equal spacing;
- `spacing = :sbp`: at the nodes of a summation-by-parts (SBP) grid, see [`sbp_nodes`](@ref); every layer then has
  at least 13 points.

Surfaces must not cross. Where a surface lies above the one above it (for example, a basin that pinches out), it is
moved down such that each cell of the layer in between has a thickness of at least `min_thickness` (in km).

# Examples
A grid with 10 cells between the topography and the base of a sedimentary basin, and 20 cells down to 30 km depth:
```julia
julia> Grid = curvilinear_grid(-20:0.5:20, -20:0.5:20, (Topo, Basin, -30.0), (10, 20); min_thickness = 0.01)
```
A 2D mesh for an SBP finite-difference code (e.g. SeismicWavePC.jl), with the columns at the SBP nodes, a grid line
along the topography and along the Moho, and a flat bottom at 40 km depth:
```julia
julia> topo(x) = 1.5 * exp(-x^2 / 50);  moho(x) = -30 + 3 * tanh(x / 10);
julia> Mesh = curvilinear_grid(sbp_nodes(-50, 50, 401), 0, (topo, moho, -40.0); dz = 0.25, spacing = :sbp)
```
The `Z` field of `Grid` holds the height of each point.
"""
function curvilinear_grid(x, y, surfaces, nz = nothing; dz = nothing, spacing = :uniform, min_thickness = 0.0)
    isnothing(nz) == isnothing(dz) && throw(ArgumentError("give either the number of cells of each layer, nz, or the target spacing dz"))
    spacing in (:uniform, :sbp) || throw(ArgumentError("spacing must be :uniform or :sbp, got $spacing"))
    nlayers = length(surfaces) - 1
    min_cells = spacing == :sbp ? 2 * length(SBP_BOUNDARY_NODES) : 1
    if isnothing(nz)
        Xs, Ys, _ = xyz_grid(x, y, 0)
        H = [surface_heights(S, Xs, Ys) for S in surfaces]
        nz = [max(min_cells, ceil(Int, maximum(H[i] .- H[i + 1]) / dz)) for i in 1:nlayers]
    end
    length(nz) == nlayers || throw(ArgumentError("give one more surface than layers: $(length(surfaces)) surfaces for $(length(nz)) layers"))
    all(>=(min_cells), nz) || throw(ArgumentError("each layer needs at least $min_cells cells, got nz = $nz"))
    # a grid whose heights are the relative positions of the points in their layer (plus the layer number below it),
    # which deform_grid stretches between the surfaces
    η(n) = spacing == :sbp ? sbp_nodes(0.0, 1.0, n + 1) : collect(range(0.0, 1.0, n + 1))
    zref = [0.0]
    for (l, n) in enumerate(reverse(collect(nz)))
        append!(zref, (l - 1) .+ η(n)[2:end])
    end
    X, Y, Z = xyz_grid(x, y, zref)
    levels = reverse(cumsum([1; reverse(collect(nz))]))
    Grid = deform_grid(CartData(X, Y, Z, (Z = Z,)), surfaces; levels, min_thickness)
    Z = NumValue(Grid.z)
    return CartData(X, Y, Z, (Z = Z,))
end

# boundary nodes, in units of the interior spacing, of the SBP operators of SeismicWavePC.jl (dual-pairing upwind
# operators of order 7 of Mattsson, Niemelä & Winters, from SummationByPartsOperators.jl)
const SBP_BOUNDARY_NODES = (0.0, 0.4796891391869102, 1.3753231251827358, 2.3753231251827374, 3.3753231251827396, 4.375323125182741)

"""
    x = sbp_nodes(a, b, n; P = (0.0, 0.4796891391869102, 1.3753231251827358, 2.3753231251827374, 3.3753231251827396, 4.375323125182741))

The `n` nodes on `[a, b]` of a summation-by-parts (SBP) finite-difference grid: equidistant with spacing `h` in the
interior, except the outermost `length(P)` nodes on either side, which lie at `a + h P[j]` and `b - h P[j]`, with
`h = (b - a) / (2P[end] + n - 2length(P) + 1)`. The default `P` are the boundary nodes of the SBP operators of
SeismicWavePC.jl; `n` must be at least `2length(P) + 1` (13). Use them as the `x` (and `y`) of [`curvilinear_grid`](@ref).
"""
function sbp_nodes(a, b, n; P = SBP_BOUNDARY_NODES)
    np = length(P)
    n >= 2np + 1 || throw(ArgumentError("an SBP grid with $np boundary nodes needs at least $(2np + 1) points, got $n"))
    h = (b - a) / (2 * P[end] + n - 2np + 1)
    x = [a + h * (j <= np ? P[j] : P[end] + (j - np)) for j in 1:n]
    for j in 1:np
        x[n + 1 - j] = b - h * P[j]
    end
    return x
end

# Height of a surface at the points (X[:,:,1], Y[:,:,1]) of a grid
surface_heights(S::Number, X, Y) = fill(Float64(S), size(X, 1), size(X, 2))
function surface_heights(S::Function, X, Y)
    x, y = X[:, :, 1], Y[:, :, 1]
    return Float64.(applicable(S, x[1], y[1]) ? S.(x, y) : S.(x))      # z(x, y), or z(x) for a 2D grid
end
function surface_heights(S::CartData, X, Y)
    is_surface(S) || throw(ArgumentError("a CartData surface must have size (nx, ny, 1), got $(size(S))"))
    check_covers(S, finite_extrema(X), finite_extrema(Y), X, Y)
    itp = surface_interpolation(S)
    return itp.(X[:, :, 1], Y[:, :, 1])
end

# Moves the points z of one column to znew, such that the fixed levels `anchors` (from top to bottom, including
# the top and bottom of the column) lie at the heights `h`, and stretches the points in between linearly
function stretch_column!(znew, z, anchors, h)
    for i in 1:(length(anchors) - 1)
        kt, kb = anchors[i], anchors[i + 1]
        for k in (kb + 1):(kt - 1)
            znew[k] = h[i + 1] + (z[k] - z[kb]) / (z[kt] - z[kb]) * (h[i] - h[i + 1])
        end
        # exactly on the surfaces, so that points on a surface are not above or below it by a rounding error
        znew[kt], znew[kb] = h[i], h[i + 1]
    end
    return znew
end

# true if the levels of the 3D array Z are not horizontal, as on a grid created by `deform_grid`
is_curvilinear(Z) = ndims(Z) == 3 && size(Z, 3) > 1 && any(k -> !allequal(view(Z, :, :, k)), axes(Z, 3))

# Fractional vertical index of the point (x, y, z) in a grid with coordinate vectors xv and yv and heights Z, whose
# columns are monotonic. Within a cell, the height varies bilinearly in x and y and linearly with the index, so that
# interpolating in (x, y, index) is the trilinear (isoparametric) interpolation of a hexahedron with vertical edges.
# Beyond the grid, the index is extrapolated linearly (and horizontally from the nearest column), so that the
# interpolation handles points outside the grid as it does on a regular grid.
function vertical_index(xv, yv, Z, x, y, z)
    x, y, z = ustrip(x), ustrip(y), ustrip(z)
    (isnan(x) || isnan(y) || isnan(z)) && return NaN
    (i, s), (j, t) = cell_coordinate(xv, x), cell_coordinate(yv, y)
    i2, j2 = min(i + 1, size(Z, 1)), min(j + 1, size(Z, 2))
    col(k) = (1 - t) * ((1 - s) * Z[i, j, k] + s * Z[i2, j, k]) + t * ((1 - s) * Z[i, j2, k] + s * Z[i2, j2, k])
    up = col(size(Z, 3)) > col(1)
    lo, hi = 1, size(Z, 3)
    while hi - lo > 1
        m = (lo + hi) ÷ 2
        (col(m) <= z) == up ? (lo = m) : (hi = m)
    end
    dz = col(hi) - col(lo)
    return dz == 0 ? float(lo) : lo + (z - col(lo)) / dz
end

# Index of the cell of the sorted vector v that contains x (or the nearest cell), and the position of x within it,
# from 0 to 1
function cell_coordinate(v, x)
    length(v) == 1 && return 1, 0.0
    rev = first(v) > last(v)
    i = clamp(rev ? searchsortedlast(v, x; rev) : searchsortedlast(v, x), 1, length(v) - 1)
    return i, clamp((x - v[i]) / (v[i + 1] - v[i]), 0, 1)
end
