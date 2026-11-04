# helper functions to make GMG geometry routines work directly on JustPIC particles and free surfaces
module JustPIC_utils

using JustPIC: Particles, MarkerChain, MarkerSurface
import JustPIC: fill_chain_from_vertices!, set_topo_from_array!

using GeophysicalModelGenerator: CartData, AbstractGeneralGrid, Trench, km, uconvert, ustrip
import GeophysicalModelGenerator: coordinate_grids
using GeophysicalModelGenerator.Interpolations: linear_interpolation
import GeophysicalModelGenerator: particle_grid
import GeophysicalModelGenerator: add_box!, add_sphere!, add_ellipsoid!, add_cylinder!, add_layer!, add_polygon!
import GeophysicalModelGenerator: add_plate!, add_slab!, add_volcano!, add_fault!, add_stripes!
import GeophysicalModelGenerator: above_surface, below_surface

# Particle coordinates in km, read lazily from the particle cell arrays; NaN for empty particle slots, so
# that every shape test of the setup routines is false there and empty slots are left unchanged
struct ParticleCoordinate{T, A, I} <: AbstractArray{T, 3}
    data::A     # coordinate data of the particle cell array
    active::I   # index data of the particle cell array: true for active particles
    scale::T    # km per unit of the particle coordinates
end
Base.size(c::ParticleCoordinate) = size(c.data)
Base.IndexStyle(::Type{<:ParticleCoordinate}) = IndexLinear()
Base.getindex(c::ParticleCoordinate, i::Int) = c.active[i] ? c.scale * c.data[i] : oftype(c.scale, NaN)

# Zero temperatures for calls that only set phases; never written to
struct Zeros{N} <: AbstractArray{Float64, N}
    size::NTuple{N, Int}
end
Base.size(z::Zeros) = z.size
Base.IndexStyle(::Type{<:Zeros}) = IndexLinear()
Base.getindex(z::Zeros, i::Int) = 0.0

struct ParticleGrid{X, Y, Z} <: AbstractGeneralGrid
    X::X
    Y::Y
    Z::Z
end

function coordinate_grids(pg::ParticleGrid; cell = false)
    cell && error("particle grids have no cell centers")
    return pg.X, pg.Y, pg.Z
end

km_per_unit(length_unit) = float(ustrip(uconvert(km, length_unit)))

# The setup routines index arrays point by point, so particle arrays on a GPU are processed in host copies
cpu(A) = A isa Array ? A : Array(A)

function particle_grid(particles::Particles; length_unit)
    s = km_per_unit(length_unit)
    active = cpu(particles.index.data)
    coord(k) = cpu(particles.coords[k].data)
    X = ParticleCoordinate(coord(1), active, s)
    # in 2D, JustPIC's second coordinate is vertical, and y is 0 (zero times the active x coordinate)
    if length(particles.coords) == 2
        return ParticleGrid(X, ParticleCoordinate(coord(1), active, zero(s)), ParticleCoordinate(coord(2), active, s))
    end
    return ParticleGrid(X, ParticleCoordinate(coord(2), active, s), ParticleCoordinate(coord(3), active, s))
end

# routine and the types of its positional arguments after the grid
function_names = (
    (:add_box!, ()), (:add_sphere!, ()), (:add_ellipsoid!, ()), (:add_cylinder!, ()), (:add_layer!, ()), (:add_polygon!, ()),
    (:add_plate!, ()), (:add_slab!, (:Trench,)), (:add_volcano!, ()), (:add_fault!, ()),
)

for (fn, types) in function_names
    args = [Symbol(:arg, i) for i in eachindex(types)]
    signature = join(("pPhases", "pT", "pg", lowercase.(string.(types))...), ", ")
    @eval begin
        """
            $($fn)($($signature); kwargs...)

        Applies `$($fn)` to the particles collected in `pg = particle_grid(particles; length_unit)`, using their own
        coordinates rather than those of a grid, so that interfaces are resolved at particle resolution.
        `pPhases` and `pT` are particle `CellArray`s, as created by `init_cell_arrays`. With `pT = nothing` only phases
        are set, and passing a thermal structure `T` is an error. Geometry parameters are given in km.
        """
        function $fn(pPhases, pT, pg::ParticleGrid, $((:($a::$t) for (a, t) in zip(args, types))...); kwargs...)
            isnothing(pT) && !isnothing(get(kwargs, :T, nothing)) &&
                error("$($fn): a thermal structure `T` requires a particle temperature array `pT`")
            Phase = cpu(pPhases.data)
            Temp = isnothing(pT) ? Zeros(size(Phase)) : cpu(pT.data)
            # without `pT`, `T = nothing` overrides the default thermal structure of the grid method
            kw = isnothing(pT) ? (; kwargs..., T = nothing) : kwargs
            invoke($fn, Tuple{Any, Any, AbstractGeneralGrid, $(types...)}, Phase, Temp, pg, $(args...); kw...)
            Phase === pPhases.data || copyto!(pPhases.data, Phase)
            isnothing(pT) || Temp === pT.data || copyto!(pT.data, Temp)
            return nothing
        end
    end
end

"""
    add_stripes!(pPhases, pg; kwargs...)

Applies `add_stripes!` to the particles collected in `pg = particle_grid(particles; length_unit)`.
"""
function add_stripes!(pPhases, pg::ParticleGrid; kwargs...)
    Phase = cpu(pPhases.data)
    invoke(add_stripes!, Tuple{Any, AbstractGeneralGrid}, Phase, pg; kwargs...)
    Phase === pPhases.data || copyto!(pPhases.data, Phase)
    return nothing
end

for (fn, sel) in ((:above_surface, true), (:below_surface, false))
    @eval begin
        """
            $($fn)(pg, DataSurface_Cart::CartData)

        Returns a `Bool` array, shaped like the data of the particle `CellArray`s, that is `true` for the active particles
        of `pg = particle_grid(particles; length_unit)` $($(sel ? "above" : "below")) the surface `DataSurface_Cart`.
        """
        function $fn(pg::ParticleGrid, DataSurface_Cart::CartData)
            xs, ys, zs = ustrip.(DataSurface_Cart.x.val), ustrip.(DataSurface_Cart.y.val), ustrip.(DataSurface_Cart.z.val)
            surface = linear_interpolation((xs[:, 1, 1], ys[1, :, 1]), zs[:, :, 1])
            compare = $sel ? (>) : (<)
            return ((x, y, z) -> !isnan(z) && compare(z, surface(x, y))).(coordinate_grids(pg)...)
        end
    end
end

"""
    fill_chain_from_vertices!(chain::MarkerChain, topo::CartData; length_unit)

Sets the free surface `chain` to the topography `topo`, a line along `x` of size `(nx, 1, 1)`, for example from
`cross_section`. `length_unit` is the length of one unit of the chain coordinates.
"""
function fill_chain_from_vertices!(chain::MarkerChain, topo::CartData; length_unit)
    size(topo.z.val)[2:3] == (1, 1) || error("`topo` must be a line along x, of size (nx, 1, 1)")
    s = km_per_unit(length_unit)
    h = linear_interpolation(topo.x.val[:, 1, 1], topo.z.val[:, 1, 1]).(s .* cpu(chain.cell_vertices))
    return fill_chain_from_vertices!(chain, h ./ s)
end

"""
    set_topo_from_array!(surf::MarkerSurface, topo::CartData; length_unit)

Sets the free surface `surf` to the topography surface `topo`, of size `(nx, ny, 1)`. `length_unit` is the length of
one unit of the surface coordinates.
"""
function set_topo_from_array!(surf::MarkerSurface, topo::CartData; length_unit)
    size(topo.z.val, 3) == 1 || error("`topo` must be a surface, of size (nx, ny, 1)")
    s = km_per_unit(length_unit)
    itp = linear_interpolation((topo.x.val[:, 1, 1], topo.y.val[1, :, 1]), topo.z.val[:, :, 1])
    return set_topo_from_array!(surf, itp.(s .* cpu(surf.xv), s .* cpu(surf.yv)') ./ s)
end

end # end of module
