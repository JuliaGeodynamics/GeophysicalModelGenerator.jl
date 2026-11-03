# helper functions to make GMG geometry routines work directly on JustPIC particles and free surfaces
module JustPIC_utils

using JustPIC: Particles, MarkerChain, MarkerSurface
import JustPIC: fill_chain_from_vertices!, set_topo_from_array!

using GeophysicalModelGenerator: CartData, km, uconvert, ustrip
using GeophysicalModelGenerator.Interpolations: linear_interpolation
import GeophysicalModelGenerator: particle_grid
import GeophysicalModelGenerator: add_box!, add_sphere!, add_ellipsoid!, add_cylinder!, add_layer!, add_polygon!
import GeophysicalModelGenerator: add_plate!, add_slab!, add_volcano!, add_fault!, add_stripes!
import GeophysicalModelGenerator: above_surface, below_surface

struct ParticleGrid{G, I, S, D}
    grid::G         # CartData of size (n, 1, 1) holding the n active particles, in km
    active::I       # linear indices of the active particles in the data of the particle CellArrays
    scratch::S      # zero temperatures, for calls that only set phases
    datasize::D     # size of the data of the particle CellArrays
end

km_per_unit(length_unit) = ustrip(uconvert(km, length_unit))

cpu(A) = A isa Array ? A : Array(A)

function particle_grid(particles::Particles; length_unit)
    s = km_per_unit(length_unit)
    active = findall(vec(cpu(particles.index.data)))
    coords = map(c -> reshape(s .* view(vec(cpu(c.data)), active), :, 1, 1), particles.coords)
    X = coords[1]
    # in 2D, JustPIC's second coordinate is vertical
    Y, Z = length(coords) == 2 ? (zero(X), coords[2]) : (coords[2], coords[3])
    return ParticleGrid(CartData(X, Y, Z, (Z = Z,)), active, zero(X), size(particles.index.data))
end

# Calls f with the values of the active particles of the particle CellArray A, shaped like the particle grid.
# CPU data is modified in place; other data is copied to the CPU and written back afterwards.
function with_active_values(f, A, pg::ParticleGrid)
    data = cpu(A.data)
    f(reshape(view(vec(data), pg.active), :, 1, 1))
    data === A.data || copyto!(A.data, data)
    return nothing
end

function_names = (
    :add_box!, :add_sphere!, :add_ellipsoid!, :add_cylinder!, :add_layer!, :add_polygon!,
    :add_plate!, :add_slab!, :add_volcano!, :add_fault!,
)

for fn in function_names
    @eval begin
        """
            $($fn)(pPhases, pT, pg; kwargs...)

        Applies `$($fn)` to the particles collected in `pg = particle_grid(particles; length_unit)`, using their own
        coordinates rather than those of a grid, so that interfaces are resolved at particle resolution.
        `pPhases` and `pT` are particle `CellArray`s, as created by `init_cell_arrays`. With `pT = nothing` only phases
        are set, and passing a thermal structure `T` is an error. Geometry parameters are given in km.
        """
        function $fn(pPhases, pT, pg::ParticleGrid; kwargs...)
            with_active_values(pPhases, pg) do Phase
                if isnothing(pT)
                    isnothing(get(kwargs, :T, nothing)) ||
                        error("$($fn): a thermal structure `T` requires a particle temperature array `pT`")
                    $fn(Phase, pg.scratch, pg.grid; kwargs..., T = nothing)
                else
                    with_active_values(Temp -> $fn(Phase, Temp, pg.grid; kwargs...), pT, pg)
                end
            end
            return nothing
        end
    end
end

"""
    add_stripes!(pPhases, pg; kwargs...)

Applies `add_stripes!` to the particles collected in `pg = particle_grid(particles; length_unit)`.
"""
function add_stripes!(pPhases, pg::ParticleGrid; kwargs...)
    with_active_values(Phase -> add_stripes!(Phase, pg.grid; kwargs...), pPhases, pg)
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
            mask = falses(pg.datasize)
            vec(mask)[pg.active] = vec(above_surface(pg.grid, DataSurface_Cart; above = $sel))
            return mask
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
