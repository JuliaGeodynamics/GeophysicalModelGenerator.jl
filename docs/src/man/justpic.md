# JustPIC

[JustPIC.jl](https://github.com/JuliaGeodynamics/JustPIC.jl) stores material properties on particles. Once `JustPIC` is loaded, the geometry routines of `GMG` accept JustPIC particles in place of a grid. Phases and temperatures are then evaluated at the position of every particle. Interfaces are therefore resolved at particle resolution, rather than at the resolution of the numerical grid, and do not show a staircase pattern.

```julia
using GeophysicalModelGenerator, JustPIC

particles = JustPIC.init_particles(backend, nxcell, max_xcell, min_xcell, grid_vel...)
pPhases, pT = JustPIC.init_cell_arrays(particles, Val(2))

pg = particle_grid(particles; length_unit = 1m)     # particle coordinates are in m

add_box!(pPhases, pT, pg; xlim = (0, 500), zlim = (-80, 0),
    phase = LithosphericPhases(Layers = [20, 80], Phases = [1, 2]), T = HalfspaceCoolingTemp(Age = 80))
add_sphere!(pPhases, nothing, pg; cen = (250, 0, -100), radius = 20, phase = ConstantPhase(3))
```

## Particle grid

`particle_grid` is a view of the particle coordinates, and the geometry routines evaluate every particle in place, leaving empty particle slots unchanged. It copies nothing for particles on the CPU; particles on a GPU are processed in host copies.

- **Units:** `length_unit` is the length of one unit of the particle coordinates, for example `1m`, or the characteristic length of a nondimensional model. Geometry parameters are given in km, as everywhere else in `GMG`.
- **2D axes:** in 2D, the second particle coordinate is the vertical one and corresponds to `z` in `GMG`.
- **Temperature:** pass `nothing` instead of a temperature array to set only phases. Passing a thermal structure `T` without a temperature array is an error.

The following routines are supported: `add_box!`, `add_sphere!`, `add_ellipsoid!`, `add_cylinder!`, `add_layer!`, `add_polygon!`, `add_plate!`, `add_slab!`, `add_volcano!`, `add_fault!`, `add_stripes!`, `above_surface` and `below_surface`. The last two return a `Bool` array shaped like the particle `CellArray`s, which is `false` for empty particle slots.

## Free surface

A `GMG` topography sets the free surface of a JustPIC `MarkerChain` (2D) or `MarkerSurface` (3D):

```julia
fill_chain_from_vertices!(chain, topo_line; length_unit = 1m)   # topo_line: CartData of size (nx, 1, 1), e.g. from cross_section
set_topo_from_array!(surf, topo; length_unit = 1m)              # topo: CartData of size (nx, ny, 1)
```

Together with `above_surface(pg, topo)` for the air phase, this keeps the particles and the free surface consistent with the same topography.
