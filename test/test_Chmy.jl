using Test, GeophysicalModelGenerator

using Chmy, Chmy.Architectures, Chmy.Grids, Chmy.Fields
using KernelAbstractions # for backend-agnostic kernels

backend = CPU()
arch = Arch(backend)

# 3D test
lx, ly, lz  = 10.0, 11.0, 12.0
nx, ny, nz  = 10,11,12
grid   = UniformGrid(arch;
                    origin=(-lx/2, -ly/2, -lz/2),
                    extent=(lx, ly, lz),
                    dims=(nx, ny, nz))

# create field
Temp_C   =  Field(backend, grid, Center(), Float64; halo=1)
Phases_C =  Field(backend, grid, Center(), Int32; halo=1)
Temp_V   =  Field(backend, grid, Vertex(), Float64; halo=1)
Phases_V =  Field(backend, grid, Vertex(), Int32; halo=1)

# grid:
CartGrid =  create_CartGrid(grid)

@test sum.(CartGrid.coord1D) == (0.0, 0.0, 0.0)

# test add_box! directly. Note that this requires you to specify a "cell" keyword for Center() locations 
add_box!(Phases_C,Temp_C,CartGrid, xlim=(0,1.0), zlim=(-2,0), phase=ConstantPhase(3), cell=true)
@test extrema(Phases_C) == (0,3)

add_box!(Phases_V,Temp_V,CartGrid, xlim=(0,1.0), zlim=(-2,0), phase=ConstantPhase(3))
@test extrema(Phases_V) == (0,3)

# multiple dispatch functions
add_box!(Phases_C,Temp_C,grid, xlim=(0,1.0), zlim=(-2,0), phase=ConstantPhase(2))
@test extrema(Phases_C) == (0,2)

add_box!(Phases_V,Temp_V,grid, xlim=(0,1.0), zlim=(-2,0), phase=ConstantPhase(2))
@test extrema(Phases_V) == (0,2)

add_sphere!(Phases_V,Temp_V,grid,  cen=(0,0,-1), radius=2.5, phase=ConstantPhase(3), T=ConstantTemp(800))
@test extrema(Phases_V) == (0,3)
@test extrema(Temp_V) == (0.0,800.0)

# test above/below surface intersection
Topo_cart   =   CartData(xyz_grid(-6:.2:6,-12:.2:13,0));
ind         =   above_surface(grid, Phases_V, Topo_cart);
Phases_V[ind] .= 4;
@test extrema(Phases_V) == (0,3)

ind         =   above_surface(grid, Phases_C, Topo_cart);
Phases_C[ind] .= 4;
@test extrema(Phases_V) == (0,4)


