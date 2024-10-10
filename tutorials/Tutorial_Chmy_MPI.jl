# # Create an initial model setup for Chmy and run it in parallel
#
# ## Aim
# In this tutorial, your will learn how to use [Chmy](https://github.com/PTsolvers/Chmy.jl) to perform a 2D diffusion simulation
# on one or multiple CPU's or GPU's.
# `Chmy` is a package that allows you to specify grids and fields and create finite difference simulations   
#
# ## 1. Load Chmy and required packages
using Chmy, Chmy.Architectures, Chmy.Grids, Chmy.Fields, Chmy.BoundaryConditions, Chmy.GridOperators, Chmy.KernelLaunch
using KernelAbstractions
using Printf
using CairoMakie
using GeophysicalModelGenerator

# In case you want to use GPU's, you need to sort out whether you have AMD or NVIDIA GPU's
# and load the package accordingly:
#=
 using AMDGPU
 AMDGPU.allowscalar(false)
 using CUDA
 CUDA.allowscalar(false)
=#

# To run this in parallel you need to load this:
using Chmy.Distributed
using MPI
MPI.Init()

# ## 2. Define computational routines
# You need to specify compute kernel for the gradients:
@kernel inbounds = true function compute_q!(q, C, χ, g::StructuredGrid, O)
    I = @index(Global, NTuple)
    I = I + O
    q.x[I...] = -χ * ∂x(C, g, I...)
    q.y[I...] = -χ * ∂y(C, g, I...)
end

# You need to specify compute kernel to update the concentration
@kernel inbounds = true function update_C!(C, q, Δt, g::StructuredGrid, O)
    I = @index(Global, NTuple)
    I = I + O
    C[I...] -= Δt * divg(q, g, I...)
end

# And a main function is required:
@views function main(backend=CPU(); nxy_l=(126, 126))
    arch = Arch(backend, MPI.COMM_WORLD, (0, 0))
    topo = topology(arch)
    me   = global_rank(topo)

    ## geometry
    dims_l = nxy_l
    dims_g = dims_l .* dims(topo)
    grid   = UniformGrid(arch; origin=(-2, -2), extent=(4, 4), dims=dims_g)
    launch = Launcher(arch, grid, outer_width=(16, 8))
    
    ##@info "mpi" me grid

    nx, ny = dims_g
    ## physics
    χ = 1.0
    ## numerics
    Δt = minimum(spacing(grid))^2 / χ / ndims(grid) / 2.1
    ## allocate fields
    C = Field(backend, grid, Center())
    P = Field(backend, grid, Center(), Int32)   # phases
    
    q = VectorField(backend, grid)
    C_v = (me==0) ? KernelAbstractions.zeros(CPU(), Float64, size(interior(C)) .* dims(topo)) : nothing
    
    ## Use the `GeophysicalModelGenerator` to set the initial conditions. Note that 
    ## you have to call this for a `Phases` and a `Temp` grid, which we call `C` here.
    add_box!(P,C,grid,  xlim=(-1.0,1.0), zlim=(-1.0,1.0), phase=ConstantPhase(4), T=ConstantTemp(400))

    ## set BC's and updates the halo:
    bc!(arch, grid, C => Neumann(); exchange=C)
    
    ## visualisation
    fig = Figure(; size=(400, 320))
    ax  = Axis(fig[1, 1]; aspect=DataAspect(), xlabel="x", ylabel="y", title="it = 0")
    plt = heatmap!(ax, centers(grid)..., interior(C) |> Array; colormap=:turbo)
    Colorbar(fig[1, 2], plt)
    ## action
    nt = 100
    for it in 1:nt
        (me==0) && @printf("it = %d/%d \n", it, nt)
        launch(arch, grid, compute_q! => (q, C, χ, grid))
        launch(arch, grid, update_C! => (C, q, Δt, grid); bc=batch(grid, C => Neumann(); exchange=C))
    end
    KernelAbstractions.synchronize(backend)
    gather!(arch, C_v, C)
    if me == 0
        fig = Figure(; size=(400, 320))
        ax  = Axis(fig[1, 1]; aspect=DataAspect(), xlabel="x", ylabel="y", title="it = 0")
        plt = heatmap!(ax, C_v; colormap=:turbo) # how to get the global grid for axes?
        Colorbar(fig[1, 2], plt)
        save("out_gather_$nx.png", fig)
    end
    return
end

# In the code above, the part that calls `GMG` is:

# ```julia 
# add_box!(P,C,grid,  xlim=(-1.0,1.0), zlim=(-1.0,1.0), phase=ConstantPhase(4), T=ConstantTemp(400))
# ```
# which works just like any of the other GMG function.

# ## 3. Run the simulation on one CPU machine or GPU card:

# Running the code on the CPU is done with this:
n = 128
main(; nxy_l=(n, n) .- 2)

# If you instead want to run this on AMD or NVIDIA GPU's do this:
## main(ROCBackend(); nxy_l=(n, n) .- 2)
## main(CUDABackend(); nxy_l=(n, n) .- 2)

# And we need to finalize the simulation with
MPI.Finalize()


# ## 4. Run the simulation on an MPI-parallel machine
# If you want to run this on multiple cores, you will need to setup the [MPI.jl]() package,
# such that `mpiexecjl` is created on the command line.
#
# You can than run it with:
# ```
# mpiexecjl -n 4 --project=. julia Tutorial_Chmy_MPI.jl
# ```

# The full file can be downloaded [here](../../../tutorials/Tutorial_Chmy_MPI.jl)

#src Note: The markdown page is generated using:
#src Literate.markdown("tutorials/Tutorial_Chmy_MPI.jl","docs/src/man",keepcomments=true, execute=false, codefence = "```julia" => "```")
