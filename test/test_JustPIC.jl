using Test, GeophysicalModelGenerator, JustPIC

# velocity grids for init_particles: each velocity component lives on its own vertex axis
# and on cell centers (with one ghost center per side) along the other axes
function staggered_grids(xvi...)
    centers(v) = collect(range(first(v) - step(v) / 2, last(v) + step(v) / 2, step = step(v)))
    return ntuple(d -> ntuple(i -> i == d ? collect(xvi[i]) : centers(xvi[i]), length(xvi)), length(xvi))
end

active_values(p, A) = vec(A)[vec(p.index.data)]

@testset "2D particles, coordinates in m" begin
    particles = JustPIC.init_particles(JustPIC.CPU, 12, 24, 6, staggered_grids(range(0.0, 100.0e3, length = 21), range(-50.0e3, 0.0, length = 11))...)
    pPhases, pT = JustPIC.init_cell_arrays(particles, Val(2))
    fill!(pPhases.data, 0.0)
    fill!(pT.data, 0.0)
    x, z = (active_values(particles, c.data) ./ 1.0e3 for c in particles.coords)

    pg = particle_grid(particles; length_unit = 1m)
    add_box!(pPhases, pT, pg; xlim = (0, 100), zlim = (-50, 0), phase = LithosphericPhases(Layers = [15, 50], Phases = [1, 2]), T = LinearTemp(Ttop = 0, Tbot = 500))
    add_sphere!(pPhases, nothing, pg; cen = (50, 0, -25), radius = 10, phase = ConstantPhase(5))
    ph, T = active_values(particles, pPhases.data), active_values(particles, pT.data)

    # interfaces are resolved at particle resolution, not at grid resolution
    inside = hypot.(x .- 50, z .+ 25) .<= 10
    @test all(ph[inside] .== 5)
    @test all(ph[.!inside .& (z .> -15)] .== 1)
    @test all(ph[.!inside .& (z .< -15)] .== 2)
    @test T ≈ -z ./ 50 .* 500
    @test all(iszero, vec(pPhases.data)[.!vec(particles.index.data)])

    @test_throws "requires a particle temperature array" add_sphere!(pPhases, nothing, pg; cen = (50, 0, -25), radius = 10, T = ConstantTemp(T = 1000))

    X, Y, _ = xyz_grid(-10:10:110, -10:10:10, 0)
    surf = CartData(X, Y, fill(-5.0, size(X)), (Z = fill(-5.0, size(X)),))
    above, below = above_surface(pg, surf), below_surface(pg, surf)
    @test size(above) == size(particles.index.data)
    @test active_values(particles, above) == (z .> -5)
    @test active_values(particles, below) == (z .< -5)
    @test !any(vec(above)[.!vec(particles.index.data)])

    # free surface from a topography line, in m
    xv = range(0.0, 100.0e3, length = 21)
    chain = JustPIC.init_markerchain(JustPIC.CPU, 6, 3, 12, xv, 0.0)
    Xl, Yl, Zl = xyz_grid(-10:5:110, 0, 0)
    fill_chain_from_vertices!(chain, CartData(Xl, Yl, Xl ./ 100, (Z = Zl,)); length_unit = 1m)
    @test Array(chain.h_vertices) ≈ xv ./ 100
end

@testset "3D particles and free surface" begin
    xv, yv = range(0.0, 10.0, length = 11), range(0.0, 10.0, length = 11)
    particles = JustPIC.init_particles(JustPIC.CPU, 8, 16, 4, staggered_grids(xv, yv, range(-10.0, 0.0, length = 11))...)
    pPhases, = JustPIC.init_cell_arrays(particles, Val(1))
    fill!(pPhases.data, 0.0)
    x, y, z = (active_values(particles, c.data) for c in particles.coords)

    pg = particle_grid(particles; length_unit = 1km)
    add_sphere!(pPhases, nothing, pg; cen = (5, 5, -5), radius = 3, phase = ConstantPhase(1))
    add_box!(pPhases, nothing, pg; xlim = (0, 10), ylim = (0, 4), zlim = (-2, 0), phase = ConstantPhase(2))
    ph = active_values(particles, pPhases.data)

    insphere = (x .- 5) .^ 2 .+ (y .- 5) .^ 2 .+ (z .+ 5) .^ 2 .<= 9
    inbox = (y .<= 4) .& (z .>= -2)
    @test all(ph[insphere] .== 1)
    @test all(ph[inbox] .== 2)
    @test all(ph[.!insphere .& .!inbox] .== 0)

    # the polygon is extruded along y only within ylim, tested per particle
    fill!(pPhases.data, 0.0)
    add_polygon!(pPhases, nothing, pg; xlim = (2, 8, 8, 2), ylim = (0.0, 4.0), zlim = (-8, -8, -2, -2), phase = ConstantPhase(3))
    @test (active_values(particles, pPhases.data) .== 3) == ((2 .<= x .<= 8) .& (y .<= 4) .& (-8 .<= z .<= -2))

    surf = JustPIC.init_marker_surface(JustPIC.CPU, xv, yv, 0.0)
    Xs, Ys, _ = xyz_grid(-1:0.5:11, -1:0.5:11, 0)
    set_topo_from_array!(surf, CartData(Xs, Ys, Xs ./ 10 .- Ys ./ 20, (Z = Xs,)); length_unit = 1km)
    @test Array(surf.topo) ≈ xv ./ 10 .- yv' ./ 20
end
