using Test, GeophysicalModelGenerator
using GeophysicalModelGenerator: surface_heights, is_curvilinear

# Surfaces on a coarser grid than the 3D grid, and with a different nx, ny
xs, ys = range(-25, 25, 41), range(-25, 25, 37)
Xs, Ys, _ = xyz_grid(xs, ys, 0)
topo(x, y) = 1.5 * exp(-(x^2 + y^2) / 60) + 0.1 * x / 25
basin(x, y) = -3 + 4 * exp(-((x - 5)^2 + y^2) / 100)    # crosses the topography near (5, 0)
Topo = CartData(Xs, Ys, topo.(Xs, Ys), (Topography = topo.(Xs, Ys),))
Basin = CartData(Xs, Ys, basin.(Xs, Ys), (Z = basin.(Xs, Ys),))

X0, Y0, Z0 = xyz_grid(-20:1:20, -15:1:15, -30:0.5:0)
Grid0 = CartData(X0, Y0, Z0, (Z = Z0,))
X2, Y2 = X0[:, :, 1], Y0[:, :, 1]
k10 = findfirst(==(-10.0), Z0[1, 1, :])

# a field that is linear in x, y and z is interpolated exactly by trilinear (isoparametric) interpolation
f(x, y, z) = 6 + 0.01x - 0.02y - 0.05z

@testset "deform_grid" begin
    G = deform_grid(Grid0, (Topo, -10.0))
    X, Y, Z = coordinate_grids(G)
    @test X == X0 && Y == Y0 && G.fields.Z == Z0           # x, y and fields are kept
    @test Z[:, :, end] == surface_heights(Topo, X0, Y0)     # exactly on the surface
    @test all(==(-10.0), Z[:, :, k10])
    @test Z[:, :, 1:k10] == Z0[:, :, 1:k10]                 # undeformed below the flat surface
    @test all(diff(Z, dims = 3) .> 0)
    @test is_curvilinear(Z) && !is_curvilinear(Z0)
    # equal spacing between the topography and the flat level in each column
    dz = diff(Z[:, :, k10:end], dims = 3)
    @test maximum(abs.(dz .- dz[:, :, 1:1])) < 1.0e-12

    # a surface given as a function, and default levels
    @test deform_grid(Grid0, topo).z.val[:, :, end] ≈ topo.(X2, Y2)
    Gd = deform_grid(Grid0, (Topo, -10.0, -20.0))
    @test all(==(-20.0), Gd.z.val[:, :, findfirst(==(-20.0), Z0[1, 1, :])])

    # crossing surfaces: the basin is moved below the topography, with a minimum cell thickness
    Gb = deform_grid(Grid0, (Topo, Basin); levels = (61, 55), min_thickness = 0.05)
    @test minimum(diff(Gb.z.val[:, :, 55:61], dims = 3)) >= 0.05 - 1.0e-12
    @test Gb.z.val[:, :, 55] == min.(surface_heights(Basin, X0, Y0), Gb.z.val[:, :, 61] .- 6 * 0.05)

    # errors
    @test_throws ArgumentError deform_grid(Grid0, (Topo, -10.0); levels = (61, 61))
    @test_throws ArgumentError deform_grid(Grid0, Topo; levels = (55,))        # topography above the top of the grid
    @test_throws ArgumentError deform_grid(Grid0, -40.0; levels = (5,))        # below the bottom
    @test_throws "z must increase" deform_grid(CartData(X0, Y0, reverse(Z0, dims = 3), (Z = Z0,)), -10.0)
    Small = CartData(xyz_grid(-5:1:5, -5:1:5, 0)..., (Z = zeros(11, 11, 1),))
    @test_throws ArgumentError deform_grid(Grid0, Small)                       # does not cover the grid
end

@testset "curvilinear_grid" begin
    C = curvilinear_grid(-20:1:20, -15:1:15, (Topo, Basin, -30.0), (6, 20); min_thickness = 0.01)
    @test size(C) == (41, 31, 27)
    @test keys(C.fields) == (:Z,)
    @test C.z.val[:, :, end] == surface_heights(Topo, X0, Y0)
    @test C.z.val[:, :, 21] == min.(surface_heights(Basin, X0, Y0), C.z.val[:, :, end] .- 6 * 0.01)
    @test all(==(-30.0), C.z.val[:, :, 1])
    @test minimum(diff(C.z.val, dims = 3)) >= 0.01 - 1.0e-12
    # equal spacing within each layer
    for r in (1:21, 21:27)
        dz = diff(C.z.val[:, :, r], dims = 3)
        @test maximum(abs.(dz .- dz[:, :, 1:1])) < 1.0e-12
    end
    @test_throws ArgumentError curvilinear_grid(-20:1:20, -15:1:15, (Topo, -30.0), (6, 20))
end

@testset "given (non-uniform) node coordinates" begin
    # nodes of a summation-by-parts (SBP) grid: equidistant, except the outer six columns on each side
    h = 0.5
    edge = h .* [0.48, 0.9, 1.04, 0.97, 1.02, 0.99]
    interior = 20 - 2 * sum(edge)
    x = cumsum([-10.0; edge; fill(interior / round(interior / h), Int(round(interior / h))); reverse(edge)])
    @test x[end] ≈ 10.0
    y = range(-8, 8, 33)
    C = curvilinear_grid(x, y, (Topo, Basin, -30.0), (6, 20); min_thickness = 0.01)
    X, Y, Z = coordinate_grids(C)
    @test X[:, 1, 1] == x && all(X .== x) && all(Y .== y')              # the columns are exactly at the given nodes
    @test Z[:, :, end] == surface_heights(Topo, X, Y)
    # deform_grid of a grid on the same nodes, and data and gravity on it
    G = deform_grid(CartData(xyz_grid(x, y, -30:0.5:0)..., (Z = xyz_grid(x, y, -30:0.5:0)[3],)), (Topo, -10.0))
    @test G.x.val[:, 1, 1] == x && G.z.val[:, :, end] == surface_heights(Topo, X[:, :, 1:1], Y[:, :, 1:1])
    Gf = addfield(C, "F", f.(X, Y, Z))
    Xr, Yr, Zr = xyz_grid(-9.7:0.9:9.7, -7.5:1.1:7.5, -20:0.7:-1)
    @test maximum(abs.(filter(!isnan, interpolate_datafields(Gf, Xr, Yr, Zr).fields.F .- f.(Xr, Yr, Zr)))) < 1.0e-12
    dg, = voxel_grav(C, fill(2670.0, size(X)); refMod = zeros(size(X, 3)))
    @test size(dg) == (length(x), length(y)) && all(isfinite, dg)
end

@testset "2D meshes for SBP finite-difference codes" begin
    # column positions: equidistant except the outermost 6 on each side
    P = (0.0, 0.4796891391869102, 1.3753231251827358, 2.3753231251827374, 3.3753231251827396, 4.375323125182741)
    a, b, N = -50.0, 50.0, 201
    x = sbp_nodes(a, b, N)
    h = (b - a) / (2P[6] + N - 11)
    @test all(x[j] ≈ a + h * P[j] for j in 1:6)
    @test all(x[j] ≈ a + h * (P[6] + (j - 6)) for j in 7:(N - 6))
    @test all(x[N + 1 - j] ≈ b - h * P[j] for j in 1:6)
    @test x[1] == a && x[end] == b
    @test_throws ArgumentError sbp_nodes(0, 1, 12)

    # topography, a sediment base (as a coarser CartData surface) and a Moho; flat bottom
    topo(x) = 1.5 * exp(-x^2 / 50) + 0.4 * sin(x / 3)
    moho(x) = -30 + 3 * tanh(x / 10)
    xs = range(-60, 60, 97)
    sed = CartData(reshape(collect(xs), :, 1, 1), zeros(97, 1, 1), reshape(@.(-2 + 1.2 * exp(-(xs - 10)^2 / 80)), :, 1, 1), (z = zeros(97, 1, 1),))
    M = curvilinear_grid(x, 0, (topo, sed, moho, -40.0); dz = 0.5, spacing = :sbp)
    X, Y, Z = coordinate_grids(M)
    nz = size(Z, 3)
    @test size(M) == (N, 1, nz) && all(==(0), Y)
    @test all(k -> X[:, 1, k] == x, 1:nz)                            # vertical columns at the SBP nodes
    @test all(diff(Z, dims = 3) .> 0)                                 # z strictly increasing in every column
    # one row exactly on every surface
    zsed = GeophysicalModelGenerator.surface_interpolation(sed).(x, 0.0)
    ksed = findfirst(k -> Z[:, 1, k] == zsed, 1:nz)
    kmoho = findfirst(k -> Z[:, 1, k] == moho.(x), 1:nz)
    @test Z[:, 1, end] == topo.(x) && all(==(-40.0), Z[:, 1, 1])
    @test !isnothing(ksed) && !isnothing(kmoho)
    # every block has at least 13 rows, distributed like the SBP nodes, and about dz spacing where it is thickest
    for (k1, k2) in ((1, kmoho), (kmoho, ksed), (ksed, nz))
        @test k2 - k1 + 1 >= 13
        i = argmax(Z[:, 1, k2] .- Z[:, 1, k1])
        @test Z[i, 1, k1:k2] ≈ sbp_nodes(Z[i, 1, k1], Z[i, 1, k2], k2 - k1 + 1)
    end
    @test maximum(diff(Z[:, 1, kmoho:ksed], dims = 2)) <= 1.1 * 0.5            # the interior SBP spacing is slightly larger than thickness / (n - 1)

    # no folding under a steep cliff; surfaces that touch are an error
    cliff(x) = 2 + 2 * tanh(4 * x)                     # slope up to 8
    @test all(diff(curvilinear_grid(x, 0, (cliff, -1.0, -10.0); dz = 0.2, spacing = :sbp).z.val, dims = 3) .> 0)
    @test_throws ArgumentError curvilinear_grid(x, 0, (x -> max(0.0, 0.2x), 0.0, -10.0); dz = 0.5)
    @test_throws ArgumentError curvilinear_grid(x, 0, (topo, -40.0), (5,); spacing = :sbp)   # fewer than 13 rows
end

@testset "data processing on curvilinear grids" begin
    G = deform_grid(Grid0, (Topo, -10.0))
    X, Y, Z = coordinate_grids(G)
    G = addfield(G, "F", f.(X, Y, Z))

    # onto a curvilinear grid, from a regular one
    V = CartData(xyz_grid(-25:1:25, -25:1:25, -40:1:5)..., (Z = zeros(51, 51, 46),))
    V = addfield(V, "F", f.(coordinate_grids(V)...))
    @test interpolate_datafields(V, X, Y, Z).fields.F ≈ G.fields.F

    # from a curvilinear grid; points above the topography are NaN
    Xr, Yr, Zr = xyz_grid(-15:1.3:15, -12:1.1:12, -12:0.37:2)
    F = interpolate_datafields(G, Xr, Yr, Zr).fields.F
    inside = Zr .< topo.(Xr, Yr) .- 0.05
    @test maximum(abs.(F[inside] .- f.(Xr[inside], Yr[inside], Zr[inside]))) < 1.0e-12
    @test all(isnan, F[Zr .> topo.(Xr, Yr) .+ 0.05])

    cs = cross_section(G; Depth_level = -1.0)              # interpolated, as the grid levels are not horizontal
    @test all(==(-1.0), ustrip.(cs.z.val))
    @test maximum(abs.(filter(!isnan, cs.fields.F .- f.(ustrip.(cs.x.val), ustrip.(cs.y.val), -1.0)))) < 1.0e-12
    cs = cross_section(G; Start = (-10, -5), End = (10, 5))
    @test maximum(abs.(filter(!isnan, cs.fields.F .- f.(ustrip.(cs.x.val), ustrip.(cs.y.val), ustrip.(cs.z.val))))) < 1.0e-12
    cs = cross_section(G; Lat_level = 0.0)                 # a row of the grid, which follows the topography
    @test ustrip.(cs.z.val[:, 1, end]) == Z[:, 16, end]

    s = extract_subvolume(G; X_level = (-5, 5), Y_level = (-5, 5), Z_level = (-5, 0.5))
    @test maximum(abs.(filter(!isnan, s.fields.F .- f.(ustrip.(s.x.val), ustrip.(s.y.val), ustrip.(s.z.val))))) < 1.0e-12
    @test_throws ArgumentError extract_subvolume(G; Interpolate = false, Z_level = (-5, 0))
    @test size(extract_subvolume(G; Interpolate = false, X_level = (-5, 5))) == (11, 31, 61)

    Xs2, Ys2, _ = xyz_grid(-15:1:15, -10:1:10, 0)
    S = CartData(Xs2, Ys2, fill(-1.0, size(Xs2)), (Z = Xs2,))
    @test interpolate_data_surface(G, S).fields.F ≈ f.(Xs2, Ys2, -1.0)

    # integer fields (phases) take the nearest point
    Gp = addfield(G, "Phase", Int64.(Z .< -5))
    P = interpolate_datafields(Gp, Xr, Yr, Zr).fields.Phase
    @test P isa Array{Int64} && P[Zr .< -5.3] == ones(Int64, count(Zr .< -5.3))
end

@testset "model setup on curvilinear grids" begin
    G = deform_grid(Grid0, (Topo, -10.0))
    Phase, Temp = zeros(Int64, size(G)), fill(1350.0, size(G))
    add_layer!(Phase, Temp, G; zlim = (-10.0, Topo), phase = ConstantPhase(1), T = LinearTemp(Ttop = 0, Tbot = 300))
    @test all(Phase[:, :, k10:end] .== 1) && all(Phase[:, :, 1:(k10 - 1)] .== 0)
    @test all(==(0.0), Temp[:, :, end]) && all(==(300.0), Temp[:, :, k10])
    @test count(above_surface(G, Topo)) == 0
    @test count(below_surface(G, Topo)) == length(Phase) - length(X2)

    # the same layer on cell centers
    Pc, Tc = zeros(Int64, size(G) .- 1), zeros(size(G) .- 1)
    add_layer!(Pc, Tc, G; zlim = (-10.0, Topo), phase = ConstantPhase(1), cell = true)
    @test all(Pc[:, :, k10:end] .== 1) && all(Pc[:, :, 1:(k10 - 1)] .== 0)

    # Paraview output of a curvilinear grid
    G = addfield(G, (; Phase, Temp))
    write_paraview(G, joinpath(@__DIR__, "curvilinear_test"); verbose = false)
    @test isfile(joinpath(@__DIR__, "curvilinear_test.vts"))
end
