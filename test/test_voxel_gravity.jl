using Test
using GeophysicalModelGenerator
using Statistics

@testset "voxel_grav" begin
    dir = mktempdir()

    @testset "spheres vs. analytical solution" begin
        # survey
        x = [-20.0, 100.0]
        y = [100.0, 200.0]
        z = [-100.0, 10.0]
        nx, ny, nz = 111, 101, 91

        # spheres
        centerX = [50.0, 10.0, 80.0]
        centerY = [150.0, 150.0, 180.0]
        centerZ = [-30.0, -10.0, -50.0]
        radius = [10.0, 3.0, 10.0]
        rho = [2650.0, 2600.0, 2750.0]

        # background density model
        background = [0.0 2700.0]

        G = 6.67408e-11
        mGal = 1.0e5

        x_vec = LinRange(x[1], x[2], nx)
        y_vec = LinRange(y[1], y[2], ny)
        z_vec = LinRange(z[1], z[2], nz)
        Y, X, Z = meshgrid(y_vec, x_vec, z_vec)

        RefMod = zeros(nz)
        RefMod[z_vec .>= 0] .= background[1]
        RefMod[z_vec .< 0] .= background[2]

        RHO = zeros(nx, ny, nz)
        RHO[Z .>= 0] .= background[1]
        RHO[Z .< 0] .= background[2]
        for i in eachindex(rho)
            d = (X .- centerX[i]) .^ 2 .+ (Y .- centerY[i]) .^ 2 .+ (Z .- centerZ[i]) .^ 2
            RHO[d .< radius[i]^2] .= rho[i]
        end

        ana = zeros(nx, ny)
        for iS in eachindex(rho), iX in 1:nx, iY in 1:ny
            d = ((x_vec[iX] - centerX[iS])^2 + (y_vec[iY] - centerY[iS])^2 + (centerZ[iS])^2)^0.5
            depth = -centerZ[iS]
            ana[iX, iY] += (4 * π * G * (radius[iS]^3) * (rho[iS] - background[2]) * depth) / (3 * (d^3)) * mGal
        end

        dg1, = voxel_grav(X, Y, Z, RHO, refMod = RefMod, outName = joinpath(dir, "Benchmark1"), printing = false)
        dg2, = voxel_grav(X, Y, Z, RHO, rhoTol = 1, refMod = "SW", outName = joinpath(dir, "Benchmark2"), printing = false)
        dg3, = voxel_grav(X, Y, Z, RHO, rhoTol = 70, refMod = "NW", outName = joinpath(dir, "Benchmark3"), printing = false)

        maxErr = maximum(abs, ana) / 20
        @test maximum(abs, ana - dg1) < maxErr
        @test maximum(abs, ana - dg2) < maxErr
        @test maximum(abs, ana - dg3) > maxErr   # rhoTol = 70 ignores the density anomalies
    end

    @testset "regular-grid check" begin
        # regular grid in m, whose spacing has rounding errors
        X, Y, Z = xyz_grid(range(-50.0e3, 50.0e3, 32), range(-50.0e3, 50.0e3, 32), range(-30.0e3, 2.0e3, 32))
        RHO = fill(2800.0, size(X))
        RHO[Z .>= 0] .= 0.0
        dg, = voxel_grav(X, Y, Z, RHO, refMod = "NE", outName = joinpath(dir, "Benchmark4"), printing = false)
        @test size(dg) == (32, 32)

        X, Y, Z = xyz_grid([-50.0e3, -40.0e3, 0.0, 50.0e3], range(-50.0e3, 50.0e3, 4), range(-30.0e3, 2.0e3, 4))
        @test_throws "Non-regular grids are not supported yet" voxel_grav(X, Y, Z, fill(2800.0, size(X)), printing = false)
    end

    # small regular grid (coordinates in km) with a single density anomaly
    X, Y, Z = xyz_grid(range(-5.0, 5.0, 6), range(-5.0, 5.0, 6), range(-3.0, 1.0, 5))
    RHO = fill(2800.0, size(X))
    RHO[Z .>= 0] .= 0.0
    RHO[3, 3, 2] = 3000.0
    outName = joinpath(dir, "Bouguer")
    dg_m, gradX_m, gradY_m = voxel_grav(1000 .* X, 1000 .* Y, 1000 .* Z, RHO, refMod = "NE", outName = outName, printing = false)

    @testset "input checks" begin
        @test_throws "Coordinate orientation looks wrong!" voxel_grav(fill(1.0, size(X)), Y, Z, RHO, printing = false)
        @test_throws "X, Y, Z, RHO must be 3D matrices of the same size." voxel_grav(X, Y, Z, RHO[:, :, 1:2], printing = false)
        @test_throws "lengthUnit should be \"m\" or \"km\"." voxel_grav(X, Y, Z, RHO, lengthUnit = "cm", printing = false)
        @test_throws "RefMod must have the same length as the third dimension" voxel_grav(X, Y, Z, RHO, refMod = [2800.0, 2800.0], printing = false)
        @test_throws "RefMod should be NE, SE, SW, NW, AVG" voxel_grav(X, Y, Z, RHO, refMod = "N", printing = false)
        @test_throws "outName must be a string." voxel_grav(X, Y, Z, RHO, refMod = "NE", outName = :Bouguer, printing = false)
    end

    @testset "units and orientation" begin
        @test isfile(outName * ".vts")
        @test size(dg_m) == (6, 6)
        @test maximum(abs, dg_m) > 0

        # coordinates in m and in km give the same anomaly
        dg_km, = voxel_grav(X, Y, Z, RHO, refMod = "NE", lengthUnit = "km", outName = outName, printing = false)
        @test dg_km ≈ dg_m

        # coordinates with x varying along the second dimension give transposed results
        P(A) = permutedims(A, (2, 1, 3))
        dg_p, gradX_p, gradY_p = voxel_grav(1000 .* P(X), 1000 .* P(Y), 1000 .* P(Z), P(RHO), refMod = "NE", outName = outName, printing = false)
        @test dg_p ≈ permutedims(dg_m)
        @test gradX_p ≈ permutedims(gradX_m)
        @test gradY_p ≈ permutedims(gradY_m)
    end

    @testset "reference models" begin
        # a reference-model vector equal to the NE corner column gives the same anomaly
        dg_ref, = voxel_grav(1000 .* X, 1000 .* Y, 1000 .* Z, RHO, refMod = RHO[end, end, :], outName = outName, printing = false)
        @test dg_ref ≈ dg_m

        # the default reference model is the mean density of each depth slice
        dg_avg, = voxel_grav(1000 .* X, 1000 .* Y, 1000 .* Z, RHO, outName = outName, printing = false)
        dg_mean, = voxel_grav(1000 .* X, 1000 .* Y, 1000 .* Z, RHO, refMod = vec(mean(RHO, dims = (1, 2))), outName = outName, printing = false)
        @test size(dg_avg) == (6, 6)
        @test dg_avg ≈ dg_mean
        @test !(dg_avg ≈ dg_m)
    end

    @testset "topography" begin
        # a topography of the wrong size is replaced by a flat one (reported when printing)
        out = mktemp() do path, io
            redirect_stdout(io) do
                voxel_grav(1000 .* X, 1000 .* Y, 1000 .* Z, RHO, refMod = "NE", Topo = zeros(2, 2), outName = outName, printing = true)
            end
            flush(io)
            read(path, String)
        end
        @test occursin("Using a flat topography", out)
        @test occursin("Assuming coordinates to be in meters", out)
    end
end
