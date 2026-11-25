using Test
using GeophysicalModelGenerator
using Statistics
using StableRNGs
using ForwardDiff
using KernelAbstractions

@testset "voxel_grav" begin
    dir = mktempdir()

    # regular grid, and a grid whose spacing varies by up to ±2a in every direction
    wiggle(r, a = 0.2) = [v + a * step(r) * sin(3i) * (1 < i < length(r)) for (i, v) in enumerate(r)]
    @testset "spheres vs. analytical solution, $grid grid" for grid in ("regular", "variable")
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
        if grid == "variable"
            # ±10%: the small sphere is only ~5 cells wide, so larger variations change its discretized volume too much
            x_vec, y_vec, z_vec = wiggle(x_vec, 0.05), wiggle(y_vec, 0.05), wiggle(z_vec, 0.05)
        end
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

        dg1, gx1, gy1 = voxel_grav(X, Y, Z, RHO, refMod = RefMod, outName = joinpath(dir, "Benchmark1"), printing = false)
        dg2, = voxel_grav(X, Y, Z, RHO, rhoTol = 1, refMod = "SW", outName = joinpath(dir, "Benchmark2"), printing = false)
        dg3, = voxel_grav(X, Y, Z, RHO, rhoTol = 70, refMod = "NW", outName = joinpath(dir, "Benchmark3"), printing = false)

        maxErr = maximum(abs, ana) / 20
        @test maximum(abs, ana - dg1) < maxErr

        # horizontal gradients [mGal/m]
        anagx, anagy = zeros(nx, ny), zeros(nx, ny)
        for iS in eachindex(rho), iX in 1:nx, iY in 1:ny
            dx, dy, h = x_vec[iX] - centerX[iS], y_vec[iY] - centerY[iS], -centerZ[iS]
            t = -3 * 4 / 3 * π * G * radius[iS]^3 * (rho[iS] - background[2]) * h / (dx^2 + dy^2 + h^2)^2.5 * mGal
            anagx[iX, iY] += t * dx
            anagy[iX, iY] += t * dy
        end
        @test maximum(abs, anagx - gx1) < maximum(abs, anagx) / 10
        @test maximum(abs, anagy - gy1) < maximum(abs, anagy) / 10
        @test maximum(abs, ana - dg2) < maxErr
        @test maximum(abs, ana - dg3) > maxErr   # rhoTol = 70 ignores the density anomalies
    end

    @testset "grid check" begin
        # regular grid in m, whose spacing has rounding errors
        X, Y, Z = xyz_grid(range(-50.0e3, 50.0e3, 32), range(-50.0e3, 50.0e3, 32), range(-30.0e3, 2.0e3, 32))
        RHO = fill(2800.0, size(X))
        RHO[Z .>= 0] .= 0.0
        dg, = voxel_grav(X, Y, Z, RHO, refMod = "NE", outName = joinpath(dir, "Benchmark4"), printing = false)
        @test size(dg) == (32, 32)

        X, Y, Z = xyz_grid([-50.0e3, -40.0e3, 0.0, 50.0e3], range(-50.0e3, 50.0e3, 4), range(-30.0e3, 2.0e3, 4))
        @test size(voxel_grav(X, Y, Z, fill(2800.0, size(X)), outName = nothing, printing = false)[1]) == (4, 4)
        X[2, 3, 1] += 1   # not a rectilinear grid
        @test_throws "X, Y, Z must form a rectilinear grid" voxel_grav(X, Y, Z, fill(2800.0, size(X)), printing = false)
        X, Y, Z = xyz_grid([-50.0e3, -40.0e3, -45.0e3, 50.0e3], range(-50.0e3, 50.0e3, 4), range(-30.0e3, 2.0e3, 4))
        @test_throws "X, Y, Z must form a rectilinear grid" voxel_grav(X, Y, Z, fill(2800.0, size(X)), printing = false)
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
        @test_throws "outName must be a string or nothing." voxel_grav(X, Y, Z, RHO, refMod = "NE", outName = :Bouguer, printing = false)
        @test_throws "R has size (2, 2)" voxel_grav_adjoint(X, Y, Z, zeros(2, 2))
    end

    @testset "no output file" begin
        name = joinpath(dir, "NoFile")
        dg, = voxel_grav(1000 .* X, 1000 .* Y, 1000 .* Z, RHO, refMod = "NE", outName = nothing, printing = false)
        @test dg ≈ dg_m
        @test !isfile(name * ".vts")
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

    @testset "cells above sea level" begin
        # the anomaly does not depend on the order of the depth levels, some of which lie above sea level
        Xa, Ya, Za = xyz_grid(range(-5.0e3, 5.0e3, 6), range(-5.0e3, 5.0e3, 6), range(-20.0e3, 4.0e3, 7))
        RHOa = fill(2800.0, size(Xa))
        RHOa[Za .>= 0] .= 0.0
        RHOa[3, 3, 3] = 3000.0
        dg_up, = voxel_grav(Xa, Ya, Za, RHOa, refMod = "NE", outName = outName, printing = false)
        R(A) = A[:, :, end:-1:1]
        dg_down, = voxel_grav(Xa, Ya, R(Za), R(RHOa), refMod = "NE", outName = outName, printing = false)
        @test dg_up ≈ dg_down
    end

    @testset "gradients are per meter" begin
        # interior gradient equals the finite difference of dg over the 2 km spacing
        dx = 2000.0
        @test gradX_m[2, 3] ≈ (dg_m[3, 3] - dg_m[1, 3]) / (2dx) rtol = 0.5
        @test gradY_m[3, 2] ≈ (dg_m[3, 3] - dg_m[3, 1]) / (2dx) rtol = 0.5

        @test voxel_grav(1000 .* X, 1000 .* Y, 1000 .* Z, RHO, refMod = "NE", outName = nothing, printing = false, gradients = false)[2] === nothing
    end

    @testset "cutoff" begin
        Xv, Yv, Zv = xyz_grid(wiggle(range(-20.0e3, 20.0e3, 21)), wiggle(range(-15.0e3, 15.0e3, 16)), range(-6.0e3, 0.0, 7))
        m = rand(StableRNG(9), size(Xv)...)
        r = rand(StableRNG(10), 21, 16)
        g(; kw...) = voxel_grav(Xv, Yv, Zv, m; refMod = zeros(7), rhoTol = 0, outName = nothing, printing = false, kw...)
        @test all(g(cutoff = 1.0e6) .≈ g())
        # cells farther than the cutoff are left out: compare with masking the density of a single station's view
        dg = g(cutoff = 8.0e3)[1]
        far = (Xv .- Xv[11, 8, 1]) .^ 2 .+ (Yv .- Yv[11, 8, 1]) .^ 2 .> 12.0e3^2
        near = (Xv .- Xv[11, 8, 1]) .^ 2 .+ (Yv .- Yv[11, 8, 1]) .^ 2 .< 4.0e3^2
        @test g(cutoff = 8.0e3)[1][11, 8] ≈ voxel_grav(Xv, Yv, Zv, m .* .!far; refMod = zeros(7), rhoTol = 0, outName = nothing, printing = false, cutoff = 8.0e3)[1][11, 8]
        @test !(dg[11, 8] ≈ voxel_grav(Xv, Yv, Zv, m .* near; refMod = zeros(7), rhoTol = 0, outName = nothing, printing = false)[1][11, 8])
        # the adjoint uses the same cutoff
        sens = voxel_grav_adjoint(Xv, Yv, Zv, r; cutoff = 8.0e3)
        @test sum(r .* dg) ≈ sum(sens .* m) rtol = 1.0e-10
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
        @test_throws "Topo has size (2, 2)" voxel_grav(1000 .* X, 1000 .* Y, 1000 .* Z, RHO, refMod = "NE", Topo = zeros(2, 2), outName = outName)

        # orientation 2 with nx != ny: Topo follows the grid permutation
        Xo = [Float64(j) for i in 1:4, j in 1:3, k in 1:4]
        Yo = [Float64(i) for i in 1:4, j in 1:3, k in 1:4]
        Zo = [Float64(k) for i in 1:4, j in 1:3, k in 1:4]
        chk = GeophysicalModelGenerator.checkInput(Xo, Yo, Zo, ones(size(Xo)), "NE", "m", 1.0e-6, fill(500.0, 4, 3), outName, false)
        @test chk[7] == fill(500.0, 3, 4)
    end

    # tests below are ported from other gravity codes:
    # Harmonica (test_prism_against_infinite_slab), Choclo (test_mirrored_points, test_sign) and
    # SimPEG (test_Jtvec, test_G_t_dot_v: forward and adjoint operators must be consistent)
    G = 6.67408e-11
    fwd(X, Y, Z, RHO) = voxel_grav(X, Y, Z, RHO, refMod = zeros(size(X, 3)), rhoTol = 0, outName = nothing, printing = false)[1]

    @testset "large prism vs. exact prism and infinite slab" begin
        # gz [mGal] of a prism below the origin, z positive downwards (Nagy et al. 2000, J. Geodesy)
        function prism(x1, x2, y1, y2, z1, z2, ρ)
            s = 0.0
            for (i, x) in enumerate((x1, x2)), (j, y) in enumerate((y1, y2)), (k, z) in enumerate((z1, z2))
                r = sqrt(x^2 + y^2 + z^2)
                s += (-1)^(i + j + k) * (x * log(y + r) + y * log(x + r) - z * atan(x * y / (z * r)))
            end
            return -G * ρ * s * 1.0e5
        end
        h, ρ = 1000.0, 300.0
        slab = 2π * G * ρ * h * 1.0e5              # Bouguer slab
        err = Float64[]
        for n in (9, 17, 33)                        # 100 m cells, growing horizontal size L
            L = (n - 1) * 100.0
            X, Y, Z = xyz_grid(range(-L / 2, L / 2, n), range(-L / 2, L / 2, n), range(-h, 0, 11))
            dg = fwd(X, Y, Z, fill(ρ, size(X)))
            c = (n + 1) ÷ 2
            @test dg[c, c] ≈ prism(-L / 2, L / 2, -L / 2, L / 2, 0.0, h, ρ) rtol = 0.04
            push!(err, abs(dg[c, c] - slab))
        end
        @test issorted(err, rev = true)             # converges towards the infinite slab

        # observation surface zRef = 200 m above the top of the slab, with air (zero density) in between,
        # and a grid spacing that varies in x, y and z
        r = range(-800.0, 800.0, 17)
        X, Y, Z = xyz_grid(wiggle(r), wiggle(r), wiggle(range(-h, 200.0, 13)))
        X[9, :, :] .= 0.0                           # station at the origin
        Y[:, 9, :] .= 0.0
        Z[:, :, 11] .= 0.0                          # top of the slab
        dg = fwd(X, Y, Z, ρ .* (Z .<= 0))
        dgz = voxel_grav(X, Y, Z .+ 1000, ρ .* (Z .<= 0), refMod = zeros(13), zRef = 1200, rhoTol = 0, outName = nothing, printing = false)[1]
        @test dg[9, 9] ≈ prism(-800.0, 800.0, -800.0, 800.0, 0.0, h, ρ) rtol = 0.04    # zRef = 0: air above is ignored
        # the cell above the top is half rock, half air: it has density ρ/2
        δ = Z[1, 1, 12]
        @test dgz[9, 9] ≈ prism(-800.0, 800.0, -800.0, 800.0, 200.0, h + 200, ρ) + prism(-800.0, 800.0, -800.0, 800.0, 200.0 - δ, 200.0, ρ / 2) rtol = 0.04
    end

    # random model on a small grid, partly above sea level
    X, Y, Z = xyz_grid(range(-5.0e3, 5.0e3, 9), range(-4.0e3, 4.0e3, 7), range(-6.0e3, 1.0e3, 8))
    RHO = rand(StableRNG(1), size(X)...) .- 0.5
    R = rand(StableRNG(2), size(X, 1), size(X, 2)) .- 0.5

    @testset "symmetry and sign" begin
        body = @. 300.0 * (X^2 + Y^2 + (Z + 3.0e3)^2 < 2.5e3^2)
        dg = fwd(X, Y, Z, body)
        @test all(>(0), dg)
        @test dg ≈ dg[end:-1:1, :]                  # body is symmetric about x = 0
        @test dg ≈ dg[:, end:-1:1]                  # and about y = 0
        @test fwd(X, Y, Z, -body) ≈ -dg
        @test fwd(X, Y, Z, RHO[end:-1:1, :, :]) ≈ fwd(X, Y, Z, RHO)[end:-1:1, :]
    end

    @testset "station heights" begin
        # model without density anomalies above z = 0, stations between 200 and 1000 m
        Xv, Yv, Zv = xyz_grid(wiggle(range(-5.0e3, 5.0e3, 9)), wiggle(range(-4.0e3, 4.0e3, 7)), range(-6.0e3, 1.0e3, 8))
        m = RHO .* (Zv .<= 0)
        zo = 600 .+ 800 .* R
        dg(z) = voxel_grav(Xv, Yv, Zv, m, refMod = zeros(8), zRef = z, rhoTol = 0, outName = nothing, printing = false)[1]
        @test dg(fill(700.0, size(R))) ≈ dg(700.0)
        dgo = dg(zo)
        @test all(dgo[i] ≈ dg(zo[i])[i] for i in (1, 17, 40, 63))
        @test dgo ≈ voxel_grav(Xv ./ 1000, Yv ./ 1000, Zv ./ 1000, m, refMod = zeros(8), zRef = zo ./ 1000, lengthUnit = "km", outName = nothing, printing = false)[1]
        P(A) = permutedims(A, (2, 1, 3))
        @test permutedims(dgo) ≈ voxel_grav(P(Xv), P(Yv), P(Zv), P(m), refMod = zeros(8), zRef = permutedims(zo), outName = nothing, printing = false)[1]
        @test_throws "zRef has size (2, 2)" dg(zeros(2, 2))
    end

    @testset "adjoint" begin
        # dot-product test: sum(R .* (G * RHO)) == sum((G' * R) .* RHO)
        sens = voxel_grav_adjoint(X, Y, Z, R)
        @test sum(R .* fwd(X, Y, Z, RHO)) ≈ sum(sens .* RHO) rtol = 1.0e-10
        @test all(iszero, sens[:, :, Z[1, 1, :] .> 0])

        # same for km coordinates and for x varying along the second dimension
        @test voxel_grav_adjoint(X ./ 1000, Y ./ 1000, Z ./ 1000, R, lengthUnit = "km") ≈ sens
        P(A) = permutedims(A, (2, 1, 3))
        @test voxel_grav_adjoint(P(X), P(Y), P(Z), permutedims(R)) ≈ P(sens)

        # variable spacing and observation surface above the top of the model
        Xv, Yv, Zv = xyz_grid(wiggle(range(-5.0e3, 5.0e3, 9)), wiggle(range(-4.0e3, 4.0e3, 7)), wiggle(range(-6.0e3, 1.0e3, 8)))
        fv(m) = voxel_grav(Xv, Yv, Zv, m, refMod = zeros(8), zRef = 2.0e3, rhoTol = 0, outName = nothing, printing = false)[1]
        @test sum(R .* fv(RHO)) ≈ sum(voxel_grav_adjoint(Xv, Yv, Zv, R, zRef = 2.0e3) .* RHO) rtol = 1.0e-10

        # stations at different heights (e.g. on topography)
        zo = 2.0e3 .+ 500 .* R
        fz(m) = voxel_grav(Xv, Yv, Zv, m, refMod = zeros(8), zRef = zo, rhoTol = 0, outName = nothing, printing = false)[1]
        @test sum(R .* fz(RHO)) ≈ sum(voxel_grav_adjoint(Xv, Yv, Zv, R, zRef = zo) .* RHO) rtol = 1.0e-10

        # sensitivity of one station: δg at that station = sum(sens .* δRHO)
        e = zero(R)
        e[3, 4] = 1
        @test sum(voxel_grav_adjoint(X, Y, Z, e) .* RHO) ≈ fwd(X, Y, Z, RHO)[3, 4]

        # automatic differentiation (as used by Turing) through voxel_grav gives the same gradient as the adjoint
        body(θ) = θ[1] .* exp.(-((X .- θ[2]) .^ 2 .+ Y .^ 2 .+ (Z .- θ[3]) .^ 2) ./ 2.0e3^2)
        θ = [300.0, 500.0, -3.0e3]
        @test ForwardDiff.gradient(θ -> sum(R .* fwd(X, Y, Z, body(θ))), θ) ≈ ForwardDiff.jacobian(body, θ)' * vec(sens)

        # misfit gradient vs. finite difference
        dg_obs = fwd(X, Y, Z, RHO)
        J(m) = sum((fwd(X, Y, Z, m) - dg_obs) .^ 2) / 2
        m = zero(RHO)
        grad = voxel_grav_adjoint(X, Y, Z, fwd(X, Y, Z, m) - dg_obs)
        δ = RHO ./ 100
        @test (J(m + δ) - J(m - δ)) / 2 ≈ sum(grad .* δ) rtol = 1.0e-8
    end

    @testset "KernelAbstractions backend" begin
        # same results with the parallel kernels (CPU threads here; GPUs work the same way)
        Xv, Yv, Zv = xyz_grid(wiggle(range(-5.0e3, 5.0e3, 9)), wiggle(range(-4.0e3, 4.0e3, 7)), wiggle(range(-6.0e3, 1.0e3, 8)))
        zo = 2.0e3 .+ 500 .* R
        g(; kw...) = voxel_grav(Xv, Yv, Zv, RHO; refMod = zeros(8), zRef = zo, rhoTol = 0.1, outName = nothing, printing = false, kw...)
        s(; kw...) = voxel_grav_adjoint(Xv, Yv, Zv, R, zRef = zo; kw...)
        @test all(g(backend = CPU()) .≈ g())
        @test g(backend = CPU(), float_type = Float32)[1] ≈ g()[1] rtol = 1.0e-4
        @test s(backend = CPU()) ≈ s()
        @test s(backend = CPU(), float_type = Float32) ≈ s() rtol = 1.0e-4
        @test_throws "requires KernelAbstractions.jl" g(backend = :gpu)
    end

    @testset "GeoData (ellipsoid) and CartData (flat Earth)" begin
        # spherical cap (5° radius, 10 km thick, 300 kg/m³) on a sphere, station 20 km above its centre, vs. numerical integration
        R0, t, H, α, ρ = 6371.0e3, 10.0, 20.0, 5.0, 300.0
        Lon, Lat, Depth = lonlatdepth_grid(-6:0.2:6, -6:0.2:6, -t:2:0)
        Grid = GeoData(Lon, Lat, Depth, (z = Depth,))
        cap = ρ .* (cosd.(Lat) .* cosd.(Lon) .>= cosd(α))
        dg, = voxel_grav(Grid, cap; refMod = zeros(6), zRef = H, ellipsoid = (a = R0, f = 0.0))
        Rs, nθ, nr = R0 + H * 1.0e3, 4000, 20
        θs = ((1:nθ) .- 0.5) .* deg2rad(α) / nθ
        rs = R0 .- t * 1.0e3 .+ ((1:nr) .- 0.5) .* t * 1.0e3 / nr
        exact = G * ρ * sum(2π * r^2 * sin(θ) * (Rs - r * cos(θ)) / (Rs^2 + r^2 - 2Rs * r * cos(θ))^1.5 for θ in θs, r in rs) * deg2rad(α) / nθ * t * 1.0e3 / nr * 1.0e5
        @test dg[31, 31] ≈ exact rtol = 1.0e-2     # 22 km wide cells, 20-30 km below the station

        # flat Earth (CartData) misses the curvature: 4% less at the centre of the cap
        X, Y, Z = xyz_grid((-6:0.2:6) .* deg2rad(1) * R0 / 1.0e3, (-6:0.2:6) .* deg2rad(1) * R0 / 1.0e3, -t:2:0)
        dgF, = voxel_grav(CartData(X, Y, Z, (z = Z,)), cap; refMod = zeros(6), zRef = H)
        @test dgF ≈ voxel_grav(X, Y, Z, cap; refMod = zeros(6), zRef = H, lengthUnit = "km", outName = nothing, printing = false)[1]
        @test 0.03 < (dg[31, 31] - dgF[31, 31]) / dg[31, 31] < 0.05

        # small area: ellipsoid and flat Earth agree
        Lon, Lat, Depth = lonlatdepth_grid(10:0.02:10.5, 45:0.02:45.5, -10:1:0)
        e2 = (1 / 298.257223563) * (2 - 1 / 298.257223563)
        w = 1 - e2 * sind(45.25)^2
        M, N = 6378137.0 * (1 - e2) / w^1.5, 6378137.0 / sqrt(w)
        X = (Lon .- 10.25) .* deg2rad(1) .* N .* cosd(45.25) ./ 1.0e3
        Y = (Lat .- 45.25) .* deg2rad(1) .* M ./ 1.0e3
        body = 300.0 .* (X .^ 2 .+ Y .^ 2 .+ (Depth .+ 5) .^ 2 .< 4)
        Grid = GeoData(Lon, Lat, Depth, (z = Depth,))
        dgG, gxG, gyG = voxel_grav(Grid, body; refMod = zeros(11), gradients = true)
        dgC, gxC, gyC = voxel_grav(CartData(X, Y, Depth, (z = Depth,)), body; refMod = zeros(11), gradients = true)
        @test maximum(abs, dgG - dgC) < 3.0e-3 * maximum(dgC)
        @test maximum(abs, gxG - gxC) < 1.0e-2 * maximum(abs, gxC)      # east = x
        @test maximum(abs, gyG - gyC) < 1.0e-2 * maximum(abs, gyC)      # north = y

        # adjoint, stations at different heights, parallel kernels
        m = rand(StableRNG(5), size(Lon)...) .- 0.5
        r = rand(StableRNG(6), size(Lon, 1), size(Lon, 2)) .- 0.5
        zo = 0.5 .+ r
        f(; kw...) = voxel_grav(Grid, m; refMod = zeros(11), zRef = zo, rhoTol = 0, kw...)[1]
        sens(; kw...) = voxel_grav_adjoint(Grid, r; zRef = zo, kw...)
        @test sum(r .* f()) ≈ sum(sens() .* m) rtol = 1.0e-10
        @test f(backend = CPU()) ≈ f()
        @test sens(backend = CPU(), float_type = Float32) ≈ sens() rtol = 1.0e-4
        @test sum(r .* voxel_grav(CartData(X, Y, Depth, (z = Depth,)), m; refMod = zeros(11))[1]) ≈ sum(voxel_grav_adjoint(CartData(X, Y, Depth, (z = Depth,)), r) .* m) rtol = 1.0e-10
        # type stability (GeoData and CartData fields are not concretely typed)
        # the return type is one of two tuples, depending on `gradients`
        M3, M1 = Tuple{Matrix{Float64}, Matrix{Float64}, Matrix{Float64}}, Tuple{Matrix{Float64}, Nothing, Nothing}
        @test (@inferred M3 voxel_grav(Grid, m; refMod = zeros(11), zRef = zo)) isa M1
        @test (@inferred M1 voxel_grav(Grid, m; refMod = zeros(11), gradients = true)) isa M3
        @test (@inferred voxel_grav_adjoint(Grid, r)) isa Array{Float64, 3}
        @test (@inferred M3 voxel_grav(CartData(X, Y, Depth, (z = Depth,)), m; refMod = zeros(11))) isa M1

        # cutoff: the windows of the serial code (in lon/lat) give the same result as the masked parallel kernels,
        # also at high latitude; the adjoint uses the same cutoff
        for lat0 in (45.0, 80.0)
            L1, L2, D3 = lonlatdepth_grid(10:0.5:20, lat0:0.25:(lat0 + 5), -20:5:0)
            Gc = GeoData(L1, L2, D3, (z = D3,))
            mc = rand(StableRNG(11), size(L1)...)
            rc = rand(StableRNG(12), size(L1, 1), size(L1, 2))
            serial = voxel_grav(Gc, mc; refMod = zeros(5), rhoTol = 0, cutoff = 100, gradients = true)
            @test all(serial .≈ voxel_grav(Gc, mc; refMod = zeros(5), rhoTol = 0, cutoff = 100, gradients = true, backend = CPU()))
            @test !(serial[1] ≈ voxel_grav(Gc, mc; refMod = zeros(5), rhoTol = 0)[1])
            @test sum(rc .* serial[1]) ≈ sum(voxel_grav_adjoint(Gc, rc; cutoff = 100) .* mc) rtol = 1.0e-10
            @test voxel_grav_adjoint(Gc, rc; cutoff = 100, backend = CPU()) ≈ voxel_grav_adjoint(Gc, rc; cutoff = 100)
        end
        @test_throws "Grid must be a regular lon/lat/depth grid" voxel_grav(GeoData(Lon, Lat, Depth .+ rand(size(Depth)...), (z = Depth,)), m)
    end

    @testset "adaptive subdivision" begin
        # 1 km thick slab directly below the station with 4 km wide cells: point masses are 57% off, subdivided cells < 1%
        L, h, ρ = 32.0e3, 1.0e3, 300.0
        X, Y, Z = xyz_grid(range(-L / 2, L / 2, 9), range(-L / 2, L / 2, 9), range(-h, 0, 11))
        g(s) = voxel_grav(X, Y, Z, fill(ρ, size(X)), refMod = zeros(11), outName = nothing, printing = false, subdivide = s)[1][5, 5]
        x1, x2, z2 = -L / 2, L / 2, h
        exact = 0.0
        for (i, x) in enumerate((x1, x2)), (j, y) in enumerate((x1, x2)), (k, z) in enumerate((0.0, z2))
            r = sqrt(x^2 + y^2 + z^2)
            exact += -(-1)^(i + j + k) * G * ρ * (x * log(y + r) + y * log(x + r) - z * atan(x * y / (z * r))) * 1.0e5
        end
        @test abs(g(0) - exact) > 0.5 * exact
        @test g(2.5) ≈ exact rtol = 1.0e-3
        # the adjoint includes the subdivision
        m = rand(StableRNG(13), size(X)...)
        r = rand(StableRNG(14), 9, 9)
        zo = 300 .* r
        dg = voxel_grav(X, Y, Z, m, refMod = zeros(11), zRef = zo, rhoTol = 0, outName = nothing, printing = false)[1]
        @test sum(r .* dg) ≈ sum(voxel_grav_adjoint(X, Y, Z, r; zRef = zo) .* m) rtol = 1.0e-10
    end

    @testset "arbitrary stations" begin
        # CartData: stations at the grid columns give the same result as the grid version
        X, Y, Z = xyz_grid(wiggle(range(-10.0, 10.0, 11)), wiggle(range(-8.0, 8.0, 9)), range(-6.0, 0.0, 7))
        Grid = CartData(X, Y, Z, (z = Z,))
        m = rand(StableRNG(15), size(X)...)
        zo = 0.2 .+ 0.3 .* rand(StableRNG(16), 11, 9)
        St = CartData(X[:, :, 1:1], Y[:, :, 1:1], reshape(zo, 11, 9, 1), (z = reshape(zo, 11, 9, 1),))
        a = voxel_grav(Grid, m; refMod = zeros(7), zRef = zo, gradients = true, cutoff = 6)
        b = voxel_grav(Grid, m; refMod = zeros(7), stations = St, gradients = true, cutoff = 6)
        @test all(vec.(a) .≈ vec.(b))

        # scattered stations: gradients equal finite differences of the anomaly (moving the stations by ±1 m)
        sx, sy, sz = 16 .* rand(StableRNG(17), 30) .- 8, 12 .* rand(StableRNG(18), 30) .- 6, 0.1 .+ rand(StableRNG(19), 30)
        S(dx, dy) = CartData(reshape(sx .+ dx, :, 1, 1), reshape(sy .+ dy, :, 1, 1), reshape(sz, :, 1, 1), (z = reshape(sz, :, 1, 1),))
        dg, gx, gy = voxel_grav(Grid, m; refMod = zeros(7), stations = S(0, 0), gradients = true, subdivide = 0)
        δ = 1.0e-3   # km
        @test gx ≈ (voxel_grav(Grid, m; refMod = zeros(7), stations = S(δ, 0), subdivide = 0)[1] - voxel_grav(Grid, m; refMod = zeros(7), stations = S(-δ, 0), subdivide = 0)[1]) / 2.0e3δ rtol = 1.0e-5
        @test gy ≈ (voxel_grav(Grid, m; refMod = zeros(7), stations = S(0, δ), subdivide = 0)[1] - voxel_grav(Grid, m; refMod = zeros(7), stations = S(0, -δ), subdivide = 0)[1]) / 2.0e3δ rtol = 1.0e-5

        # adjoint, cutoff and parallel kernels with scattered stations
        r = rand(StableRNG(20), size(dg)...)
        f(; kw...) = voxel_grav(Grid, m; refMod = zeros(7), stations = S(0, 0), rhoTol = 0, cutoff = 6, gradients = true, kw...)
        @test sum(r .* f()[1]) ≈ sum(voxel_grav_adjoint(Grid, r; stations = S(0, 0), cutoff = 6) .* m) rtol = 1.0e-10
        @test all(f(backend = CPU()) .≈ f())
        @test voxel_grav_adjoint(Grid, r; stations = S(0, 0), cutoff = 6, backend = CPU()) ≈ voxel_grav_adjoint(Grid, r; stations = S(0, 0), cutoff = 6)

        # GeoData: the same checks on the ellipsoid
        Lon, Lat, Depth = lonlatdepth_grid(10:0.05:10.5, 45:0.05:45.4, -6:1:0)
        GridG = GeoData(Lon, Lat, Depth, (z = Depth,))
        mg = rand(StableRNG(21), size(Lon)...)
        StG = GeoData(Lon[:, :, 1:1], Lat[:, :, 1:1], fill(0.3, 11, 9, 1), (z = fill(0.3, 11, 9, 1),))
        ag = voxel_grav(GridG, mg; refMod = zeros(7), zRef = 0.3, gradients = true, cutoff = 10)
        bg = voxel_grav(GridG, mg; refMod = zeros(7), stations = StG, gradients = true, cutoff = 10)
        @test all(vec.(ag) .≈ vec.(bg))
        slon, slat = 10.05 .+ 0.4 .* rand(StableRNG(22), 20), 45.05 .+ 0.3 .* rand(StableRNG(23), 20)
        SG = GeoData(reshape(slon, :, 1, 1), reshape(slat, :, 1, 1), fill(0.5, 20, 1, 1), (z = fill(0.5, 20, 1, 1),))
        fg(; kw...) = voxel_grav(GridG, mg; refMod = zeros(7), stations = SG, rhoTol = 0, cutoff = 10, gradients = true, kw...)
        rg = rand(StableRNG(24), 20, 1, 1)
        @test sum(rg .* fg()[1]) ≈ sum(voxel_grav_adjoint(GridG, rg; stations = SG, cutoff = 10) .* mg) rtol = 1.0e-10
        @test all(fg(backend = CPU()) .≈ fg())
    end

    @testset "gravity corrections" begin
        # GRS80 normal gravity at the equator and the pole, and the free-air gradient (Hofmann-Wellenhof & Moritz, 2006)
        @test normal_gravity(0.0) ≈ 978032.67715
        @test normal_gravity(90.0) ≈ 983218.63685
        @test normal_gravity(45.0) ≈ 980619.92 atol = 0.01
        @test normal_gravity(45.0, 0.001) - normal_gravity(45.0) ≈ -0.3086 atol = 1.0e-3
        # Bouguer slab, also for oceans
        @test bouguer_correction(1.0) ≈ 2π * G * 2670 * 1.0e3 * 1.0e5
        @test bouguer_correction([-4.0, 0.0]) ≈ [2π * G * (1040 - 2670) * 4.0e3 * 1.0e5, 0.0]
        # complete Bouguer correction at the centre of a plateau and an ocean basin, 200 km wide ≈ slab of finite width R
        X, Y, _ = xyz_grid(-100:4:100, -100:4:100, 0)
        for h in (1.0, -4.0)
            topo = fill(h, size(X))
            te, = topographic_effect(CartData(X, Y, topo, (z = topo,)); cutoff = Inf, dz = 0.25)
            @test te[26, 26] ≈ bouguer_correction(h) * (1 - abs(h) / 200) rtol = 0.05
        end
        # GeoData: stations given explicitly give the same as the default stations on the topography
        Lon, Lat, _ = lonlatdepth_grid(9:0.1:10, 44:0.1:45, 0)
        topo = 0.5 .+ 0.5 .* sind.(40 .* Lon) .* cosd.(40 .* Lat)
        T = GeoData(Lon, Lat, topo, (z = topo,))
        te, = topographic_effect(T; cutoff = 50)
        @test vec(te) ≈ vec(topographic_effect(T; cutoff = 50, stations = T)[1])
        @test all(>(0), te)
    end
end
