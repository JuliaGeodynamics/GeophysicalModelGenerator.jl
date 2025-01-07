using GeophysicalModelGenerator
using Test

@testset "Ridge Thermal Structure Test" begin
    # Grid parameters
    nx, ny, nz = 128, 128, 64
    x = range(-1000, 1000, length=nx)
    y = range(-1000, 1000, length=ny)
    z = range(-660, 0, length=nz)
    Grid = CartData(xyz_grid(x, y, z))

    # Matrices
    Phases = fill(2, nx, ny, nz)
    Temp = fill(1350.0, nx, ny, nz)

    # Ridge segments
    segments = [
        ((-500.0, -1000.0), (-500.0, 0.0)), 
        ((-250.0, 0.0), (-250.0, 200.0)),
        ((-750.0, 200.0), (-750.0, 1000.0))
    ]

    # Parameters
    params = SpreadingRateTemp(
        Tsurface=0.0,
        Tmantle=1350.0,
        Adiabat=0.0,
        SpreadingVel=3.0,
        AgeRidge=0.0,
        maxAge=100,
        segments=segments
    )

    # Call function
    compute_thermal_structure(Temp, x, y, z, Phases, params)

    # Test verification
    @test minimum(Temp) >= 0.0  # Minimum temperature
    @test maximum(Temp) <= 1350.0  # Maximum temperature

end

