using GeophysicalModelGenerator
using Test

@testset "Ridge Thermal Structure Tests" begin
    # Grid parameters
    nx, ny, nz = 128, 128, 57
    x = range(-1000, 1000, length=nx)
    y = range(-1000, 1000, length=ny)
    z = range(-660, 0, length=nz)
    Grid = CartData(xyz_grid(x, y, z))

    # Phases and temperature
    Phases = fill(2, nx, ny, nz)
    Temp = fill(1350.0, nx, ny, nz)

    # Ridge Segments
    segments = [
        ((-500.0, -1000.0), (-500.0, 0.0)),  # Segment 1
        ((-250.0, 0.0), (-250.0, 200.0)),    # Segment 2
        ((-750.0, 200.0), (-750.0, 1000.0))  # Segment 3
    ]

    # Lithospheric phases
    lith = LithosphericPhases(Layers=[15 55], Phases=[1 2], Tlab=1250)
    add_box!(Phases, Temp, Grid; xlim=(-1000.0,0.0), ylim=(-500.0, 500.0), zlim=(-80.0, 0.0), phase = lith, T=SpreadingRateTemp(SpreadingVel=3), segments=segments)

    # Add and save results
    Grid = addfield(Grid, (; Phases, Temp))
    #write_paraview(Grid, "Ridge_Thermal_Structure_test_2")

    # Test verifications
    @test minimum(Temp) >= 0.0  # Minimum temperature
    @test maximum(Temp) <= 1350.0  # Maximum temperature
    @test all(Temp .>= 0.0)  # Ensure no negative temperatures
    @test all(Temp .<= 1350.0)  # Ensure no temperatures above max

    # Check if phases are correctly assigned in expected regions
    @test Phases[1, 1, 1] == 2  # Example: Verify a point's phase
    @test Phases[end, end, end] == 2  # Example: Verify another point's phase
end
