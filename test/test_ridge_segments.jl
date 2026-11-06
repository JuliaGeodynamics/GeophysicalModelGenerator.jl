using GeophysicalModelGenerator
using Test

@testset "Ridge Thermal Structure Tests" begin
    # Grid parameters
    x = -1000.0:100.0:1000.0
    y = -1000.0:50.0:1000.0
    z = -100.0:10.0:0.0
    nx, ny, nz = length(x), length(y), length(z)
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

    ix, iy, iz = findfirst(==(-500.0), x), findfirst(==(-250.0), y), findfirst(==(-10.0), z)
    @test Temp[ix, iy, iz] ≈ 1350.0 atol = 1.0e-3   # on the ridge (segment 1)
    @test Temp[ix + 3, iy, iz] < 1000.0           # 300 km off the ridge

    @test minimum(Temp) >= 0.0  # Minimum temperature
    @test maximum(Temp) <= 1350.0  # Maximum temperature
    @test all(≥(0),    Temp) # Ensure no negative temperatures
    @test all(≤(1350), Temp) # Ensure no temperatures above max

    # Check if phases are correctly assigned in expected regions
    @test first(Phases) == 2  # Example: Verify a point's phase
    @test last(Phases)  == 2  # Example: Verify another point's phase
end
