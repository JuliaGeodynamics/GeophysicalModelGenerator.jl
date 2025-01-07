using GeophysicalModelGenerator

# Grid parameters
nx, ny, nz = 512, 512, 128
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
    ((-750.0, 200.0), (-750.0, 1000.0))  # Segment 3  #It is possible to add more segments
]

lith = LithosphericPhases(Layers=[15 55], Phases=[1 2], Tlab=1250)
add_box!(Phases, Temp, Grid; xlim=(-1000.0,0.0), ylim=(-500.0, 500.0), zlim=(-80.0, 0.0), phase = lith, T=SpreadingRateTemp(SpreadingVel=3), segments=segments);

# Add and save results
Grid = addfield(Grid, (; Phases, Temp))
write_paraview(Grid, "Ridge_Thermal_Structure_test_2")
