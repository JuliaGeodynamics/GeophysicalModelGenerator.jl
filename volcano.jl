using GeophysicalModelGenerator

nx,ny,nz = 128, 128, 128 
x = range(-100, 100, nx);
y = range(-100, 100, ny);
z = range(-110,  50, nz);
Grid = CartData(xyz_grid(x,y,z));

# Now we create an integer array that will hold the `Phases` information (which usually refers to the material or rock type in the simulation)
Phases = fill(0, nx, ny, nz);

# In many (geodynamic) models, one also has to define the temperature, so lets define it as well
Temp = fill(1350.0, nx, ny, nz);

lith = LithosphericPhases(Layers=[15 45 100], Phases=[1 2 3])

# And an an inclined part:
add_box!(Phases, Temp, Grid; 
    xlim=(-100, 100), 
    ylim=(-400, 400.0), 
    zlim=(-110.0, 0.0), 
    phase = lith, 
    Origin = (0,0,0), 
    T = HalfspaceCoolingTemp(Age=20));
    

# # Add them to the `CartData` dataset:
Grid = addfield(Grid,(;Phases, Temp))
# heatmap(x, z, Grid.fields.Temp[:,64,:])

# Which looks like
# write_paraview(Grid,"Grid3D_FreeSubduction");

center     = (0, 0, 0)
height     = 10 # [kmm]
radius     = 15 # [km]
crater     = 0.0
base       = 0.0
background = nothing


add_volcano!(Phases, Temp, Grid;
    volcanic_phase  = 1,
    center     = (0, 0, 0),
    height     = 10,
    radius     = 15,
    crater     = 0.0,
    base       = 0.0,
    background = nothing,
    T = HalfspaceCoolingTemp(Age=20)
)

Grid = addfield(Grid,(;Phases, Temp))

write_paraview(Grid,"Grid3D_FreeSubduction");

Grid.fields.Temp[Grid.fields.Temp.==0] .= NaN
heatmap(x, z, Grid.fields.Temp[:,64,:])
heatmap(x, z, Grid.fields.Phases[:,64,:])

# heatmap(x, z, Grid.z.val[:,64,:])
# heatmap(x, z, depth[:,64,:])
# heatmap(x, z, ind[:,64,:])

# Grid.fields.Temp[:,:,1]
# lines(Grid.fields.Temp[64,64,:], depth[64,64,:])
