# test LaMEM i/o routines
using Test, GeophysicalModelGenerator

# Load LaMEM input file:
Grid        =   ReadLaMEM_InputFile("test_files/SaltModels.dat")
@test Grid.X[10] ≈ -2.40625

# Transfer into ParaviewData struct:
Phases      =   zeros(Int32,   size(Grid.X));
Temp        =   zeros(Float64, size(Grid.X));
Model3D     =   CartData(Grid, (Phases=Grid.Z,));
@test  Value(Model3D.y[100])==-1.9375km

# Read Partitioning file: 
PartitioningFile = "test_files/ProcessorPartitioning_4cpu_1.2.2.bin"
Nprocx,Nprocy,Nprocz, xc,yc,zc, nNodeX,nNodeY,nNodeZ = GetProcessorPartitioning(PartitioningFile)
@test Nprocz==2
@test yc[2]==0.0

# Save serial output
Save_LaMEMMarkersParallel(Model3D, verbose=false)

# Save parallel output
Save_LaMEMMarkersParallel(Model3D, PartitioningFile=PartitioningFile, verbose=false)

# Test creating model setups
Grid        =   ReadLaMEM_InputFile("test_files/Subduction2D_FreeSlip_Particles_Linear_DirectSolver.dat")
Phases      =   zeros(Int32,   size(Grid.X));

# constant T
Temp        =   ones(Float64, size(Grid.X))*1350;
AddBox!(Phases,Temp,Grid, xlim=(0,500), zlim=(-50,0), phase=ConstantPhase(3), DipAngle=10, T=ConstantTemp(1000))
@test sum(Temp) == 1.1905107e9

# Add a layer above the slab with a different phase but no thermal structure
AddBox!(Phases,Temp,Grid, xlim=(0,500), zlim=(0,20), phase=ConstantPhase(1), DipAngle=10, Origin=(0,0,0))
@test sum(Temp) == 1.1905107e9

# Linear T
Temp        =   ones(Float64, size(Grid.X))*1350;
AddBox!(Phases,Temp,Grid, xlim=(0,500), zlim=(-50,0), phase=ConstantPhase(3), DipAngle=10, T=LinearTemp(Tbot=1350, Ttop=200))
@test sum(Temp) == 1.1881296265169694e9

# Halfspace cooling T structure
Phases      =   zeros(Int32,   size(Grid.X));
Temp        =   ones(Float64, size(Grid.X))*1350;
AddBox!(Phases,Temp,Grid, xlim=(0,500), zlim=(-500,0), phase=LithosphericPhases(Layers=[15 15 250], Phases=[1 2 3 0], Tlab=1250), DipAngle=10, T=HalfspaceCoolingTemp(Age=20, Adiabat=0.3))
@test sum(Temp) == 1.1942982365477426e9

# Mid-oceanic ridge cooling temperature structure
Phases      =   zeros(Int32,   size(Grid.X));
Temp        =   ones(Float64, size(Grid.X))*1350;
AddBox!(Phases,Temp,Grid, xlim=(0,500), zlim=(-500,-20), phase=LithosphericPhases(Layers=[15 15 250], Phases=[1 2 3 0],Tlab=1250), DipAngle=10, T=SpreadingRateTemp(MORside="right", SpreadingVel=3))
@test sum(Temp) == 1.189394358568891e9

Model3D     =   CartData(Grid, (Phases=Phases,Temp=Temp));
Write_Paraview(Model3D,"LaMEM_ModelSetup")                  # Save model to paraview    


# Test writing a LaMEM topography file
X,Y,Z = XYZGrid(-20:20,-10:10,0);
Z = cos.(2*pi.*X./5).*cos.(2*pi.*Y./10)

Topo = CartData(X,Y,Z,(Topography=Z,))
@test Save_LaMEMTopography(Topo, "test_topo.dat")==nothing
rm("test_topo.dat")


# Test adding geometric primitives
Grid    = ReadLaMEM_InputFile("test_files/GeometricPrimitives.dat")
Phases  = zeros(Int32,size(Grid.X));
Temp    = zeros(Float64,size(Grid.X));
AddSphere!(Phases,Temp,Grid, cen=(0,0,-6), radius=2, phase=ConstantPhase(1), T=ConstantTemp(800))
@test Phases[55,55,55] == 1
@test Phases[56,56,56] == 0
@test Temp[44,52,21]   == 800.0
@test Temp[44,52,20]   == 0.0

AddEllipsoid!(Phases,Temp,Grid, cen=(-2,-1,-7), axes=(1,2,3), StrikeAngle=90, DipAngle=45, phase=ConstantPhase(2), T=ConstantTemp(600))
@test Phases[11,37,28] == 2
@test Phases[10,37,28] == 0
@test Temp[31,58,18]   == 600.0
@test Temp[31,59,18]   == 0.0

AddCylinder!(Phases,Temp,Grid, base=(0,0,-5), cap=(3,3,-2), radius=1.5, phase=ConstantPhase(3), T=ConstantTemp(400))
@test Phases[55,65,75] == 3
@test Phases[54,65,75] == 0
@test Temp[55,46,45]   == 400.0
@test Temp[55,45,45]   == 800.0