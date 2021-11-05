# test LaMEM i/o routines
using Test, GeophysicalModelGenerator

# Load LaMEM input file:
Grid        =   ReadLaMEM_InputFile("test_files/SaltModels.dat")
@test Grid.X[10] ≈ -2.40625

# Transfer into ParaviewData struct:
Phases      =   zeros(Int32,   size(Grid.X));
Temp        =   zeros(Float64, size(Grid.X));
Model3D     =   ParaviewData(Grid, (Phases=Grid.Z,));
@test  Model3D.y.val[100]==-1.9375

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

Model3D     =   ParaviewData(Grid, (Phases=Phases,Temp=Temp));
Write_Paraview(Model3D,"LaMEM_ModelSetup")                  # Save model to paraview    



