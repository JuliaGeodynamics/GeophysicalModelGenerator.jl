# test LaMEM i/o routines
using Test, GeophysicalModelGenerator

# Load LaMEM input file:
Grid        =   ReadLaMEM_InputFile("test_files/SaltModels.dat")
@test Grid.X[10] ≈ -2.40625

# Transfer into CartData struct:
Phases      =   zeros(Int32, size(Grid.X));
Model3D     =   CartData(Grid, (Phases=Grid.Z,));
@test  Model3D.y.val[100]==-1.9375

# Read Partitioning file: 
PartitioningFile = "test_files/ProcessorPartitioning_4cpu_1.2.2.bin"
Nprocx,Nprocy,Nprocz, xc,yc,zc, nNodeX,nNodeY,nNodeZ = GetProcessorPartitioning(PartitioningFile)
@test Nprocz==2
@test yc[2]==0.0

# Save serial
Save_LaMEMMarkersParallel(Model3D, verbose=false)

# Save parallel
Save_LaMEMMarkersParallel(Model3D, PartitioningFile=PartitioningFile, verbose=false)

