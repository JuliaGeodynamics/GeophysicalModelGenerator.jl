# test LaMEM I/O routines
using Test, GeophysicalModelGenerator

# Load LaMEM input file with grid refinement:
Grid = read_LaMEM_inputfile("test_files/non-uniform_grid.dat")
@test Grid.X[2358] ≈ 1.833333333333335
@test Grid.Y[7741] ≈ -0.45
@test Grid.Z[5195] ≈ -8.897114711471147

# Non-uniform grid in z-direction (taken from LaMEM test suite)
Grid = read_LaMEM_inputfile("test_files/Subduction_VEP.dat")
@test maximum(diff(Grid.z_vec)) ≈ 2.626262626262701
@test minimum(diff(Grid.z_vec)) ≈ 1.3333333333333321
@test maximum(diff(Grid.x_vec)) ≈ 10.4166666666667
@test minimum(diff(Grid.x_vec)) ≈ 10.4166666666667

# Add command-line arguments
args = "-coord_z -660,-300,5,25 -nel_z 36,90,5 -nel_x 32"
coord_z = GeophysicalModelGenerator.ParseValue_LaMEM_InputFile("test_files/Subduction_VEP.dat", "coord_z", Float64, args = args)
nel_z = GeophysicalModelGenerator.ParseValue_LaMEM_InputFile("test_files/Subduction_VEP.dat", "nel_z", Int64, args = args)
@test coord_z ≈ [-660.0, -300.0, 5.0, 25.0]
@test nel_z == [36, 90, 5]

Grid = read_LaMEM_inputfile("test_files/Subduction_VEP.dat", args = args)
@test Grid.nel_z == 131

# Load LaMEM input file:
Grid = read_LaMEM_inputfile("test_files/SaltModels.dat")
@test Grid.X[10] ≈ -2.40625

# Transfer into ParaviewData struct:
Phases = zeros(Int32, size(Grid.X));
Temp = zeros(Float64, size(Grid.X));
Model3D = CartData(Grid, (Phases = Grid.Z,));
@test  Value(Model3D.y[100]) == -1.9375km

# Read Partitioning file:
PartitioningFile = "test_files/ProcessorPartitioning_4cpu_1.2.2.bin"
Nprocx, Nprocy, Nprocz, xc, yc, zc, nNodeX, nNodeY, nNodeZ = get_processor_partitioning(PartitioningFile)
@test Nprocz == 2
@test yc[2] == 0.0

# Save serial output
save_LaMEM_markers_parallel(Model3D, verbose = false)

# Save parallel output
save_LaMEM_markers_parallel(Model3D, PartitioningFile = PartitioningFile, verbose = false)

# Get partitioning without generating partitioning file
n_ranks = 8; nx=128; ny=64; nz=32;
xcoords_ans=[-35.35034129014803, -19.867446943762857, -4.3845525973776835, 11.098341749007488, 26.58123609539267];
ycoords_ans=[-24.510578171895816, 5.122856149231547, 34.75629047035891];
zcoords_ans=[-6.4, 6.4];
nProcX_ans = 4; nProcY_ans = 2; nProcZ_ans = 1;
nNodeX_ans = 129; nNodeY_ans = 65; nNodeZ_ans = 33;
P         = setup_model_domain([extrema(xcoords_ans)...], [extrema(ycoords_ans)...], [extrema(zcoords_ans)...], nx, ny, nz, n_ranks)
@test isapprox(nProcX_ans, P.nProcX; atol=1e-8)
@test isapprox(nProcY_ans, P.nProcY; atol=1e-8)
@test isapprox(nProcZ_ans, P.nProcZ; atol=1e-8)
@test isapprox(xcoords_ans, P.xc; atol=1e-8)
@test isapprox(ycoords_ans, P.yc; atol=1e-8)
@test isapprox(zcoords_ans, P.zc; atol=1e-8)
@test isapprox(nNodeX_ans, P.nNodeX; atol=1e-8)
@test isapprox(nNodeY_ans, P.nNodeY; atol=1e-8)
@test isapprox(nNodeZ_ans, P.nNodeZ; atol=1e-8)

n_ranks = 128; nx=256; ny=1; nz=128;
xcoords_ans=[-35.35034129014803, -31.479617703551735, -27.608894116955444, -23.73817053035915, -19.867446943762857, -15.996723357166562, -12.125999770570267, -8.255276183973974, -4.384552597377681, -0.5138290107813872, 3.3568945758149065, 7.227618162411201, 11.098341749007494, 14.969065335603787, 18.839788922200082, 22.710512508796374, 26.58123609539267];
ycoords_ans=[-0.1, 0.1];
zcoords_ans=[-6.4, -4.8, -3.2, -1.6, 0.0, 1.6, 3.2, 4.8, 6.4];
nProcX_ans = 16; nProcY_ans = 1; nProcZ_ans = 8;
nNodeX_ans = 257; nNodeY_ans = 2; nNodeZ_ans = 129;
P         = setup_model_domain([extrema(xcoords_ans)...], [extrema(ycoords_ans)...], [extrema(zcoords_ans)...], nx, ny, nz, n_ranks)
@test isapprox(nProcX_ans, P.nProcX; atol=1e-8)
@test isapprox(nProcY_ans, P.nProcY; atol=1e-8)
@test isapprox(nProcZ_ans, P.nProcZ; atol=1e-8)
@test isapprox(xcoords_ans, P.xc; atol=1e-8)
@test isapprox(ycoords_ans, P.yc; atol=1e-8)
@test isapprox(zcoords_ans, P.zc; atol=1e-8)
@test isapprox(nNodeX_ans, P.nNodeX; atol=1e-8)
@test isapprox(nNodeY_ans, P.nNodeY; atol=1e-8)
@test isapprox(nNodeZ_ans, P.nNodeZ; atol=1e-8)

n_ranks = 2048; nx=512; ny=2048; nz=512;
xcoords_ans=[-35.35034129014803, -27.60889411695544, -19.867446943762857, -12.12599977057027, -4.3845525973776835, 3.356894575814903, 11.098341749007488, 18.839788922200068, 26.58123609539267 ];
ycoords_ans=[-24.510578171895816, -22.658488526825355, -20.806398881754898, -18.954309236684434, -17.102219591613977, -15.250129946543517, -13.398040301473054, -11.545950656402596,  -9.693861011332135,  -7.841771366261674];
zcoords_ans=[-6.4, -4.800000000000001, -3.2, -1.6, 0.0, 1.6, 3.2, 4.800000000000001, 6.4];
nProcX_ans = 8; nProcY_ans = 32; nProcZ_ans = 8;
nNodeX_ans = 513; nNodeY_ans = 2049; nNodeZ_ans = 513;
P         = setup_model_domain([-35.35034129014803,26.58123609539267], [-24.510578171895816,34.75629047035891], [-6.4,6.4], nx, ny, nz, n_ranks)
@test isapprox(nProcX_ans, P.nProcX; atol=1e-8)
@test isapprox(nProcY_ans, P.nProcY; atol=1e-8)
@test isapprox(nProcZ_ans, P.nProcZ; atol=1e-8)
@test isapprox(xcoords_ans, P.xc; atol=1e-8)
@test isapprox(ycoords_ans, P.yc[1:10]; atol=1e-8)
@test isapprox(zcoords_ans, P.zc; atol=1e-8)
@test isapprox(nNodeX_ans, P.nNodeX; atol=1e-8)
@test isapprox(nNodeY_ans, P.nNodeY; atol=1e-8)
@test isapprox(nNodeZ_ans, P.nNodeZ; atol=1e-8)

n_ranks = 32768; nx=2048; ny=512; nz=1024;
xcoords_ans=[-35.35034129014803, -34.38266039349895, -33.41497949684988, -32.4472986002008, -31.47961770355173, -30.511936806902664, -29.54425591025359, -28.576575013604515, -27.60889411695544, -26.64121322030637 ];
ycoords_ans=[ -24.510578171895816, -20.806398881754898, -17.102219591613977, -13.398040301473054,  -9.693861011332135,  -5.989681721191215,  -2.285502431050293,   1.418676859090624,   5.122856149231547,   8.82703543937247];
zcoords_ans=[-6.4, -6.0, -5.6000000000000005, -5.2, -4.800000000000001, -4.4, -4.0, -3.6, -3.2, -2.8000000000000003];
nProcX_ans = 64; nProcY_ans = 16; nProcZ_ans = 32;
nNodeX_ans = 2049; nNodeY_ans = 513; nNodeZ_ans = 1025;
P         = setup_model_domain([-35.35034129014803,26.58123609539267], [-24.510578171895816,34.75629047035891], [-6.4,6.4], nx, ny, nz, n_ranks)
@test isapprox(nProcX_ans, P.nProcX; atol=1e-8)
@test isapprox(nProcY_ans, P.nProcY; atol=1e-8)
@test isapprox(nProcZ_ans, P.nProcZ; atol=1e-8)
@test isapprox(xcoords_ans, P.xc[1:10]; atol=1e-8)
@test isapprox(ycoords_ans, P.yc[1:10]; atol=1e-8)
@test isapprox(zcoords_ans, P.zc[1:10]; atol=1e-8)
@test isapprox(nNodeX_ans, P.nNodeX; atol=1e-8)
@test isapprox(nNodeY_ans, P.nNodeY; atol=1e-8)
@test isapprox(nNodeZ_ans, P.nNodeZ; atol=1e-8)

# Test creating model setups
Grid = read_LaMEM_inputfile("test_files/Subduction2D_FreeSlip_Particles_Linear_DirectSolver.dat")
Phases = zeros(Int32, size(Grid.X));

# constant T
Temp = ones(Float64, size(Grid.X)) * 1350;
add_box!(Phases, Temp, Grid, xlim = (0, 500), zlim = (-50, 0), phase = ConstantPhase(3), DipAngle = 10, T = ConstantTemp(1000))
@test sum(Temp) == 1.1905107e9

# Add a layer above the slab with a different phase but no thermal structure
add_box!(Phases, Temp, Grid, xlim = (0, 500), zlim = (0, 20), phase = ConstantPhase(1), DipAngle = 10, Origin = (0, 0, 0))
@test sum(Temp) == 1.1905107e9

# Linear T
Temp = ones(Float64, size(Grid.X)) * 1350;
add_box!(Phases, Temp, Grid, xlim = (0, 500), zlim = (-50, 0), phase = ConstantPhase(3), DipAngle = 10, T = LinearTemp(Tbot = 1350, Ttop = 200))
@test sum(Temp) == 1.1881296265169694e9

# Halfspace cooling T structure
Phases = zeros(Int32, size(Grid.X));
Temp = ones(Float64, size(Grid.X)) * 1350;
add_box!(Phases, Temp, Grid, xlim = (0, 500), zlim = (-500, 0), phase = LithosphericPhases(Layers = [15 15 250], Phases = [1 2 3 0], Tlab = 1250), DipAngle = 10, T = HalfspaceCoolingTemp(Age = 20, Adiabat = 0.3))
@test sum(Temp) == 1.1942982365477426e9

# Mid-oceanic ridge cooling temperature structure
Phases = zeros(Int32, size(Grid.X));
Temp = ones(Float64, size(Grid.X)) * 1350;
add_box!(Phases, Temp, Grid, xlim = (0, 500), zlim = (-500, -20), phase = LithosphericPhases(Layers = [15 15 250], Phases = [1 2 3 0], Tlab = 1250), DipAngle = 10, T = SpreadingRateTemp(MORside = "right", SpreadingVel = 3))
@test sum(Temp) == 1.189394358568891e9

Model3D = CartData(Grid, (Phases = Phases, Temp = Temp));
write_paraview(Model3D, "LaMEM_ModelSetup")                  # Save model to paraview


# Test writing a LaMEM topography file
X, Y, Z = xyz_grid(-20:20, -10:10, 0);
Z = cos.(2 * pi .* X ./ 5) .* cos.(2 * pi .* Y ./ 10)

Topo = CartData(X, Y, Z, (Topography = Z,))
@test save_LaMEM_topography(Topo, "test_topo.dat") == nothing
rm("test_topo.dat")


# Test adding geometric primitives
Grid = read_LaMEM_inputfile("test_files/GeometricPrimitives.dat")
Phases = zeros(Int32, size(Grid.X));
Temp = zeros(Float64, size(Grid.X));
add_sphere!(Phases, Temp, Grid, cen = (0, 0, -6), radius = 2.0, phase = ConstantPhase(1), T = ConstantTemp(800))
@test Phases[55, 55, 55] == 1
@test Phases[56, 56, 56] == 0
@test Temp[44, 52, 21] == 800.0
@test Temp[44, 52, 20] == 0.0

add_ellipsoid!(Phases, Temp, Grid, cen = (-2, -1, -7), axes = (1, 2, 3), StrikeAngle = 90, DipAngle = 45, phase = ConstantPhase(2), T = ConstantTemp(600))
@test Phases[11, 37, 28] == 2
@test Phases[10, 37, 28] == 0
@test Temp[31, 58, 18] == 600.0
@test Temp[31, 59, 18] == 0.0

add_cylinder!(Phases, Temp, Grid, base = (0, 0, -5), cap = (3, 3, -2), radius = 1.5, phase = ConstantPhase(3), T = ConstantTemp(400))
@test Phases[55, 65, 75] == 3
@test Phases[54, 65, 75] == 0
@test Temp[55, 46, 45] == 400.0
@test Temp[55, 45, 45] == 800.0

# for debugging:
#data = CartData(Grid, (;Phases, Temp));
#write_paraview(data,"test")

# test adding generic volcano topography
Grid = read_LaMEM_inputfile("test_files/SaltModels.dat");
Topo = make_volc_topo(Grid, center = [0.0, 0.0], height = 0.4, radius = 1.5, crater = 0.5, base = 0.1);
@test Topo.fields.Topography[13, 13] ≈ 0.279583654km
Topo = make_volc_topo(Grid, center = [0.0, 0.0], height = 0.8, radius = 0.5, crater = 0.0, base = 0.4, background = Topo.fields.Topography);
@test Topo.fields.Topography[13, 13] ≈ 0.279583654km
@test Topo.fields.Topography[16, 18] ≈ 0.619722436km
