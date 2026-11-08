# test setting geometries in the different grid types
using Test, GeophysicalModelGenerator, GeoParams, Statistics

# GeoData
Lon3D, Lat3D, Depth3D = lonlatdepth_grid(1.0:1:10.0, 11.0:1:20.0, (-20:1:-10) * km);
Data = zeros(size(Lon3D));
Temp = ones(Float64, size(Data)) * 1350;
Phases = zeros(Int32, size(Data));
Grid = GeoData(Lon3D, Lat3D, Depth3D, (DataFieldName = Data,))

add_box!(Phases, Temp, Grid, xlim = (2, 4), zlim = (-15, -10), phase = ConstantPhase(3), DipAngle = 10, T = LinearTemp(Tbot = 1350, Ttop = 200))
@test sum(Temp[1, 1, :]) ≈ 14850.0

add_ellipsoid!(Phases, Temp, Grid, cen = (4, 15, -17), axes = (1, 2, 3), StrikeAngle = 90, DipAngle = 45, phase = ConstantPhase(2), T = ConstantTemp(600))
@test sum(Temp[1, 1, :]) ≈ 14850.0

Phases = zeros(Int32, size(Data));
empty_bounds = []
add_box!(Phases, Temp, Grid, empty_bounds, phase = ConstantPhase(69))
@test (unique(Phases)) == [0]

not_empty_bounds = [[2, 4], [minimum(Grid.lat.val), maximum(Grid.lat.val)], [-15, -10]]
add_box!(Phases, Temp, Grid, not_empty_bounds, phase = ConstantPhase(69))
@test maximum(Phases) == 69

# CartData
X, Y, Z = xyz_grid(1.0:1:10.0, 11.0:1:20.0, -20:1:-10);
Data = zeros(size(X));
Temp = ones(Float64, size(Data)) * 1350;
Phases = zeros(Int32, size(Data));
Grid = CartData(X, Y, Z, (DataFieldName = Data,))

add_box!(Phases, Temp, Grid, xlim = (2, 4), zlim = (-15, -10), phase = ConstantPhase(3), DipAngle = 10, T = LinearTemp(Tbot = 1350, Ttop = 200))
@test sum(Temp[1, 1, :]) ≈ 14850.0

add_ellipsoid!(Phases, Temp, Grid, cen = (4, 15, -17), axes = (1, 2, 3), StrikeAngle = 90, DipAngle = 45, phase = ConstantPhase(2), T = ConstantTemp(600))
@test sum(Temp[1, 1, :]) ≈ 14850.0

# CartGrid
Grid = create_CartGrid(size = (10, 20, 30), x = (0.0, 10), y = (0.0, 10), z = (2.0, 10))
Temp = ones(Float64, Grid.N...) * 1350;
Phases = zeros(Int32, Grid.N...);

add_box!(Phases, Temp, Grid, xlim = (2, 4), zlim = (4, 8), phase = ConstantPhase(3), DipAngle = 10, T = LinearTemp(Tbot = 1350, Ttop = 200))
@test maximum(Phases) == 3

add_stripes!(Phases, Grid, stripAxes = (1, 1, 1), stripeWidth = 0.2, stripeSpacing = 1, Origin = nothing, StrikeAngle = 0, DipAngle = 10, phase = ConstantPhase(3), stripePhase = ConstantPhase(4))
@test maximum(Phases) == 4

# stripes are centered on multiples of stripeSpacing and have the full stripeWidth
Xs = reshape(collect(-10.0:0.01:10.0), :, 1, 1)
Phs = fill(1, size(Xs))
add_stripes!(Phs, CartData(Xs, zero(Xs), zero(Xs), (Z = Xs,)); stripAxes = (1, 0, 0), stripeWidth = 0.4, stripeSpacing = 2, Origin = (0, 0, 0), DipAngle = 0, phase = ConstantPhase(1), stripePhase = ConstantPhase(2))
@test [Phs[findfirst(≈(x), vec(Xs))] for x in (-2.15, -0.15, 0.15, 1.85, 2.15, 1.0)] == [2, 2, 2, 2, 2, 1]
@test count(Phs .== 2) / length(Phs) ≈ 0.2 atol = 0.01

# Create a CartData structure from it
Data = CartData(Grid, (T = Temp, Phases = Phases))

@test NumValue(Data.x[3, 3, 2]) ≈ 2.2222222222222223

# Doing the same for vertical cross-sections
Grid2D = create_CartGrid(size = (10, 30), x = (0.0, 10), z = (2.0, 10))
Temp2D = ones(Float64, Grid2D.N...) * 1350;
Phases2D = zeros(Int32, Grid2D.N...);

Data2D = CartData(Grid2D, (T = Temp2D, Phases = Phases2D))

@test NumValue(Data.x[3, 1, 2]) ≈ 2.2222222222222223


# LithosphericPhases
LP = LithosphericPhases(Layers = [5 10 6], Phases = [0 1 2 3], Tlab = nothing);
X, Y, Z = xyz_grid(-5:1:5, -5:1:5, -20:1:5);
Phases = zeros(Int32, size(X));
Temp = zeros(Int32, size(X));
Phases = compute_phase(Phases, Temp, X, Y, Z, LP);

@test Phases[1, 1, end] == 3
@test Phases[1, 1, 7] == 1

Phases = compute_phase(Phases, Temp, X, Y, Z, LP, Ztop = 5);

@test Phases[1, 1, end - 4] == 0
@test Phases[1, 1, 5] == 2

LP = LithosphericPhases(Layers = [0.5 1.0 1.0], Phases = [0 1 2], Tlab = nothing);
Grid = read_LaMEM_inputfile("test_files/SaltModels.dat");
Phases = zeros(Int32, size(Grid.X));
Temp = zeros(Int32, size(Grid.X));
Phases = compute_phase(Phases, Temp, Grid, LP);

@test Phases[1, 1, 25] == 1
@test Phases[1, 1, 73] == 0

# Create Grid & nondimensionalize it
CharDim = GEO_units();
Grid = create_CartGrid(size = (10, 20, 30), x = (0.0km, 10km), y = (0.0km, 10km), z = (-10.0km, 2.0km), CharDim = CharDim)
@test Grid.Δ[2] ≈ 0.0005263157894736842


# test 1D-explicit thermal solver for AddBox -----------
nel = 96
Grid = create_CartGrid(size = (nel, nel, nel), x = (-200.0, 200.0), y = (-200.0, 200.0), z = (-200.0, 0))
Temp = zeros(Float64, Grid.N...);
Phases = zeros(Int64, Grid.N...);

# 1) horizontally layer lithosphere; UpperCrust,LowerCrust,Mantle
add_box!(
    Phases, Temp, Grid, xlim = (-100, 100), zlim = (-100, 0), Origin = (0.0, 0.0, 0.0),
    phase = LithosphericPhases(Layers = [20 15 65], Phases = [1 2 3], Tlab = nothing),
    DipAngle = 0.0, T = LithosphericTemp(nz = 201)
)

@test sum(Temp[Int64(nel / 2), Int64(nel / 2), :]) ≈ 36141.716490778184

# 2) inclined lithosphere; UpperCrust,LowerCrust,Mantle
Temp = zeros(Float64, Grid.N...);
Phases = zeros(Int64, Grid.N...);

add_box!(
    Phases, Temp, Grid, xlim = (-100, 100), zlim = (-100, 0), Origin = (0.0, 0.0, 0.0),
    phase = LithosphericPhases(Layers = [20 15 65], Phases = [1 2 3], Tlab = nothing),
    DipAngle = 30.0, T = LithosphericTemp(nz = 201)
)

@test sum(Temp[Int64(nel / 2), Int64(nel / 2), :]) ≈ 41890.41004596295

# 3) inclined lithosphere with respect to the default origin; UpperCrust,LowerCrust,Mantle
Temp = zeros(Float64, Grid.N...);
Phases = zeros(Int64, Grid.N...);

add_box!(
    Phases, Temp, Grid, xlim = (-100, 100), zlim = (-100, 0),
    phase = LithosphericPhases(Layers = [20 15 65], Phases = [1 2 3], Tlab = nothing),
    DipAngle = 30.0, T = LithosphericTemp(nz = 201)
)

@test sum(Temp[Int64(nel / 2), Int64(nel / 2), :]) ≈ 41308.02652464832

# 4) inclined lithosphere with only two layers
Temp = zeros(Float64, Grid.N...);
Phases = zeros(Int64, Grid.N...);

ρM = 3.0e3            # Density [ kg/m^3 ]
CpM = 1.0e3           # Specific heat capacity [ J/kg/K ]
kM = 2.3              # Thermal conductivity [ W/m/K ]
HM = 0.0              # Radiogenic heat source per mass [H] = W/kg; [H] = [Q/rho]
ρUC = 2.7e3           # Density [ kg/m^3 ]
CpUC = 1.0e3          # Specific heat capacity [ J/kg/K ]
kUC = 3.0             # Thermal conductivity [ W/m/K ]
HUC = 617.0e-12       # Radiogenic heat source per mass [H] = W/kg; [H] = [Q/rho]

rheology = (
    # Name              = "UpperCrust",
    SetMaterialParams(;
        Phase = 1,
        Density = ConstantDensity(; ρ = ρUC),
        HeatCapacity = ConstantHeatCapacity(; Cp = CpUC),
        Conductivity = ConstantConductivity(; k = kUC),
        RadioactiveHeat = ConstantRadioactiveHeat(; H_r = HUC * ρUC),     # [H] = W/m^3
    ),
    # Name              = "LithosphericMantle",
    SetMaterialParams(;
        Phase = 2,
        Density = ConstantDensity(; ρ = ρM),
        HeatCapacity = ConstantHeatCapacity(; Cp = CpM),
        Conductivity = ConstantConductivity(; k = kM),
        RadioactiveHeat = ConstantRadioactiveHeat(; H_r = HM * ρM),       # [H] = W/m^3
    ),
);

add_box!(
    Phases, Temp, Grid, xlim = (-100, 100), zlim = (-100, 0),
    phase = LithosphericPhases(Layers = [20 80], Phases = [1 2], Tlab = nothing),
    DipAngle = 30.0, T = LithosphericTemp(rheology = rheology, nz = 201)
)

@test sum(Temp[Int64(nel / 2), Int64(nel / 2), :]) ≈ 40506.25905747755

# 5) using flux lower boundary conditions
Temp = zeros(Float64, Grid.N...);
Phases = zeros(Int64, Grid.N...);

add_box!(
    Phases, Temp, Grid, xlim = (-100, 100), zlim = (-100, 0),
    phase = LithosphericPhases(Layers = [20 15 65], Phases = [1 2 3], Tlab = nothing),
    DipAngle = 30.0, T = LithosphericTemp(lbound = "flux", nz = 201)
)

@test sum(Temp[Int64(nel / 2), Int64(nel / 2), :]) ≈ 37366.43817255413


# Test the McKenzie thermal structure

# Create CartGrid struct
x = LinRange(0.0, 1200.0, 64);
y = LinRange(0.0, 1200.0, 64);
z = LinRange(-660, 50, 64);
Cart = CartData(xyz_grid(x, y, z));

# initialize phase and temperature matrix
Phase = ones(Int32, size(Cart));
Temp = ones(Float64, size(Cart)) * 1350;

# Create thermal structures
TsHC = HalfspaceCoolingTemp(Tsurface = 20.0, Tmantle = 1350, Age = 120, Adiabat = 0.4)
TsMK = McKenzie_subducting_slab(Tsurface = 20.0, Tmantle = 1350.0, v_cm_yr = 4.0, Adiabat = 0.0)

@test TsMK.v_cm_yr == 4.0
@test TsMK.it == 36

# Add a box with a McKenzie thermal structure

# horizontal
Temp = ones(Float64, size(Cart)) * 1350;
add_box!(Phase, Temp, Cart; xlim = (0.0, 600.0), ylim = (0.0, 600.0), zlim = (-80.0, 0.0), phase = ConstantPhase(5), T = TsMK);
@test sum(Temp) ≈ 3.509714275430339e8

# inclined slab
Temp = ones(Float64, size(Cart)) * 1350;
add_box!(Phase, Temp, Cart; xlim = (0.0, 600.0), ylim = (0.0, 600.0), zlim = (-80.0, 0.0), StrikeAngle = 0, DipAngle = 45, phase = ConstantPhase(5), T = TsMK);
@test sum(Temp) ≈ 3.511828440428276e8

# Origin must be [xmin,ymin,zmax] of the box for a McKenzie thermal structure
Temp = ones(Float64, size(Cart)) * 1350;
add_box!(Phase, Temp, Cart; xlim = (0.0, 600.0), ylim = (0.0, 600.0), zlim = (-80.0, 0.0), Origin = (0.0, 0.0, 0.0), phase = ConstantPhase(5), T = TsMK);
@test sum(Temp) ≈ 3.509714275430339e8
@test_throws "Origin is not set up correctly for McKenzie_subducting_slab" add_box!(Phase, Temp, Cart; xlim = (0.0, 600.0), ylim = (0.0, 600.0), zlim = (-80.0, 0.0), Origin = (10.0, 0.0, 0.0), phase = ConstantPhase(5), T = TsMK)


# horizontal slab, constant T
T_slab = LinearWeightedTemperature(0, 1, 600.0, :X, ConstantTemp(1000), ConstantTemp(2000));
Temp = ones(Float64, size(Cart)) * 1350;
add_box!(Phase, Temp, Cart; xlim = (0.0, 600.0), ylim = (0.0, 600.0), zlim = (-80.0, 0.0), phase = ConstantPhase(5), T = T_slab);
@test sum(Temp) ≈ 3.549127111111111e8

# horizontal slab, halfspace and McKenzie
T_slab = LinearWeightedTemperature(crit_dist = 600, F1 = TsHC, F2 = TsMK);
Temp = ones(Float64, size(Cart)) * 1350;
add_box!(Phase, Temp, Cart; xlim = (0.0, 600.0), ylim = (0.0, 600.0), zlim = (-80.0, 0.0), phase = ConstantPhase(5), T = T_slab);
@test sum(Temp) ≈ 3.495550378910777e8


Data_Final = addfield(Cart, "Temp", Temp)


# test polygon structure

x = LinRange(0.0, 1200.0, 64);
y = LinRange(0.0, 1200.0, 64);
z = LinRange(-660, 50, 64);
Cart = CartData(xyz_grid(x, y, z));

# initialize phase and temperature matrix
Phase = ones(Int32, (length(x), length(y), length(z)));
Temp = ones(Float64, (length(x), length(y), length(z))) * 1350;

add_box!(Phase, Temp, Cart; xlim = (0.0, 600.0), ylim = (0.0, 600.0), zlim = (-80.0, 0.0), phase = ConstantPhase(5), T = T = ConstantTemp(120.0));

# add accretionary prism
add_polygon!(Phase, Temp, Cart; xlim = (500.0, 200.0, 500.0), ylim = (100.0, 400.0), zlim = (-5.0, -5.0, -65.0), phase = ConstantPhase(8), T = LinearTemp(Ttop = 20, Tbot = 30))
@test maximum(Phase) == 8
@test minimum(Temp) ≈ 20.22486772486773
@test sum(Phase) == 293264

# vector limits give the same polygon as tuples
Phase_v = ones(Int32, size(Phase));
Temp_v = fill(1350.0, size(Temp));
add_box!(Phase_v, Temp_v, Cart; xlim = (0.0, 600.0), ylim = (0.0, 600.0), zlim = (-80.0, 0.0), phase = ConstantPhase(5), T = ConstantTemp(120.0));
add_polygon!(Phase_v, Temp_v, Cart; xlim = [500.0, 200.0, 500.0], ylim = [100.0, 400.0], zlim = [-5.0, -5.0, -65.0], phase = ConstantPhase(8), T = LinearTemp(Ttop = 20, Tbot = 30))
@test Phase_v == Phase
@test Temp_v == Temp

# Test the Bending slab geometry

# Create CartGrid struct
x = LinRange(0.0, 1200.0, 128);
y = LinRange(0.0, 1200.0, 128);
z = LinRange(-660, 50, 128);
Cart = CartData(xyz_grid(x, y, z));
X, Y, Z = xyz_grid(x, y, z);

# initialize phase and temperature matrix
Phase = ones(Int32, size(Cart));
Temp = fill(1350.0, size(Cart));

t1 = Trench(Start = (400.0, 400.0), End = (800.0, 800.0), θ_max = 45, direction = 1.0, n_seg = 50, Length = 600.0, Thickness = 80.0, Lb = 500.0, d_decoupling = 100.0, type_bending = :Ribe)
@test t1.θ_max == 45.0
@test t1.Thickness == 80.0
@test t1.Length == 600.0
@test t1.Lb == 500.0

phase = LithosphericPhases(Layers = [5 7 88], Phases = [2 3 4], Tlab = nothing)
TsHC = HalfspaceCoolingTemp(Tsurface = 20.0, Tmantle = 1350, Age = 30, Adiabat = 0.4)
temp = TsHC;

add_slab!(Phase, Temp, Cart, t1, phase = phase, T = TsHC)
@test Temp[84, 84, 110] ≈ 1044.9174678275472
@test extrema(Phase) == (1, 4)

# with weak zone
t1 = Trench(Start = (400.0, 400.0), End = (800.0, 800.0), θ_max = 45, direction = 1.0, n_seg = 50, Length = 600.0, Thickness = 80.0, Lb = 500.0, d_decoupling = 100.0, WeakzoneThickness = 10, WeakzonePhase = 9)
Phase = ones(Int32, size(Cart));
Temp = fill(1350.0, size(Cart));
add_slab!(Phase, Temp, Cart, t1, phase = phase, T = TsHC)
@test extrema(Phase) == (1, 9)

#Data_Final      =   CartData(X,Y,Z,(Phase=Phase,Temp=Temp))
#write_paraview(Data_Final, "Data_Final");

Phase = ones(Int32, size(Cart));
Temp = fill(1350.0, size(Cart));
TsMK = McKenzie_subducting_slab(Tsurface = 20.0, Tmantle = 1350.0, v_cm_yr = 4.0, Adiabat = 0.0)
temp = TsMK

Phase = ones(Int32, size(Cart));
Temp = fill(1350.0, size(Cart));
TsHC = HalfspaceCoolingTemp(Tsurface = 20.0, Tmantle = 1350, Age = 120, Adiabat = 0.4)
TsMK = McKenzie_subducting_slab(Tsurface = 20.0, Tmantle = 1350.0, v_cm_yr = 4.0, Adiabat = 0.0)
T_slab = LinearWeightedTemperature(crit_dist = 600, F1 = TsHC, F2 = TsMK);
phase = LithosphericPhases(Layers = [5 7 88], Phases = [2 3 4], Tlab = nothing)

t1 = Trench(Start = (400.0, 400.0), End = (800.0, 800.0), θ_max = 90.0, direction = 1.0, n_seg = 50, Length = 600.0, Thickness = 80.0, Lb = 500.0, d_decoupling = 100.0, type_bending = :Ribe, WeakzoneThickness = 10, WeakzonePhase = 9)

add_slab!(Phase, Temp, Cart, t1, phase = phase, T = T_slab)
@test Temp[84, 84, 110] ≈ 624.4897312106756

Data_Final = CartData(X, Y, Z, (Phase = Phase, Temp = Temp))

# 2D slab:
nx, nz = 512, 128
x = range(-1000, 1000, nx);
z = range(-660, 0, nz);
Grid2D = CartData(xyz_grid(x, 0, z))
Phases = zeros(Int64, nx, 1, nz);
Temp = fill(1350.0, nx, 1, nz);
add_box!(Phases, Temp, Grid2D; xlim = (-800.0, 0.0), zlim = (-80.0, 0.0), phase = ConstantPhase(1), T = HalfspaceCoolingTemp(Age = 40));

trench = Trench(Start = (0.0, -100.0), End = (0.0, 100.0), Thickness = 80.0, θ_max = 30.0, Length = 300, Lb = 150);
add_slab!(Phases, Temp, Grid2D, trench, phase = ConstantPhase(2), T = HalfspaceCoolingTemp(Age = 40));

T_slab = LinearWeightedTemperature(F1 = HalfspaceCoolingTemp(Age = 40), F2 = McKenzie_subducting_slab(Tsurface = 0, v_cm_yr = 4, Adiabat = 0.0), crit_dist = 600)
add_slab!(Phases, Temp, Grid2D, trench, phase = ConstantPhase(2), T = T_slab);

@test sum(Temp) ≈ 8.570245073330998e7
@test extrema(Phases) == (0, 2)

# Add them to the `CartData` dataset:
#Grid2D = CartData(Grid2D.x.val, Grid2D.y.val, Grid2D.z.val ,(;Phases, Temp))
#write_paraview(Grid2D,"Grid2D_SubductionCurvedMechanical");


# More sophisticated 2D example with overriding plate
nx, nz = 512, 128
x = range(-1000, 1000, nx);
z = range(-660, 0, nz);
Grid2D = CartData(xyz_grid(x, 0, z))
Phases = zeros(Int64, nx, 1, nz);
Temp = fill(1350.0, nx, 1, nz);
lith = LithosphericPhases(Layers = [15 20 55], Phases = [3 4 5], Tlab = 1250)

# Lets add the overriding plate. Note that we add this twice with a different thickness to properly represent the transition around the trench
add_box!(Phases, Temp, Grid2D; xlim = (200.0, 1000.0), zlim = (-150.0, 0.0), phase = lith, T = HalfspaceCoolingTemp(Age = 80));
add_box!(Phases, Temp, Grid2D; xlim = (0.0, 200.0), zlim = (-60.0, 0.0), phase = lith, T = HalfspaceCoolingTemp(Age = 80));

# The horizontal part of the oceanic plate is as before
v_spread_cm_yr = 3      #spreading velocity
lith = LithosphericPhases(Layers = [15 55], Phases = [1 2], Tlab = 1250)
add_box!(Phases, Temp, Grid2D; xlim = (-800.0, 0.0), zlim = (-150.0, 0.0), phase = lith, T = SpreadingRateTemp(SpreadingVel = v_spread_cm_yr));

# Yet, now we add a trench as well.
AgeTrench_Myrs = 800.0e3 / (v_spread_cm_yr / 1.0e2) / 1.0e6    #plate age @ trench

# We want to add a smooth transition from a halfspace cooling thermal profile to a slab that is heated by the surrounding mantle below a decoupling depth `d_decoupling`.
T_slab = LinearWeightedTemperature(F1 = HalfspaceCoolingTemp(Age = AgeTrench_Myrs), F2 = McKenzie_subducting_slab(Tsurface = 0, v_cm_yr = v_spread_cm_yr, Adiabat = 0.0), crit_dist = 600)

# # in this case, we have a more reasonable slab thickness:
trench = Trench(
    Start = (0.0, -100.0), End = (0.0, 100.0), Thickness = 90.0, θ_max = 30.0, Length = 600, Lb = 200,
    WeakzoneThickness = 15, WeakzonePhase = 6, d_decoupling = 125
);
add_slab!(Phases, Temp, Grid2D, trench, phase = lith, T = T_slab);

# Lithosphere-asthenosphere boundary:
ind = findall(Temp .> 1250 .&& (Phases .== 2 .|| Phases .== 5));
Phases[ind] .= 0;

@test sum(Temp) ≈ 8.289431571656203e7
@test extrema(Phases) == (0, 6)
#Grid2D = CartData(Grid2D.x.val,Grid2D.y.val,Grid2D.z.val, (;Phases, Temp))
#write_paraview(Grid2D,"Grid2D_SubductionCurvedOverriding");

# 2D volcano
nx, nz = 512, 128
x = range(-100.0e0, 100.0e0, nx);
z = range(-60.0e0, 5.0e0, nz);
Grid2D = CartData(xyz_grid(x, 0, z))
Phases = zeros(Int64, nx, 1, nz);
Temp = fill(1350.0, nx, 1, nz);
lith = LithosphericPhases(Layers = [15 20 55], Phases = [3 4 5], Tlab = 1250)

add_box!(Phases, Temp, Grid2D; xlim = (-100.0, 100.0), zlim = (-60.0e0, 0.0), phase = lith, T = HalfspaceCoolingTemp(Age = 80));

add_volcano!(
    Phases, Temp, Grid2D;
    volcanic_phase = 1,
    center = (mean(Grid2D.x.val), 0.0),
    height = 3,
    radius = 5,
    base = 0.0,
    background = nothing,
    T = HalfspaceCoolingTemp(Age = 20)
)

@test any(Phases[256, 1, :] .== 1) == true

# with a background topography, 2D Phase and Temp arrays give the same result as 3D ones
volcano!(P, T) = add_volcano!(P, T, Grid2D; volcanic_phase = 1, center = (mean(Grid2D.x.val), 0.0), height = 3, radius = 5, base = 0.0, background = fill(-1.0, nx, 1), T = HalfspaceCoolingTemp(Age = 20))
Phases2D, Temp2D = Phases[:, 1, :], Temp[:, 1, :]
volcano!(Phases, Temp)
volcano!(Phases2D, Temp2D)
@test Phases2D == Phases[:, 1, :]
@test Temp2D == Temp[:, 1, :]

# 3D volcano
# Create CartGrid struct
x = LinRange(0.0, 100.0, 64);
y = LinRange(0.0, 100.0, 64);
z = LinRange(-60, 5.0e0, 64);
Cart = CartData(xyz_grid(x, y, z));

# initialize phase and temperature matrix
Phase = zeros(Int32, size(Cart));
Temp = fill(1350.0, size(Cart));
lith = LithosphericPhases(Layers = [15 20 55], Phases = [3 4 5], Tlab = 1250)

add_box!(Phase, Temp, Cart; xlim = (0.0, 100.0), ylim = (0.0, 100.0), zlim = (-60.0, 0.0), phase = lith, T = HalfspaceCoolingTemp(Age = 80));

add_volcano!(
    Phase, Temp, Cart;
    volcanic_phase = 1,
    center = (mean(Cart.x.val), mean(Cart.y.val), 0.0),
    height = 3,
    radius = 5,
    base = 0.0,
    background = nothing,
    T = HalfspaceCoolingTemp(Age = 20)
)

@test any(Phase[32, 32, :] .== 1) == true

#3D fault
# Create CartGrid struct
x = LinRange(0.0, 100.0, 64);
y = LinRange(0.0, 100.0, 64);
z = LinRange(-60, 5.0e0, 64);
Cart = CartData(xyz_grid(x, y, z));

# initialize phase and temperature matrix
Phase = zeros(Int32, size(Cart));
Temp = fill(1350.0, size(Cart));
lith = LithosphericPhases(Layers = [15 20 55], Phases = [3 4 5], Tlab = 1250)

add_box!(Phase, Temp, Cart; xlim = (0.0, 100.0), ylim = (0.0, 100.0), zlim = (-60.0, 0.0), phase = lith, T = HalfspaceCoolingTemp(Age = 80));

add_fault!(
    Phase, Temp, Cart;
    Start = (0.0, 0.0), End = (100, 100),
    Fault_thickness = 1.0,
    Depth_extent = (-30.0, 0.0),
    DipAngle = -10.0e0,
    phase = ConstantPhase(1),
    T = ConstantTemp(1200),
)

@test any(Phase[32, 32, :] .== 1) == true
@test any(Temp[32, 32, :] .== 1200) == true

# Q1Data
Grid = Q1Data(xyz_grid(1.0:1:10.0, 11.0:1:20.0, -20:1:-10))
PhasesC = zeros(Int64, size(Grid));  # at cell
TempC = ones(Float64, size(Grid)) * 1350;
PhasesV = zeros(Int64, size(Grid.x));  # at vertex
TempV = ones(Float64, size(Grid.x)) * 1350;

# Add data to vertex fields:
add_box!(PhasesV, TempV, Grid, xlim = (2, 4), zlim = (-15, -10), phase = ConstantPhase(3), DipAngle = 10, T = LinearTemp(Tbot = 1350, Ttop = 200))
@test sum(TempV[1, 1, :]) ≈ 14850.0

add_ellipsoid!(PhasesV, TempV, Grid, cen = (4, 15, -17), axes = (1, 2, 3), StrikeAngle = 90, DipAngle = 45, phase = ConstantPhase(2), T = ConstantTemp(600))
@test sum(TempV[1, 1, :]) ≈ 14850.0


# Add data to cell fields:
add_box!(PhasesC, TempC, Grid, xlim = (2, 4), zlim = (-15, -10), phase = ConstantPhase(3), DipAngle = 10, T = LinearTemp(Tbot = 1350, Ttop = 200), cell = true)
@test sum(TempC[1, 1, :]) ≈ 13051.985346405856

add_ellipsoid!(PhasesC, TempC, Grid, cen = (4, 15, -17), axes = (1, 2, 3), StrikeAngle = 90, DipAngle = 45, phase = ConstantPhase(2), T = ConstantTemp(1600), cell = true)
@test all(extrema(TempC) .≈ (253.34427030131295, 1600.0))


# ---------------------------------------------------------------------------
# Thermal structures must be anchored at the *geometric* top of the region, not
# at the topmost grid point inside it (which depends on the grid resolution).
# See PR #206 (Ttop/Tsurface wrongly assumed at z = 0) and its follow-up fix.
erfc = GeophysicalModelGenerator.erfc
halfspace_T(depth_km, Tsurf, Tmantle, Age_Myr) = (Tsurf - Tmantle) * erfc(depth_km * 1.0e3 / (2 * sqrt(1.0e-6 * Age_Myr * 1.0e6 * 3600 * 24 * 365.25))) + Tmantle

# 1) Box with its top at z = 0, on a grid WITH a node at z = 0 and on one WITHOUT:
#    the temperature at the same physical depth must be the analytic halfspace value in both cases
x = -500.0:100.0:500.0; y = -10.0:20.0:10.0
for z in (-1000.0:20.0:50.0, range(-1000.0, 50.0, length = 65))     # 2nd grid: 16.4 km spacing, no node at 0
    Cart = CartData(xyz_grid(x, y, z))
    Phase = zeros(Int64, size(Cart.x)); Temp = zeros(Float64, size(Cart.x))
    add_box!(Phase, Temp, Cart; xlim = (-500.0, 500.0), zlim = (-1000.0, 0.0), phase = ConstantPhase(3), T = HalfspaceCoolingTemp(Tsurface = 0, Tmantle = 1350, Age = 100))
    k = findlast(z .<= 0.0)                                             # topmost node inside the box
    @test Temp[6, 1, k] ≈ halfspace_T(-z[k], 0, 1350, 100)
    @test Temp[6, 1, k - 5] ≈ halfspace_T(-z[k - 5], 0, 1350, 100)
end

# 2) Box whose top is at depth (the case PR #206 addressed): Tsurface applies at the box top, not at z = 0
z = range(-1000.0, 50.0, length = 65)
Cart = CartData(xyz_grid(x, y, z))
Phase = zeros(Int64, size(Cart.x)); Temp = zeros(Float64, size(Cart.x))
add_box!(Phase, Temp, Cart; xlim = (-500.0, 500.0), zlim = (-1000.0, -10.0), phase = ConstantPhase(3), T = HalfspaceCoolingTemp(Tsurface = 0, Tmantle = 1350, Age = 100))
k = findlast(z .<= -10.0)
@test Temp[6, 1, k] ≈ halfspace_T(-10.0 - z[k], 0, 1350, 100)

# 3) Polygon with its top at depth and a linear gradient: Ttop at the polygon top, Tbot at its bottom
Phase = zeros(Int64, size(Cart.x)); Temp = zeros(Float64, size(Cart.x))
add_polygon!(Phase, Temp, Cart; xlim = (-400.0, 400.0, 400.0, -400.0), ylim = (-10.0, 10.0), zlim = (-10.0, -10.0, -210.0, -210.0), phase = ConstantPhase(2), T = LinearTemp(Ttop = 20, Tbot = 420))
k = findlast(z .<= -10.0)
@test Temp[6, 1, k] ≈ 20 + (-10.0 - z[k]) / 200.0 * 400        # 2 C/km below the polygon top

# add_box!/add_plate! with ridge segments use absolute coordinates; bounds wrappers forward all keywords
let
    Grid = CartData(xyz_grid(-1000:100.0:1000, -1000:50.0:1000, -100:10.0:0))
    segs = [((-500.0, -1000.0), (-500.0, 0.0)), ((-250.0, 0.0), (-250.0, 200.0)), ((-750.0, 200.0), (-750.0, 1000.0))]
    Tr = SpreadingRateTemp(SpreadingVel = 3)
    ix, iy, iz = findfirst(==(-500.0), -1000:100.0:1000), findfirst(==(-250.0), -1000:50.0:1000), findfirst(==(-10.0), -100:10.0:0)
    Ph = fill(2, size(Grid)); Tb = fill(1350.0, size(Grid))
    add_box!(Ph, Tb, Grid; xlim = (-1000.0, 0.0), ylim = (-500.0, 500.0), zlim = (-80.0, 0.0), T = Tr, segments = segs)
    Ph = fill(2, size(Grid)); Tp = fill(1350.0, size(Grid))
    add_plate!(Ph, Tp, Grid, [[-1000.0, 0.0, 0.0, -1000.0], [-500.0, -500.0, 500.0, 500.0], [-80.0, 0.0]]; T = Tr, segments = segs)
    @test Tb[ix, iy, iz] ≈ 1350.0 atol = 1.0e-3     # on the ridge the age is ~0: mantle temperature
    @test Tb[ix + 3, iy, iz] < 1000.0               # off the ridge the lithosphere is cold
    @test Tb[ix + 3, iy, iz] ≈ Tp[ix + 3, iy, iz]

    Ph = fill(2, size(Grid)); Tw = fill(1350.0, size(Grid))
    add_box!(Ph, Tw, Grid, [[-1000.0, 0.0], [-500.0, 500.0], [-80.0, 0.0]]; T = Tr, segments = segs)
    @test Tw == Tb
end

# add_layer! accepts all three limits
let
    Grid = CartData(xyz_grid(0:10.0, 0:10.0, -10:1.0:0))
    Ph = zeros(Int, size(Grid)); Tm = zeros(size(Grid))
    add_layer!(Ph, Tm, Grid; xlim = (0.0, 5.0), ylim = (0.0, 5.0), zlim = (-5.0, 0.0), phase = ConstantPhase(1))
    @test count(==(1), Ph) == 6 * 6 * 6
    @test_throws "at least one" add_layer!(Ph, Tm, Grid)
end

# Slab construction: the whole slab is filled and both Start/End orientations produce a slab
let
    @test size(compute_slab_surface(Trench(Length = 300.0, n_seg = 7, Thickness = 50.0, Lb = 100.0))[1], 1) == 8
    @test size(compute_slab_surface(Trench(Length = 1.0, n_seg = 10, Thickness = 50.0, Lb = 1.0))[1], 1) == 11
    Grid = CartData(xyz_grid(0:20.0:1200, 0:20.0:1200, -400:20.0:0))
    count_slab(t) = (Ph = zeros(Int, size(Grid)); Tm = zeros(size(Grid)); add_slab!(Ph, Tm, Grid, t); count(!=(0), Ph))
    kw = (θ_max = 45.0, Length = 400.0, Thickness = 80.0, Lb = 300.0, n_seg = 2)
    @test count_slab(Trench(; Start = (400.0, 400.0), End = (800.0, 800.0), kw...)) > 0
    n_fwd = count_slab(Trench(; Start = (800.0, 400.0), End = (400.0, 0.0), kw...))
    n_rev = count_slab(Trench(; Start = (400.0, 0.0), End = (800.0, 400.0), kw...))
    @test n_fwd > 0
    @test n_rev > 0
    @test count_slab(Trench(; Start = (400.0, 400.0), End = (800.0, 800.0), kw..., n_seg = 1)) > 0
end

# add_fault!: the surface trace stays at Start/End and a dipping fault shifts in the dip direction with depth
let
    Grid = CartData(xyz_grid(0:10.0:1000, 0:10.0:200, -100:10.0:0))
    x = 0:10.0:1000
    Ph = zeros(Int, size(Grid)); Tm = zeros(size(Grid))
    add_fault!(Ph, Tm, Grid; Start = (500.0, 0.0), End = (500.0, 100.0), Fault_thickness = 10.0, DipAngle = 60.0)
    @test Ph[51, 1, end] != 0                      # trace passes through x = 500 at the surface
    @test all(485 .<= x[findall(!=(0), Ph[:, 1, end])] .<= 515)
    xdeep = x[findall(!=(0), Ph[:, 1, 1])]
    @test all(>(500), xdeep)                       # z = -100: offset of 100*tan(60) ~ 173 km in +x
    @test mean(xdeep) ≈ 500 + 100 * tand(60) atol = 10

    Ph = zeros(Int, size(Grid)); Tm = zeros(size(Grid))
    add_fault!(Ph, Tm, Grid; Start = (0.0, 100.0), End = (1000.0, 100.0), Fault_thickness = 10.0, DipAngle = 60.0)
    @test all(Ph[:, 1, end] .== 0)                 # E-W fault dips in y: surface trace is at y = 100
    @test count(!=(0), Ph[:, 11, end]) > 0
    @test count(!=(0), Ph[:, 11, 1]) == 0
end

# add_layer! with curved bounds: the bounds are functions z(x, y), and depth is measured from a curved top
let
    x, z = 0:1.0:100, -40:1.0:10
    Grid = CartData(xyz_grid(x, 0, z))
    bottom(x, y) = -20 + 5cospi(x / 50)
    Ph = zeros(Int, size(Grid)); Tm = zeros(size(Grid))
    add_layer!(Ph, Tm, Grid; zlim = (bottom, 0.0), phase = ConstantPhase(1))
    @test all(Ph[i, 1, k] == (bottom(x[i], 0) <= z[k] <= 0) for i in eachindex(x), k in eachindex(z))

    # crust 10 km thick below a sloping top that rises above z = 0; halfspace cooling starts at the top
    top(x, y) = 2 - 0.2x
    lith = LithosphericPhases(Layers = [10, 50], Phases = [1, 2, 3])
    Ph = zeros(Int, size(Grid)); Tm = fill(1350.0, size(Grid))
    add_layer!(Ph, Tm, Grid; zlim = (-40.0, top), phase = lith, T = HalfspaceCoolingTemp(Tsurface = 20))
    for i in (1, 51, 101)
        k = findlast(<=(top(x[i], 0)), z)
        @test Ph[i, 1, k] == 1 && Ph[i, 1, k - 9] == 1 && Ph[i, 1, k - 11] == 2
        @test Tm[i, 1, k] ≈ 20
        @test Ph[i, 1, (k + 1):end] == zeros(Int, length(z) - k)
    end
    @test Tm[1, 1, 21] ≈ Tm[51, 1, 11]   # 22 km below the top in both columns
end

# add_fault! with DetachmentDepth: a listric fault flattens into a detachment
let
    x, z = 0:0.5:200, -30:1.0:0
    Grid = CartData(xyz_grid(x, 0:10.0:20, z))
    Ph = zeros(Int, size(Grid)); Tm = zeros(size(Grid))
    add_fault!(Ph, Tm, Grid; Start = (50.0, 0.0), End = (50.0, 20.0), Fault_thickness = 2.0, DipAngle = 30.0, DetachmentDepth = 15.0, phase = ConstantPhase(1))
    xfault(k) = x[findall(==(1), Ph[:, 2, k])]
    @test mean(xfault(length(z))) ≈ 50 atol = 0.5                          # trace at the surface
    @test mean(xfault(findfirst(==(-5), z))) ≈ 50 - 15log(2 / 3) * tand(30) atol = 0.5
    @test length(xfault(findfirst(==(-15), z))) > 100                     # horizontal detachment
    @test isempty(xfault(findfirst(==(-20), z)))
    # with a very deep detachment, the fault is planar
    P1, P2 = zeros(Int, size(Grid)), zeros(Int, size(Grid))
    add_fault!(P1, Tm, Grid; Start = (50.0, 0.0), End = (50.0, 20.0), Fault_thickness = 2.1, DipAngle = -30.0)
    add_fault!(P2, Tm, Grid; Start = (50.0, 0.0), End = (50.0, 20.0), Fault_thickness = 2.1, DipAngle = -30.0, DetachmentDepth = 1.0e9)
    @test P1 == P2
    @test_throws "DetachmentDepth must be positive" add_fault!(Ph, Tm, Grid; DetachmentDepth = -1)
    @test_throws "nonzero DipAngle" add_fault!(Ph, Tm, Grid; DetachmentDepth = 10)
end

# add_fault! with a curved Trace: the fault follows the arc and dips perpendicular to it
let
    x, y, z = 0:1.0:100, 0:1.0:100, -20:1.0:0
    Grid = CartData(xyz_grid(x, y, z))
    X, Y, Z = Grid.x.val, Grid.y.val, Grid.z.val
    r, a = hypot.(X .- 50, Y .- 50), atand.(Y .- 50, X .- 50)
    arc = [(50 + 40cosd(α), 50 + 40sind(α)) for α in 0:5:180]   # counterclockwise: the right side is outside
    Ph = zeros(Int, size(Grid)); Tm = zeros(size(Grid))
    add_fault!(Ph, Tm, Grid; Trace = arc, Fault_thickness = 2.0, DipAngle = -30.0, phase = ConstantPhase(1))
    surface, deep = (Ph .== 1) .& (Z .== 0), (Ph .== 1) .& (Z .== -20)
    @test all(abs.(r[surface] .- 40) .<= 1 / cosd(30) + 40(1 - cosd(2.5)))   # half-width of the dipping fault plus the sagitta of the chords
    @test all(0 .<= a[surface] .<= 180)                          # the fault ends at the ends of the trace
    @test mean(r[deep]) ≈ 40 - 20tand(30) atol = 0.5             # dips towards the center of the arc
    @test count(surface) > 200

    # a trace of two points equals Start and End
    P1, P2 = zeros(Int, size(Grid)), zeros(Int, size(Grid))
    add_fault!(P1, Tm, Grid; Start = (10.0, 20.0), End = (80.0, 70.0), DipAngle = 20.0, DetachmentDepth = 15.0)
    add_fault!(P2, Tm, Grid; Trace = [(10.0, 20.0), (80.0, 70.0)], DipAngle = 20.0, DetachmentDepth = 15.0)
    @test P1 == P2 && count(==(1), P1) > 0
    @test_throws "at least two points" add_fault!(Ph, Tm, Grid; Trace = [(1.0, 2.0)])
    @test_throws "must differ" add_fault!(Ph, Tm, Grid; Trace = [(1.0, 2.0), (1.0, 2.0), (3.0, 4.0)])
end

# add_slab! with a curved trench: the slab follows the arc and dips perpendicular to it
let
    Grid = CartData(xyz_grid(0:5.0:600, 0:5.0:600, -200:5.0:0))
    X, Y, Z = Grid.x.val, Grid.y.val, Grid.z.val
    r, a = hypot.(X .- 300, Y .- 100), atand.(Y .- 100, X .- 300)
    # counterclockwise arc around (300, 100): direction = 1 bends the slab to the left of the trench, inside the arc
    arc = [(300 + 300cosd(α), 100 + 300sind(α)) for α in 20:2:160]
    kw = (Thickness = 50.0, θ_max = 45.0, Length = 200.0, Lb = 100.0, direction = 1.0)
    Ph = zeros(Int, size(Grid)); Tm = zeros(size(Grid))
    add_slab!(Ph, Tm, Grid, Trench(; Trace = arc, kw...); phase = ConstantPhase(1))
    slab = Ph .== 1
    @test all(20 - 3 .<= a[slab] .<= 160 + 3)                 # the slab ends perpendicular to the end segments of the trench
    @test all(r[slab] .<= 300 + 1)
    @test mean(r[slab .& (Z .== -100)]) < mean(r[slab .& (Z .== -40)]) < 300
    @test occursin("71 points", sprint(show, Trench(; Trace = arc)))

    # a straight trench given as a trace places the slab like Start and End, up to rounding at its boundary
    P1, P2 = zeros(Int, size(Grid)), zeros(Int, size(Grid))
    add_slab!(P1, Tm, Grid, Trench(; Start = (100.0, 100.0), End = (500.0, 400.0), kw...))
    add_slab!(P2, Tm, Grid, Trench(; Trace = [(100.0, 100.0), (500.0, 400.0)], kw...))
    @test count(P1 .!= P2) <= 1.0e-3 * count(==(1), P1)
end

# Trench with a dip profile (Bending), a part behind the trench (Behind) and a thermal thickness (ThermalThickness)
let
    Top = compute_slab_surface(Trench(Bending = [(20, 150), (80, 100), (80, 150)], n_seg = 40, direction = 1.0))[1]
    @test -diff(Top[(end - 1):end, 2])[1] / diff(Top[(end - 1):end, 1])[1] ≈ tand(80) rtol = 1.0e-6   # 80° over the last part
    @test sum(hypot.(diff(Top[:, 1]), diff(Top[:, 2]))) ≈ 400                                      # length along the top
    Top, Bottom = compute_slab_surface(Trench(Bending = [(20, 150), (80, 250)], n_seg = 40, Thickness = 60.0))
    @test all(i -> hypot((Top[i, :] .- Bottom[i, :])...) ≈ 60, axes(Top, 1))                          # thickness

    Top = compute_slab_surface(Trench(Bending = [(30, 300)], Behind = 50.0, direction = -1.0))[1]
    @test Top[1, :] == [50.0, 0.0]                                       # starts 50 km behind the trench

    # Behind and ThermalThickness in a model: trench along x = 100, slab dipping east
    x, z = 0:2.0:300, -300:2.0:0
    Grid = CartData(xyz_grid(x, 0:20.0:100, z))
    t = Trench(
        Start = (100.0, 0.0), End = (100.0, 100.0), direction = -1.0, Thickness = 50.0, Bending = [(45, 300)],
        Behind = 40.0, ThermalThickness = 150.0
    )
    Ph = zeros(Int, size(Grid)); Tm = fill(1350.0, size(Grid))
    add_slab!(Ph, Tm, Grid, t; phase = ConstantPhase(1), T = LinearTemp(Ttop = 0, Tbot = 1300))
    @test any(==(1), Ph[x .< 100, :, :])                                 # slab behind the trench
    @test all(==(0), Ph[x .< 60, :, :])
    col = findfirst(==(80.0), x)                                         # behind the trench, the slab is horizontal
    @test count(==(1), Ph[col, 3, :]) ≈ 26 atol = 1                      # phases within Thickness
    @test count(<(1300), Tm[col, 3, :]) ≈ 76 atol = 1                    # temperature within ThermalThickness
    @test Tm[col, 3, findfirst(==(-100.0), z)] ≈ 100 / 150 * 1300 rtol = 0.01
    @test occursin("[Bending]", sprint(show, t)) && occursin("[ThermalThickness]", sprint(show, t)) && !occursin("θ_max", sprint(show, t))

    @test_throws "lengths of Bending must be positive" compute_slab_surface(Trench(Bending = [(30, 0)]))
    @test_throws "Behind must not be negative" compute_slab_surface(Trench(Behind = -1.0))
end

# phase = nothing sets the temperature only
let
    Grid = CartData(xyz_grid(0:1.0:10, 0:1.0:10, -10:1.0:0))
    Ph = ones(Int, size(Grid)); Ph[:, :, 1] .= 2; P0 = copy(Ph); Tm = zeros(size(Grid))
    add_box!(Ph, Tm, Grid; xlim = (0.0, 10.0), zlim = (-10.0, 0.0), phase = nothing, T = ConstantTemp(500))
    @test Ph == P0 && all(==(500), Tm)
    add_box!(Ph, Tm, Grid; xlim = (0.0, 10.0), zlim = (-10.0, 0.0), phase = nothing, T = LithosphericTemp())
    @test Ph == P0
end

# BezierSlab: the slab lies below a Bézier surface, perpendicular to it, and ends at its edges
let
    x, z = 0:2.0:400, -200:2.0:20
    Grid = CartData(xyz_grid(x, 0:20.0:200, z))
    # planar surface dipping 30° east from a trench along x = 100, with 100 km of plate behind it; bilinear net
    c, s = cosd(30), sind(30)
    Points = [(100.0 + h, y, zz) for y in (0.0, 200.0), (h, zz) in ((-100.0, 0.0), (0.0, 0.0))]
    plane = [(100.0 + h, y, zz) for y in (0.0, 200.0), (h, zz) in ((0.0, 0.0), (250c, -250s))]
    lith = LithosphericPhases(Layers = [10, 70], Phases = [2, 3])
    Ph = zeros(Int, size(Grid)); Tm = zeros(size(Grid))
    add_slab!(Ph, Tm, Grid, BezierSlab(Points = plane, Thickness = 80.0); phase = lith, T = HalfspaceCoolingTemp(Age = 50))
    column(xc) = Ph[findfirst(==(xc), x), 6, :]
    ztop = -100tand(30)                                                    # 100 km from the trench
    @test column(200.0)[findfirst(==(-60.0), z)] == 2                      # 2.3 km normal to the top (z = -57.7)
    @test column(200.0)[findfirst(==(-80.0), z)] == 3
    @test column(200.0)[findfirst(==(-56.0), z)] == 0                      # above the slab
    @test 2count(!=(0), column(200.0)) ≈ 80 / c atol = 4                  # thickness perpendicular to the slab
    @test all(==(0), Ph[x .< 100 - 80s - 2, :, :])                        # the slab ends perpendicular to the surface at the trench
    @test Tm[findfirst(==(200.0), x), 6, findfirst(==(-58.0), z)] < 100   # cold just below the top

    Ph = zeros(Int, size(Grid))
    add_slab!(Ph, Tm, Grid, BezierSlab(Points = Points, Thickness = 50.0))
    @test count(==(1), Ph[findfirst(==(50.0), x), 6, :]) == 26              # flat plate behind the trench, z = -50 ... 0
    @test count(==(1), Ph[findfirst(==(50.0), x), end, :]) == 26             # including the edge of the surface at y = 200

    # a thickness that decreases along the trench, to 0 at y = 200
    Ph = zeros(Int, size(Grid))
    Tm = zeros(size(Grid))
    add_slab!(Ph, Tm, Grid, BezierSlab(Points = plane, Thickness = (x, y, z) -> 80 * (1 - y / 200)); T = LinearTemp())
    nslab(j) = count(==(1), Ph[findfirst(==(200.0), x), j, :])
    @test nslab(1) > 2nslab(6) > 0 && nslab(11) == 0
    @test !any(isnan, Tm)                                                   # where the thickness reaches 0
    # thermal structures see the depth of the surface as y: a blend with dir = :Y follows it
    Tm = zeros(size(Grid))
    blend = LinearWeightedTemperature(F1 = ConstantTemp(0), F2 = ConstantTemp(1000), crit_dist = 200, dir = :Y)
    add_slab!(Ph, Tm, Grid, BezierSlab(Points = plane, Thickness = 80.0); T = blend)
    i, k = findfirst(==(200.0), x), findfirst(==(-60.0), z)                # just below the surface at 57.7 km depth
    @test Tm[i, 6, k] ≈ 1000 * (100tand(30) + 2.3sind(30) * tand(30)) / 200 rtol = 0.02
    # the thermal structure is compressed with the slab: the base of the slab is equally warm where it is thinner
    Tm = zeros(size(Grid))
    add_slab!(Ph, Tm, Grid, BezierSlab(Points = plane, Thickness = (x, y, z) -> 80 * (1 - y / 400)); T = LinearTemp(Ttop = 0, Tbot = 1000))
    base(j) = maximum(Tm[findfirst(==(200.0), x), j, :])
    @test base(1) ≈ base(6) rtol = 0.05
    @test_throws "lies above" add_slab!(Ph, Tm, Grid, BezierSlab(Points = reverse(plane, dims = 1)))
    @test_throws "at least 2 × 2" add_slab!(Ph, Tm, Grid, BezierSlab(Points = plane[:, 1:1]))
end

# Tear: a vertical slab along the trench x = 100, y ∈ [0, 400], necking at 100 km depth from y = 100, torn from y = 250 on
let
    x, y, z = 50:1.0:120, 0:10.0:400, -200:5.0:0
    Grid = CartData(xyz_grid(x, y, z))
    kw = (Start = (100.0, 0.0), End = (100.0, 400.0), direction = -1.0, Thickness = 40.0, Bending = [(90, 1.0), (90, 250.0)], n_seg = 251)
    function slab(t; T = LinearTemp(Ttop = 0, Tbot = 1000), G = Grid, T0 = -1.0)
        phases, temp = zeros(Int, size(G)), fill(T0, size(G))
        add_slab!(phases, temp, G, t; phase = ConstantPhase(1), T)
        return phases, temp
    end
    torn(Ph) = dropdims(.!any(==(1), Ph, dims = 1), dims = 1)    # (y, z) without slab along x
    tear(Start, Tip) = Trench(; kw..., Tear = Tear(; Depth = 100, Width = 40, Start, Tip))
    Ph, Tm = slab(tear(100, 250))
    gap = torn(Ph)[:, z .<= 0]
    @test 240 <= y[findfirst(any, eachrow(gap))] <= 250                     # the tip, where the slab is thinner than the grid
    # beyond the tip, the gap spans |z + Depth| ≤ Width √(log A) with A = (y - Start) / (Tip - Start), within a cell
    @test all(j -> abs(maximum(abs, z[gap[j, :]] .+ 100) - 40sqrt(log((y[j] - 100) / 150))) <= 5, findall(>=(250), y))
    P0, T0 = slab(Trench(; kw...))
    @test Ph[:, y .<= 100, :] == P0[:, y .<= 100, :] && Tm[:, y .<= 100, :] == T0[:, y .<= 100, :]   # no tear before Start
    @test Ph[:, y .> 100, :] != P0[:, y .> 100, :]
    @test all(==(-1), Tm[:, y .>= 250, z .<= 0][:, gap[y .>= 250, :]])        # the gap keeps the temperature before
    # the neck warms as it thins, T = T0 + φ (T_slab - T0): with T0 the base temperature of the slab, the temperature
    # within the neck changes between cells along z by at most (Tbot - Ttop) A max|∂φ/∂z| Δz = 1000 A √(2/e) / Width × 5
    _, Tneck = slab(tear(100, 250); T0 = 1000.0)
    k = findall(-180 .<= z .<= -20)
    for j in findall(>(100), y)
        @test maximum(abs, diff(Tneck[:, j, k], dims = 2)) <= 1000 * (y[j] - 100) / 150 * sqrt(2 / ℯ) / 40 * 5
    end
    @test occursin("[Tear]", sprint(show, tear(100, 250)))

    # the reversed tear, and map points that project onto the trench at the same distances
    Pr, _ = slab(tear(250, 100))
    @test 100 <= y[findlast(any, eachrow(torn(Pr)[:, z .<= 0]))] <= 110
    @test slab(tear((300.0, 100.0), (-50, 250))) == (Ph, Tm)

    # a corner trace: the distances run along it
    trace = [(100.0, 0.0), (100.0, 200.0), (300.0, 200.0)]
    G2 = CartData(xyz_grid(0:5.0:400, 0:5.0:400, -200:5.0:0))
    kc = (; Trace = trace, Thickness = 40.0, θ_max = 60.0, Length = 300.0, Lb = 100.0, direction = 1.0)
    t(Start, Tip) = Trench(; kc..., Tear = Tear(; Depth = 100, Width = 40, Start, Tip))
    @test slab(t(100, 300); G = G2) == slab(t((150, 100), (200, 150)); G = G2)
    @test slab(t(100, 300); G = G2) != slab(Trench(; kc...); G = G2)

    # ThermalThickness and the weak zone thin with the slab: in the gap, neither temperature nor weak zone is set
    tw = Trench(; kw..., ThermalThickness = 80.0, WeakzoneThickness = 10.0, WeakzonePhase = 3, d_decoupling = 150.0, Tear = Tear(; Depth = 100, Width = 40, Start = 100, Tip = 250))
    Pw, Tw = slab(tw)
    j, k = findfirst(==(300.0), y), findall(-115 .<= z .<= -85)
    @test all(==(0), Pw[:, j, k]) && all(==(-1), Tw[:, j, k])
    @test any(==(3), Pw[:, j, findfirst(==(-50.0), z)])

    @test_throws "Depth of a Tear must be positive" Tear(; Depth = -100, Width = 40, Start = 0, Tip = 100)
    @test_throws "Width of a Tear must be positive" Tear(; Depth = 100, Width = 0, Start = 0, Tip = 100)
    @test_throws "Start and Tip of a Tear must differ" Tear(; Depth = 100, Width = 40, Start = 100, Tip = 100)
    @test_throws "beyond the ends of the trench" slab(tear(100, (50.0, 500.0)))
    @test_throws "different distances along the trench" slab(tear(100, (0.0, 100.0)))
    @test_throws "must be (x, y)" slab(tear(100, (0.0, 100.0, 0.0)))
end

# Tear of a BezierSlab: a planar slab dipping 45° from the trench x = 100, y ∈ [0, 400]
let
    x, y, z = 50:5.0:400, 0:10.0:400, -300:5.0:20
    Grid = CartData(xyz_grid(x, y, z))
    plane = [(100.0 + h, yy, zz) for yy in (0.0, 400.0), (h, zz) in ((0.0, 0.0), (300 / sqrt(2), -300 / sqrt(2)))]
    lith = LithosphericPhases(Layers = [10, 50], Phases = [2, 3])
    function slab(tear; Thickness = 60.0, G = Grid, T = HalfspaceCoolingTemp(Age = 50), T0 = 0)
        phases, temp = zeros(Int, size(G)), fill(eltype(G.x.val)(T0), size(G))
        add_slab!(phases, temp, G, BezierSlab(; Points = plane, Thickness, Tear = tear); phase = lith, T)
        return phases, temp
    end
    tear(Start, Tip) = Tear(; Depth = 100, Width = 40, Start, Tip)
    # the slab of the plane from the closest point of the surface, at distance s down the dip and depth zd, with d
    # along the upward normal, and the tear factor φ there; with T = LinearTemp(Ttop = 0, Tbot = 1000) and the
    # temperature T0 = 1000 before, the temperature is T0 + φ (T_slab - T0), with T_slab of the distance d / φ.
    # Points within 1e-6 km of the boundary of the slab are left out.
    function expected(t)
        inside, Texp, clear = falses(size(Grid)), fill(1000.0, size(Grid)), falses(size(Grid))
        for I in eachindex(Grid.x.val)
            xx, yy, zz = Grid.x.val[I], Grid.y.val[I], Grid.z.val[I]
            s, d = (xx - 100 - zz) / sqrt(2), (xx - 100 + zz) / sqrt(2)
            zd = d / sqrt(2) - zz
            φ = isnothing(t) ? 1.0 : clamp(1 - max(0, (yy - t.Start) / (t.Tip - t.Start)) * exp(-((zd - t.Depth) / t.Width)^2), 0, 1)
            margin = min(s, 300 - s, yy, 400 - yy, d + 60φ, -d)
            inside[I], clear[I] = φ > 0 && margin >= 0, abs(margin) > 1.0e-6
            inside[I] && (Texp[I] = 1000 + φ * (1000 * -d / (60φ) - 1000))
        end
        return inside, Texp, clear
    end
    for t in (nothing, tear(100, 250), tear(250, 100))
        P, T = slab(t; T = LinearTemp(Ttop = 0, Tbot = 1000), T0 = 1000)
        inside, Texp, clear = expected(t)
        @test (P .!= 0)[clear] == inside[clear]
        @test T[inside .& clear] ≈ Texp[inside .& clear]
        @test all(T[.!inside .& clear] .== 1000)                               # the gap keeps the temperature before
    end
    Ph, Tm = slab(tear(100, 250))
    P0, T0 = slab(nothing)
    @test Ph[:, y .< 100, :] == P0[:, y .< 100, :] && Tm[:, y .< 100, :] == T0[:, y .< 100, :]
    Pm, Tmm = slab(tear((0.0, 100.0), (300, 250)))
    @test Pm == Ph && Tmm ≈ Tm
    # a thickness function is multiplied by the tear
    Pf, _ = slab(tear(100, 250); Thickness = (x, y, z) -> 60.0)
    @test Pf == Ph
    @test_throws "beyond the ends of the trench" slab(tear(100, (300.0, 450.0)))
    # Float32 grid
    G32 = CartData(xyz_grid(Float32.(x), Float32.(y), Float32.(z)))
    @test slab(tear(100, 250); G = G32)[1] == Ph
end

# bezier_curve: the curve passes through its end points, and a straight one is evenly spaced
let
    c = bezier_curve([(0.0, 0.0), (1.0, 2.0), (3.0, 2.0), (4.0, 0.0)], 5)
    @test first(c) == (0.0, 0.0) && last(c) == (4.0, 0.0)
    @test c[3] == (2.0, 1.5)                                  # (P0 + 3P1 + 3P2 + P3) / 8
    @test bezier_curve([(0.0, 0.0), (4.0, 2.0)], 3) == [(0.0, 0.0), (2.0, 1.0), (4.0, 2.0)]
    @test bezier_curve([[0.0, 0.0], [2.0, 2.0], [4.0, 0.0]], 3)[2] == [2.0, 1.0]
    @test_throws "at least two control points" bezier_curve([(0.0, 0.0)])
    @test_throws "at least two output points" bezier_curve([(0.0, 0.0), (1.0, 1.0)], 1)
end

# add_layer! with a surface CartData as bound gives the same layer as the function it samples
let
    x, z = 0:1.0:100, -40:1.0:10
    Grid = CartData(xyz_grid(x, 0:1.0:10, z))
    f(x, y) = -10 - 0.5x + 0.25y
    Surf = CartData(xyz_grid(-10:1.0:110, -10:1.0:20, 0))
    Surf = CartData(Surf.x.val, Surf.y.val, f.(Surf.x.val, Surf.y.val), (; Z = f.(Surf.x.val, Surf.y.val)))
    P1, P2, Tm = zeros(Int, size(Grid)), zeros(Int, size(Grid)), zeros(size(Grid))
    add_layer!(P1, Tm, Grid; zlim = (-40.0, f), phase = ConstantPhase(1))
    add_layer!(P2, Tm, Grid; zlim = (-40.0, Surf), phase = ConstantPhase(1))
    @test P1 == P2 && count(==(1), P1) > 0
    Small = CartData(xyz_grid(0:1.0:50, 0:1.0:10, 0))
    @test_throws "the surface must cover xlim and ylim" add_layer!(P2, Tm, Grid; zlim = (-40.0, Small))
    P3 = copy(P2)
    add_layer!(P3, Tm, Grid; xlim = (0, 50), zlim = (-40.0, Small), phase = ConstantPhase(2))   # within the surface
    @test count(==(2), P3) > 0
end

# GaussianTemp adds to the current temperature
let
    Grid = CartData(xyz_grid(-50:1.0:50, -50:1.0:50, -100:1.0:0))
    Z = Grid.z.val
    Ph = zeros(Int, size(Grid)); T0 = 1000 .- 10 .* Z; Tm = copy(T0)
    add_sphere!(Ph, Tm, Grid; cen = (0.0, 0.0, -50.0), radius = 40.0, T = GaussianTemp(Amplitude = 200, Width = 10))
    @test Tm[51, 51, 51] ≈ T0[51, 51, 51] + 200
    @test Tm[61, 51, 51] ≈ T0[61, 51, 51] + 200exp(-0.5)
    @test Tm[1, 1, 1] == T0[1, 1, 1]
end

# GaussianTemp is centered in the region and may be constant along an axis
let
    Grid = CartData(xyz_grid(-50:1.0:50, -50:1.0:50, -100:1.0:0))
    Ph = zeros(Int, size(Grid)); Tm = fill(1350.0, size(Grid))
    add_sphere!(Ph, Tm, Grid; cen = (0.0, 0.0, -50.0), radius = 40.0, T = GaussianTemp(Amplitude = 200, Width = 10))
    @test Tm[51, 51, 51] ≈ 1550
    @test Tm[61, 51, 51] ≈ 1350 + 200exp(-0.5)
    @test Tm[51, 51, 81] ≈ 1350 + 200exp(-4.5)
    Tm = fill(1350.0, size(Grid))
    add_box!(Ph, Tm, Grid; xlim = (-50.0, 50.0), ylim = (-50.0, 50.0), zlim = (-100.0, 0.0), T = GaussianTemp(Amplitude = -100, Width = (10, 10, Inf)))
    @test all(≈(1250), Tm[51, 51, :])
    @test Tm[61, 51, 1] ≈ Tm[61, 51, end] ≈ 1350 - 100exp(-0.5)
    @test_throws "Width must be a number or have 3 entries" add_sphere!(Ph, Tm, Grid; radius = 10, T = GaussianTemp(Width = (1, 2)))
end

# Thermal structures are anchored at the top of ellipsoids, cylinders and faults (#135); extents lie between grid nodes
let
    Grid = CartData(xyz_grid(-10:1.0:10, -10:1.0:10, -100:1.0:0))
    Ph = zeros(Int, size(Grid)); Tm = zeros(size(Grid))
    add_cylinder!(Ph, Tm, Grid; base = (0.0, 0.0, -80.5), cap = (0.0, 0.0, -20.5), radius = 5.0, T = LinearTemp(Ttop = 20, Tbot = 1000))
    @test Tm[11, 11, 80] ≈ 20 + 0.5 / 60 * 980
    @test Tm[11, 11, 21] ≈ 20 + 59.5 / 60 * 980

    Ph = zeros(Int, size(Grid)); Tm = zeros(size(Grid))
    add_ellipsoid!(Ph, Tm, Grid; cen = (0.0, 0.0, -50.5), axes = (5.0, 5.0, 30.0), T = LinearTemp(Ttop = 20, Tbot = 1000))
    @test Tm[11, 11, 80] ≈ 20 + 0.5 / 60 * 980
    @test Tm[11, 11, 21] ≈ 20 + 59.5 / 60 * 980
end

# 2D Phase and Temp arrays of size (nx, nz) give the same result as (nx, 1, nz) arrays
Grid2D = CartData(xyz_grid(range(-100.0, 100.0, 64), 0, range(-60.0, 5.0, 32)))
trench2D = Trench(Start = (0.0, -100.0), End = (0.0, 100.0), Thickness = 80.0, θ_max = 45.0, Length = 300, Lb = 200, direction = -1.0)
setups2D = (
    (P, T) -> add_slab!(P, T, Grid2D, trench2D; phase = ConstantPhase(1), T = HalfspaceCoolingTemp()),
    (P, T) -> add_polygon!(P, T, Grid2D; xlim = (20.0, 80.0, 80.0), zlim = (-5.0, -5.0, -40.0), phase = ConstantPhase(8), T = LinearTemp()),
    (P, T) -> add_stripes!(P, Grid2D; stripAxes = (1, 0, 0), stripeWidth = 2, stripeSpacing = 10, phase = ConstantPhase(8), stripePhase = ConstantPhase(5)),
)
P3, T3, P2, T2 = zeros(Int64, 64, 1, 32), zeros(64, 1, 32), zeros(Int64, 64, 32), zeros(64, 32)
for (setup, phase) in zip(setups2D, (1, 8, 5))
    setup(P3, T3)
    setup(P2, T2)
    @test vec(P2) == vec(P3) && vec(T2) == vec(T3)
    @test phase in P2
end
@test_throws DimensionMismatch add_stripes!(zeros(Int64, 64, 31), Grid2D; stripAxes = (1, 0, 0), phase = ConstantPhase(0))

# SpreadingRateTemp: MORside = "back" measures the plate age from the back (maximum y) edge, mirroring "front"
GridMOR = CartData(xyz_grid(0:10:100, 0:10:100, -100:10:0))
T_front, T_back, P_MOR = zeros(size(GridMOR)), zeros(size(GridMOR)), zeros(Int64, size(GridMOR))
add_box!(P_MOR, T_front, GridMOR; xlim = (0, 100), ylim = (0, 100), zlim = (-100, 0), T = SpreadingRateTemp(MORside = "front"))
add_box!(P_MOR, T_back, GridMOR; xlim = (0, 100), ylim = (0, 100), zlim = (-100, 0), T = SpreadingRateTemp(MORside = "back"))
@test T_back == reverse(T_front, dims = 2)
@test T_front != reverse(T_front, dims = 2)

# The setup routines evaluate points one by one and allocate nothing per point
GridAlloc = CartData(xyz_grid(range(0, 100, 40), range(0, 100, 40), range(-60, 0, 40)))
PAlloc, TAlloc = zeros(Int64, size(GridAlloc)), zeros(size(GridAlloc))
lithAlloc = LithosphericPhases(Layers = [10, 20], Phases = [1, 2, 3], Tlab = 1250)
trenchAlloc = Trench(Start = (40.0, 0.0), End = (40.0, 100.0), Thickness = 20.0, θ_max = 45.0, Length = 60.0, Lb = 40.0, n_seg = 20)
surfAlloc = CartData(xyz_grid(range(-10, 110, 25), range(-10, 110, 25), 0))
surfAlloc = CartData(surfAlloc.x.val, surfAlloc.y.val, -5 .+ 0.05 .* surfAlloc.x.val, (; Z = surfAlloc.z.val))
setupsAlloc = (
    () -> add_box!(PAlloc, TAlloc, GridAlloc; xlim = (10, 90), ylim = (10, 90), zlim = (-40, 0), DipAngle = 20, StrikeAngle = 30, phase = lithAlloc, T = HalfspaceCoolingTemp(Age = 80)),
    () -> add_sphere!(PAlloc, TAlloc, GridAlloc; cen = (50, 50, -30), radius = 15, T = LinearTemp()),
    () -> add_layer!(PAlloc, TAlloc, GridAlloc; zlim = (-20, 0), phase = lithAlloc, T = HalfspaceCoolingTemp()),
    () -> add_slab!(PAlloc, TAlloc, GridAlloc, trenchAlloc; phase = ConstantPhase(4), T = HalfspaceCoolingTemp()),
    () -> add_layer!(PAlloc, TAlloc, GridAlloc; zlim = ((x, y) -> -30 + 5sind(x), (x, y) -> -0.1x), phase = lithAlloc, T = HalfspaceCoolingTemp()),
    () -> add_fault!(PAlloc, TAlloc, GridAlloc; Start = (20, 0), End = (20, 100), DipAngle = 30, DetachmentDepth = 20),
    () -> add_fault!(PAlloc, TAlloc, GridAlloc; Trace = [(50 + 30cosd(α), 50 + 30sind(α)) for α in 0:10:180], DipAngle = -30),
    () -> add_slab!(PAlloc, TAlloc, GridAlloc, Trench(Trace = [(50 + 40cosd(α), 40sind(α)) for α in 30:10:150], Thickness = 20.0, Length = 60.0, Lb = 40.0, n_seg = 20); T = HalfspaceCoolingTemp()),
    () -> add_layer!(PAlloc, TAlloc, GridAlloc; zlim = (-40, surfAlloc), T = LinearTemp()),
    () -> add_sphere!(PAlloc, TAlloc, GridAlloc; cen = (50, 50, -30), radius = 15, T = GaussianTemp()),
    () -> add_sphere!(PAlloc, TAlloc, GridAlloc; cen = (50, 50, -30), radius = 15, T = GaussianTemp(Width = (5, 5, Inf))),
    () -> add_box!(PAlloc, TAlloc, GridAlloc; xlim = (0, 100), zlim = (-40, 0), T = SpreadingRateTemp(), segments = [((20.0, 0.0), (20.0, 40.0)), ((60.0, 40.0), (60.0, 70.0)), ((30.0, 70.0), (30.0, 100.0))]),
)
for setup in setupsAlloc
    setup()
    @test @allocated(setup()) < 50_000   # 64_000 points
end
# a BezierSlab allocates samples of its surface, but nothing per point
let
    slab = BezierSlab(
        Points = [(x + h, y, zz) for (x, y) in ((20.0, 0.0), (50.0, 50.0), (20.0, 100.0)), (h, zz) in ((-20.0, 0.0), (10.0, 0.0), (30.0, -30.0))],
        Thickness = 15.0, Tear = Tear(; Depth = 10, Width = 5, Start = 20, Tip = 60)
    )
    function allocs(G)
        P, T = zeros(Int64, size(G)), zeros(size(G))
        add_slab!(P, T, G, slab; T = HalfspaceCoolingTemp())
        return @allocated add_slab!(P, T, G, slab; T = HalfspaceCoolingTemp())
    end
    @test allocs(CartData(xyz_grid(range(0, 100, 80), range(0, 100, 80), range(-60, 0, 80)))) - allocs(GridAlloc) < 50_000   # 512_000 and 64_000 points
end

# add_layer! measures phases and temperatures from the top of the layer, for a numeric top as for a function
let
    Grid = CartData(xyz_grid(0:10.0:100, 0, -120:1.0:0))
    lith = LithosphericPhases(Layers = [20, 50], Phases = [1, 2, 3])
    P1, P2, T1, T2 = zeros(Int, size(Grid)), zeros(Int, size(Grid)), zeros(size(Grid)), zeros(size(Grid))
    add_layer!(P1, T1, Grid; zlim = (-100, -10), phase = lith, T = HalfspaceCoolingTemp(Adiabat = 0.5))
    add_layer!(P2, T2, Grid; zlim = (-100, (x, y) -> -10), phase = lith, T = HalfspaceCoolingTemp(Adiabat = 0.5))
    @test P1 == P2 && T1 ≈ T2
    @test count(==(1), P1[1, 1, :]) == 20   # crust from -10 to -30

    # numeric bounds in either order
    P3 = zeros(Int, size(Grid))
    add_layer!(P3, T2, Grid; xlim = (100, 0), zlim = (-10, -100), phase = lith, T = HalfspaceCoolingTemp(Adiabat = 0.5))
    @test P3 == P1
    @test_throws "lies below its bottom" add_layer!(P3, T2, Grid; zlim = ((x, y) -> -10, -100))

    P4, P5 = zeros(Int, size(Grid)), zeros(Int, size(Grid))
    add_plate!(P4, T2, Grid; xlim = (0, 50, 50, 0), ylim = (-1, -1, 1, 1), zlim = (-50, 0))
    add_plate!(P5, T2, Grid; xlim = (0, 50, 50, 0), ylim = (-1, -1, 1, 1), zlim = (0, -50))
    @test P4 == P5 && count(==(1), P4) > 0
end

# LinearWeightedTemperature in add_slab! blends over the length along the slab down to the decoupling depth
let
    x, z = -100:2.0:400, -300:2.0:0
    Grid = CartData(xyz_grid(x, 0, z))
    # a planar slab dipping 60° from the trench at x = 0
    t = Trench(Start = (0.0, -100.0), End = (0.0, 100.0), Thickness = 40.0, Bending = [(60, 1.0), (60, 300.0)], n_seg = 301, d_decoupling = 100.0)
    Ph = zeros(Int, size(Grid)); Tm = zeros(size(Grid))
    add_slab!(Ph, Tm, Grid, t; T = LinearWeightedTemperature(F1 = ConstantTemp(0), F2 = ConstantTemp(1000)))
    ls = Grid.x.val .* cosd(60) .- Grid.z.val .* sind(60)
    blend = (Ph .== 1) .& (Tm .< 999)
    @test maximum(abs.(Tm[blend] .- 1000 .* ls[blend] ./ (100 / sind(60)))) < 10
end

# setup routines accept Float32 grids
let
    G64 = CartData(xyz_grid(0:5.0:200, 0:5.0:200, -100:5.0:0))
    G32 = CartData(xyz_grid(0:5.0f0:200, 0:5.0f0:200, -100:5.0f0:0))
    arc = [(100 + 80cosd(α), 80sind(α)) for α in 20:10:160]
    setups = (
        (P, T, G) -> add_box!(P, T, G; xlim = (20, 150), zlim = (-50, 0), DipAngle = 10),
        (P, T, G) -> add_plate!(P, T, G; xlim = (20, 150, 150, 20), ylim = (20, 20, 150, 150), zlim = (-50, 0)),
        (P, T, G) -> add_slab!(P, T, G, Trench(Start = (100.0, 0.0), End = (100.0, 200.0), Thickness = 30.0, Length = 100.0, Lb = 60.0)),
        (P, T, G) -> add_slab!(P, T, G, Trench(Trace = arc, Thickness = 30.0, Length = 100.0, Lb = 60.0)),
        (P, T, G) -> add_slab!(P, T, G, Trench(Trace = arc, Thickness = 30.0, Length = 100.0, Lb = 60.0, Tear = Tear(; Depth = 40, Width = 20, Start = 50, Tip = (150, 100)))),
    )
    for setup in setups
        P64, P32 = zeros(Int, size(G64)), zeros(Int, size(G32))
        setup(P64, zeros(size(G64)), G64)
        setup(P32, zeros(Float32, size(G32)), G32)
        @test P32 == P64 && count(==(1), P64) > 0
    end
end

# a curved trench skips the points beyond the reach of the slab without changing the result
let
    Grid = CartData(xyz_grid(-200:10.0:800, -200:10.0:800, -300:10.0:0))
    arc = [(300 + 300cosd(α), 100 + 300sind(α)) for α in 20:2:160]
    t = Trench(; Trace = arc, Thickness = 50.0, θ_max = 45.0, Length = 200.0, Lb = 100.0, direction = 1.0, WeakzoneThickness = 10.0, WeakzonePhase = 3, ThermalThickness = 80.0)
    Ph = zeros(Int, size(Grid)); Tm = fill(1350.0, size(Grid))
    add_slab!(Ph, Tm, Grid, t; phase = ConstantPhase(1), T = LinearTemp(Ttop = 0, Tbot = 1300))
    # the same counts and temperatures as when the offset from the trace is evaluated at every point
    @test count(==(1), Ph) == 5221 && count(==(3), Ph) == 1026
    @test sum(Tm) ≈ 4.2121549043863046e8
end

# invalid inputs of slabs and thermal structures
let
    Grid = CartData(xyz_grid(0:10.0:100, 0:10.0:100, -100:10.0:0))
    Ph, Tm = zeros(Int, size(Grid)), zeros(size(Grid))
    vertical = [(50.0, y, z) for y in (0.0, 100.0), z in (0.0, -50.0, -100.0)]
    @test_throws "is vertical" add_slab!(Ph, Tm, Grid, BezierSlab(Points = vertical))
    @test_throws "at least one (angle, length) pair" add_slab!(Ph, Tm, Grid, Trench(Bending = Tuple{Float64, Float64}[]))
    segs = [((0.0, 0.0), (0.0, 100.0))]
    @test_throws "segments only apply to SpreadingRateTemp" add_box!(Ph, Tm, Grid; xlim = (0, 100), zlim = (-50, 0), T = GaussianTemp(), segments = segs)
    @test_throws "segments only apply to SpreadingRateTemp, got Nothing" add_box!(Ph, Tm, Grid; xlim = (0, 100), zlim = (-50, 0), segments = segs)
    @test_throws "segments only apply to SpreadingRateTemp, got Nothing" add_plate!(Ph, Tm, Grid; xlim = (0, 100, 100), ylim = (0, 0, 100), zlim = (-50, 0), segments = segs)
end

# a layer outside the grid sets nothing, whether its bounds are numbers or functions
let
    Grid = CartData(xyz_grid(0:10.0:100, 0:10.0:100, -100:10.0:0))
    Ph, Tm = zeros(Int, size(Grid)), zeros(size(Grid))
    add_layer!(Ph, Tm, Grid; xlim = (200, 300), zlim = (-50, -10))
    add_layer!(Ph, Tm, Grid; xlim = (200, 300), zlim = (-50, (x, y) -> -10))
    @test all(iszero, Ph) && all(iszero, Tm)
end

# the decoupling depth of a slab sets the crit_dist of a LinearWeightedTemperature along the slab (dir = :X) only
let
    Grid = CartData(xyz_grid(0:10.0:400, 0:10.0:100, -200:10.0:0))
    t = Trench(Start = (100.0, 0.0), End = (100.0, 100.0), Thickness = 80.0, θ_max = 45.0, Length = 300.0, Lb = 100.0, d_decoupling = 50.0)
    function slab_temperature(crit_dist)
        Ph, Tm = zeros(Int, size(Grid)), zeros(size(Grid))
        add_slab!(Ph, Tm, Grid, t; T = LinearWeightedTemperature(; crit_dist, dir = :Z, F1 = ConstantTemp(0), F2 = ConstantTemp(1000)))
        return Tm
    end
    @test slab_temperature(-10.0) != slab_temperature(-1000.0)    # the depth below the slab top, d ≤ 0, over crit_dist
end

# Reversed limits, inclusive sphere boundary, segment order, argument checks
let
    G = CartData(xyz_grid(0:1.0:10, 0:1.0:10, -10:1.0:0))
    P() = zeros(Int, size(G)), zeros(size(G))
    lith = LithosphericPhases(Layers = [2 3], Phases = [1 2 3])
    (P1, T1), (P2, T2) = P(), P()
    add_box!(P1, T1, G; zlim = (-10, 0), phase = lith); add_box!(P2, T2, G; zlim = (0, -10), phase = lith)
    @test P1 == P2
    add_box!(P1, T1, G; xlim = (10, 0), zlim = (-10, 0), T = McKenzie_subducting_slab())
    @test all(abs.(T1) .< 1.0e4)
    (P1, T1), (P2, T2) = P(), P()
    add_fault!(P1, T1, G; Start = (5.0, 0.0), End = (5.0, 10.0), Depth_extent = (0, -5)); add_fault!(P2, T2, G; Start = (5.0, 0.0), End = (5.0, 10.0), Depth_extent = (-5, 0))
    @test sum(P1) > 0 && P1 == P2
    (P1, T1), (P2, T2) = P(), P()
    poly = (xlim = [2.0, 8.0, 8.0, 2.0], zlim = [-8.0, -8.0, -2.0, -2.0])
    add_polygon!(P1, T1, G; poly..., ylim = (8, 2)); add_polygon!(P2, T2, G; poly..., ylim = (2, 8))
    @test sum(P1) > 0 && P1 == P2
    Ps, Ts = P()
    add_sphere!(Ps, Ts, G; cen = (5.0, 5.0, -5.0), radius = 1.0)
    @test sum(Ps) == 7
    seg = [((0, 0), (100, 0)), ((100, 50), (200, 50))]
    @test GeophysicalModelGenerator.determine_region(50, 10, seg) == 1
    @test GeophysicalModelGenerator.determine_region(150, 40, seg) == 2
    @test_throws "segments must not be empty" add_box!(Ps, Ts, G; T = SpreadingRateTemp(), segments = typeof(seg[1])[])
    @test_throws "base and cap must differ" add_cylinder!(Ps, Ts, G; base = (1.0, 1.0, 1.0), cap = (1.0, 1.0, 1.0), radius = 1.0)
    @test_throws "stripeWidth must be" add_stripes!(Ps, G; stripAxes = (1, 0, 0), stripeWidth = 1.0, stripeSpacing = 1.0)
    @test_throws "stripAxes needs" add_stripes!(Ps, G; stripAxes = (0, 0, 0), stripeWidth = 0.1, stripeSpacing = 1.0)
    t = Trench(Start = (100.0, 0.0), End = (100.0, 100.0), Thickness = 80.0, θ_max = 45.0, Length = 300.0, Lb = 100.0, d_decoupling = 50.0, type_bending = :Customize)
    @test_throws "type_bending must be" add_slab!(Ps, Ts, G, t)
end

# LithosphericTemp gives a finite geotherm for any order of the layer phases and rejects a region of zero height;
# add_slab! rejects a d_decoupling that the slab does not reach or that is not positive
let
    G = CartData(xyz_grid(0:10.0:100, 0:10.0:100, -100:10.0:0)); P = zeros(Int, size(G)); Tm = zeros(size(G))
    add_box!(P, Tm, G; xlim = (0, 100), zlim = (-100, 0), phase = LithosphericPhases(Layers = [20 15 65], Phases = [3 2 1]), T = LithosphericTemp(nz = 51))
    @test all(isfinite, Tm)
    @test_throws "nonzero height" GeophysicalModelGenerator.prepare_temperature(LithosphericTemp(), nothing, (zeros(2), zeros(2), zeros(2), ones(Int, 2)))
    Tw = LinearWeightedTemperature(0, 1, 600.0, :X, ConstantTemp(1000), ConstantTemp(2000))
    tr(d) = Trench(Start = (50.0, 0.0), End = (50.0, 100.0), Thickness = 20.0, θ_max = 45.0, Length = 300.0, Lb = 100.0, d_decoupling = d)
    @test_throws "d_decoupling must be greater" add_slab!(P, Tm, G, tr(0.0); T = Tw)
    @test_throws "does not reach d_decoupling" add_slab!(P, Tm, G, tr(1.0e5); T = Tw)
    t = Trench(Start = (50.0, 0.0), End = (50.0, 100.0), Thickness = 20.0, Length = 100.0, Lb = 50.0, ThermalThickness = 40.0)
    @test_throws "ThermalThickness cannot be combined with LithosphericTemp" add_slab!(P, Tm, G, t; T = LithosphericTemp())
end

# Ridge segments: fracture zones continue the transform faults between consecutive segments, and each point is
# aged from the segment between the fracture zones that enclose it, even when another segment is nearer
let
    region(x, y, segs) = GeophysicalModelGenerator.determine_region(x, y, segs)
    offset = [((-500.0, -1000.0), (-500.0, 0.0)), ((-250.0, 0.0), (-250.0, 200.0)), ((-750.0, 200.0), (-750.0, 1000.0))]
    @test region(-400.0, 100.0, offset) == 2
    @test region(-100.0, -50.0, offset) == 1
    @test region(-600.0, 250.0, offset) == 3
    @test region(5000.0, 1500.0, offset) == 3     # beyond the segment ends
    @test region(-500.0, -5000.0, offset) == 1
    oblique = [((-600.0, -1000.0), (-500.0, 0.0)), ((-300.0, 0.0), (-280.0, 200.0)), ((-750.0, 200.0), (-670.0, 1000.0))]
    @test region(-450.0, 100.0, oblique) == 2
    @test region(-100.0, -50.0, oblique) == 1
    @test region(5000.0, 1500.0, oblique) == 3
    @test region(-1000.0, -3000.0, oblique) == 1
    alongx = [((0.0, 0.0), (100.0, 0.0)), ((100.0, 50.0), (200.0, 50.0))]   # ridge along x, spreading along y
    @test region(90.0, 40.0, alongx) == 1
    @test region(150.0, 40.0, alongx) == 2
    @test_throws "must lie on opposite sides of the transform fault" add_box!(zeros(Int, 2, 2, 2), zeros(2, 2, 2), CartData(xyz_grid(0:1.0, 0:1.0, -1:1.0:0)); xlim = (0, 1), zlim = (-1, 0), T = SpreadingRateTemp(), segments = [((0.0, 0.0), (0.0, 100.0)), ((0.0, 100.0), (50.0, 200.0))])

    # across the fracture zone y = 0 the plate age jumps from 100 km off segment 1 to 150 km off segment 2,
    # as at the same distances from a single straight ridge
    s = SpreadingRateTemp(SpreadingVel = 3)
    Tseg, Tref = zeros(2), zeros(2)
    compute_thermal_structure(Tseg, [-400.0, -400.0], [-1.0, 1.0], [-20.0, -20.0], zeros(Int, 2), s, offset; ztop = 0.0)
    compute_thermal_structure(Tref, [100.0, 150.0], [0.0, 0.0], [-20.0, -20.0], zeros(Int, 2), s, [((0.0, -1000.0), (0.0, 1000.0))]; ztop = 0.0)
    @test Tseg ≈ Tref
    @test Tseg[1] > Tseg[2]
end

# A surface bound of a layer must have no NaN heights; a surface with a single row in y is constant in y
let
    G = CartData(xyz_grid(0:10.0:100, 0:10.0:20, -50:5.0:0))
    P, Tm = zeros(Int, size(G)), zeros(size(G))
    S = CartData(xyz_grid(0:10.0:100, 0:10.0:20, 0))
    S = CartData(S.x.val, S.y.val, -20 .- 0.1 .* S.x.val, (; Z = S.z.val))
    Sn = CartData(S.x.val, S.y.val, copy(S.z.val), (; Z = S.z.val))
    Sn.z.val[3, 2, 1] = NaN
    @test_throws "fill or remove them" add_layer!(P, Tm, G; zlim = (Sn, 0.0), phase = ConstantPhase(1))
    S1 = CartData(xyz_grid(0:10.0:100, 0, 0))
    S1 = CartData(S1.x.val, S1.y.val, -20 .- 0.1 .* S1.x.val, (; Z = S1.z.val))
    P1 = zeros(Int, size(G))
    add_layer!(P, Tm, G; zlim = (S, 0.0), phase = ConstantPhase(1))
    add_layer!(P1, Tm, G; zlim = (S1, 0.0), phase = ConstantPhase(1))
    @test P1 == P && count(==(1), P) > 0
end

# every width of GaussianTemp must be positive
let
    G = CartData(xyz_grid(0:10.0:100, 0:10.0:100, -100:10.0:0)); P = zeros(Int, size(G)); Tm = zeros(size(G))
    @test_throws "every Width must be positive" add_sphere!(P, Tm, G; cen = (50, 50, -50), radius = 30, T = GaussianTemp(Width = 0))
    @test_throws "every Width must be positive" add_sphere!(P, Tm, G; cen = (50, 50, -50), radius = 30, T = GaussianTemp(Width = (10, -1, Inf)))
end

@testset "add_plate! measures depth below the top, like add_box!" begin
    G = CartData(xyz_grid(0:1.0:2, 0:1.0:2, -100:10.0:-20))
    kw = (; zlim = (-100, -20), phase = LithosphericPhases(Layers = [15, 30], Phases = [1, 2, 3]), T = HalfspaceCoolingTemp(Adiabat = 0.5))
    Pb, Tb = zeros(Int32, size(G)), zeros(size(G))
    Pp, Tp = zeros(Int32, size(G)), zeros(size(G))
    add_box!(Pb, Tb, G; xlim = (0, 2), ylim = (0, 2), kw...)
    add_plate!(Pp, Tp, G; xlim = (0, 2, 2, 0), ylim = (0, 0, 2, 2), kw...)
    @test Pp == Pb
    @test Tp ≈ Tb
    @test Pp[1, 1, end] == 1
end

@testset "add_box! measures depth from its top for any Origin" begin
    Grid = CartData(xyz_grid(0:5:100, 0:5:50, -200:5:0))
    lith = LithosphericPhases(Layers = [20, 60], Phases = [1, 2, 3])
    function box(Origin)
        P, T = zeros(Int64, size(Grid)), zeros(size(Grid))
        add_box!(P, T, Grid; xlim = (20, 80), ylim = (0, 50), zlim = (-150, -50), Origin, phase = lith, T = HalfspaceCoolingTemp(Age = 30))
        return P, T
    end
    P_top, T_top = box(nothing)
    P_mid, T_mid = box((50, 25, -100))
    @test P_mid == P_top && T_mid == T_top
    @test sort(unique(P_top[P_top .> 0])) == [1, 2, 3]              # layers counted from the top at -50 km
end

@testset "LinearTemp spans the local thickness of a curved layer" begin
    G = CartData(xyz_grid(0:2.0:100, 0, -100:1.0:0))
    X, Z = G.x.val, G.z.val
    top(x, y) = -20 + 5sinpi(x / 25)
    P, T = zeros(Int, size(G)), zeros(size(G))
    add_layer!(P, T, G; zlim = (-70, top), phase = ConstantPhase(1), T = LinearTemp(Ttop = 100, Tbot = 900))
    for i in axes(P, 1)
        k = findall(==(1), P[i, 1, :])
        H = top(X[i, 1, 1], 0) + 70
        @test T[i, 1, last(k)] - 100 <= 800 / H + 1.0e-9     # within one cell below the local top
        @test Z[i, 1, first(k)] == -70 && T[i, 1, first(k)] ≈ 900
    end
end

@testset "listric fault does not extend beyond its up-dip side" begin
    G = CartData(xyz_grid(-100:2.0:100, 0:10.0:20, -30:1.0:0))
    P, T = zeros(Int, size(G)), zeros(size(G))
    # dips towards +x; far on the -x side, the listric surface would approach z = -D everywhere
    add_fault!(P, T, G; Start = (0, 0), End = (0, 20), Fault_thickness = 2, DipAngle = 5, DetachmentDepth = 10, phase = ConstantPhase(1))
    @test all(==(0), P[G.x.val .<= -10])
    @test P[51, 2, end] == 1                         # the trace at the surface
    @test any(==(1), P[(G.x.val .> 10) .& (G.z.val .< -9)])   # the detachment
end

@testset "invalid geometries raise errors" begin
    G = CartData(xyz_grid(0:5.0:400, 0:5.0:400, -200:5.0:0))
    P, T = zeros(Int, size(G)), zeros(size(G))
    # the edge v = 0 collapses to a point
    pts = [(j == 0 ? 100.0 : 100.0 + 100j, j == 0 ? 200.0 : 100.0i, -50.0j) for i in 1:4, j in 0:3]
    @test_throws "the Bézier surface is degenerate" add_slab!(P, T, G, BezierSlab(Points = pts, Thickness = 30.0))
    @test_throws "must not be closed" add_fault!(P, T, G; Trace = [(0, 0), (100, 0), (100, 100), (0, 0)])
    @test_throws "must not turn by more than 90°" add_fault!(P, T, G; Trace = [(0, 0), (100, 0), (0, 10)])
    @test_throws "must not turn by more than 90°" add_slab!(P, T, G, Trench(Trace = [(0.0, 0.0), (100.0, 0.0), (0.0, 10.0)]))
    @test_throws "angles of Bending must lie in (0, 90]" compute_slab_surface(Trench(Bending = [(0, 100)]))
    @test_throws "angles of Bending must lie in (0, 90]" compute_slab_surface(Trench(Bending = [(30, 100), (95, 100)]))
    @test_throws "lengths of Bending must be positive" compute_slab_surface(Trench(Bending = [(30, -100)]))
end

@testset "depth below the top of a sphere, polygon and plate" begin
    G = CartData(xyz_grid(-50:2.0:50, -50:2.0:50, -100:1.0:0))
    X, Y, Z = G.x.val, G.y.val, G.z.val
    HS = HalfspaceCoolingTemp(Tsurface = 20, Age = 30)
    shapes = [
        (P, T) -> add_sphere!(P, T, G; cen = (0, 0, -50), radius = 30, T = HS),
        (P, T) -> add_polygon!(P, T, G; xlim = (-30, 30, 30, -30), ylim = (-20, 20), zlim = (-80, -80, -20, -20), T = HS),
        (P, T) -> add_plate!(P, T, G; xlim = (-30, 30, 30, -30), ylim = (-30, -30, 30, 30), zlim = (-80, -20), T = HS),
    ]
    for f in shapes
        P, T = zeros(Int, size(G)), fill(1000.0, size(G))
        f(P, T)
        i = findfirst(i -> X[i] == 0 && Y[i] == 0 && Z[i] == -20, eachindex(X))
        @test P[i] == 1 && T[i] ≈ 20
        @test all(T[(P .== 1) .& (Z .< -20)] .> 20)
    end
end

@testset "GaussianTemp in a cylinder" begin
    G = CartData(xyz_grid(-50:2.0:50, -50:2.0:50, -100:1.0:0))
    P, T = zeros(Int, size(G)), fill(1000.0, size(G))
    add_cylinder!(P, T, G; base = (0, 0, -90), cap = (0, 0, -10), radius = 30, T = GaussianTemp(Amplitude = 100, Width = (10, 10, Inf)))
    axis = (G.x.val .== 0) .& (G.y.val .== 0) .& (P .== 1)
    @test count(axis) == 81 && all(T[axis] .≈ 1100)
    @test maximum(T) ≈ 1100
end

@testset "phase = nothing keeps the phases of a sphere" begin
    G = CartData(xyz_grid(-50:2.0:50, -50:2.0:50, -100:1.0:0))
    P, T = fill(7, size(G)), zeros(size(G))
    add_sphere!(P, T, G; cen = (0, 0, -50), radius = 30, phase = nothing, T = ConstantTemp(500))
    @test all(==(7), P)
    @test count(==(500), T) > 0 && all(t -> t == 0 || t == 500, T)
end

# a crater makes a truncated cone: flat top of radius `crater`, linear flanks down to `base` at `radius` (see #146)
@testset "add_volcano! with a crater" begin
    x, z = range(0.0, 50, 251), range(-5.0, 5.0, 101)
    Grid = CartData(xyz_grid(x, 0:50:50, z))
    for crater in (0.0, 2.0)
        Phases, Temp = zeros(Int64, size(Grid)), zeros(size(Grid))
        add_volcano!(Phases, Temp, Grid; volcanic_phase = 1, center = (25.0, 0.0), height = 3, radius = 5, crater, base = 0.0, background = nothing, T = HalfspaceCoolingTemp(Age = 20))
        top = [maximum(z[k] for k in eachindex(z) if Phases[i, 1, k] == 1; init = 0.0) for i in eachindex(x)]
        cone(r) = r <= crater ? 3.0 : r <= 5 ? 3.0 * (5 - r) / (5 - crater) : 0.0
        @test maximum(abs.(top .- cone.(abs.(x .- 25)))) <= step(z) + 1.0e-12
    end
end
