# test setting geometries in the different grid types
using Test, GeophysicalModelGenerator, GeoParams

# GeoData
Lon3D,Lat3D,Depth3D =   LonLatDepthGrid(1.0:1:10.0, 11.0:1:20.0, (-20:1:-10)*km);
Data                =   zeros(size(Lon3D));
Temp                =   ones(Float64, size(Data))*1350;
Phases              =   zeros(Int32,   size(Data));
Grid                =   GeoData(Lon3D,Lat3D,Depth3D,(DataFieldName=Data,))   

AddBox!(Phases,Temp,Grid, xlim=(2,4), zlim=(-15,-10), phase=ConstantPhase(3), DipAngle=10, T=LinearTemp(Tbot=1350, Ttop=200))
@test sum(Temp[1,1,:]) ≈ 14850.0

AddEllipsoid!(Phases,Temp,Grid, cen=(4,15,-17), axes=(1,2,3), StrikeAngle=90, DipAngle=45, phase=ConstantPhase(2), T=ConstantTemp(600))
@test sum(Temp[1,1,:]) ≈ 14850.0

# CartData 
X,Y,Z               =   XYZGrid(1.0:1:10.0, 11.0:1:20.0, -20:1:-10);
Data                =   zeros(size(X));
Temp                =   ones(Float64, size(Data))*1350;
Phases              =   zeros(Int32,   size(Data));
Grid                =   CartData(X,Y,Z,(DataFieldName=Data,))   

AddBox!(Phases,Temp,Grid, xlim=(2,4), zlim=(-15,-10), phase=ConstantPhase(3), DipAngle=10, T=LinearTemp(Tbot=1350, Ttop=200))
@test sum(Temp[1,1,:]) ≈ 14850.0

AddEllipsoid!(Phases,Temp,Grid, cen=(4,15,-17), axes=(1,2,3), StrikeAngle=90, DipAngle=45, phase=ConstantPhase(2), T=ConstantTemp(600))
@test sum(Temp[1,1,:]) ≈ 14850.0

# CartGrid
Grid                =   CreateCartGrid(size=(10,20,30),x=(0.,10), y=(0.,10), z=(2.,10))
Temp                =   ones(Float64, Grid.N...)*1350;
Phases              =   zeros(Int32,  Grid.N...);

AddBox!(Phases,Temp,Grid, xlim=(2,4), zlim=(4,8), phase=ConstantPhase(3), DipAngle=10, T=LinearTemp(Tbot=1350, Ttop=200))
@test maximum(Phases) == 3

addStripes!(Phases, Grid,stripAxes=(1,1,1),stripeWidth=0.2,stripeSpacing=1,Origin=nothing, StrikeAngle=0, DipAngle=10,phase = ConstantPhase(3),stripePhase = ConstantPhase(4))
@test maximum(Phases) == 4

# Create a CartData structure from it
Data = CartData(Grid, (T=Temp, Phases=Phases))

@test NumValue(Data.x[3,3,2]) ≈ 2.2222222222222223

# Doing the same for vertical cross-sections
Grid2D              =   CreateCartGrid(size=(10,30),x=(0.,10), z=(2.,10))
Temp2D              =   ones(Float64, Grid2D.N...)*1350;
Phases2D            =   zeros(Int32,  Grid2D.N...);

Data2D = CartData(Grid2D, (T=Temp2D, Phases=Phases2D))

@test NumValue(Data.x[3,1,2]) ≈ 2.2222222222222223


# LithosphericPhases
LP                  =   LithosphericPhases(Layers=[5 10 6], Phases=[0 1 2 3], Tlab=nothing);
X,Y,Z               =   XYZGrid(-5:1:5,-5:1:5,-20:1:5);
Phases              =   zeros(Int32,   size(X));
Temp                =   zeros(Int32,   size(X));
Phases              =   Compute_Phase(Phases, Temp, X, Y, Z, LP);

@test Phases[1,1,end]   == 3
@test Phases[1,1,7]     == 1

Phases              =   Compute_Phase(Phases, Temp, X, Y, Z, LP, Ztop=5);

@test Phases[1,1,end-4] == 0
@test Phases[1,1,5]     == 2

LP                  =   LithosphericPhases(Layers=[0.5 1.0 1.0], Phases=[0 1 2], Tlab=nothing);
Grid                =   ReadLaMEM_InputFile("test_files/SaltModels.dat");
Phases              =   zeros(Int32,   size(Grid.X));
Temp                =   zeros(Int32,   size(Grid.X));
Phases              =   Compute_Phase(Phases, Temp, Grid, LP);

@test Phases[1,1,25]    == 1
@test Phases[1,1,73]    == 0

# Create Grid & nondimensionalize it
CharDim     =   GEO_units();
Grid        =   CreateCartGrid(size=(10,20,30),x=(0.0km,10km), y=(0.0km, 10km), z=(-10.0km, 2.0km), CharDim=CharDim)
@test Grid.Δ[2] ≈ 0.0005263157894736842


# test 1D-explicit thermal solver for AddBox -----------
nel         =   96
Grid        =   CreateCartGrid(size=(nel,nel,nel),x=(-200.,200.), y=(-200.,200.), z=(-200.,0))
Temp        =   zeros(Float64, Grid.N...);
Phases      =   zeros(Int64,  Grid.N...);

# 1) horizontally layer lithosphere; UpperCrust,LowerCrust,Mantle
AddBox!(Phases,Temp,Grid, xlim=(-100,100), zlim=(-100,0), Origin=(0.0,0.0,0.0),
    phase=LithosphericPhases(Layers=[20 15 65], Phases = [1 2 3], Tlab=nothing), 
    DipAngle=0.0, T=LithosphericTemp(nz=201))

@test sum(Temp[Int64(nel/2),Int64(nel/2),:]) ≈ 36131.638045729735

# 2) inclined lithosphere; UpperCrust,LowerCrust,Mantle
Temp    =   zeros(Float64, Grid.N...);
Phases  =   zeros(Int64,  Grid.N...);

AddBox!(Phases,Temp,Grid, xlim=(-100,100), zlim=(-100,0), Origin=(0.0,0.0,0.0),
    phase=LithosphericPhases(Layers=[20 15 65], Phases = [1 2 3], Tlab=nothing), 
    DipAngle=30.0, T=LithosphericTemp(nz=201))

@test sum(Temp[Int64(nel/2),Int64(nel/2),:]) ≈ 41912.18172533137

# 3) inclined lithosphere with respect to the default origin; UpperCrust,LowerCrust,Mantle
Temp    =   zeros(Float64, Grid.N...);
Phases  =   zeros(Int64,  Grid.N...);

AddBox!(Phases,Temp,Grid, xlim=(-100,100), zlim=(-100,0),
    phase=LithosphericPhases(Layers=[20 15 65], Phases = [1 2 3], Tlab=nothing), 
    DipAngle=30.0, T=LithosphericTemp(nz=201))

@test sum(Temp[Int64(nel/2),Int64(nel/2),:]) ≈ 41316.11499878003

# 4) inclined lithosphere with only two layers
Temp    =   zeros(Float64, Grid.N...);
Phases  =   zeros(Int64,  Grid.N...);

ρM=3.0e3            # Density [ kg/m^3 ]
CpM=1.0e3           # Specific heat capacity [ J/kg/K ]
kM=2.3              # Thermal conductivity [ W/m/K ]
HM=0.0              # Radiogenic heat source per mass [H] = W/kg; [H] = [Q/rho]
ρUC=2.7e3           # Density [ kg/m^3 ]
CpUC=1.0e3          # Specific heat capacity [ J/kg/K ]
kUC=3.0             # Thermal conductivity [ W/m/K ]
HUC=617.0e-12       # Radiogenic heat source per mass [H] = W/kg; [H] = [Q/rho]

rheology = (
        # Name              = "UpperCrust",
        SetMaterialParams(;
            Phase               =   1,
            Density             =   ConstantDensity(; ρ=ρUC),
            HeatCapacity        =   ConstantHeatCapacity(; Cp=CpUC),
            Conductivity        =   ConstantConductivity(; k=kUC),
            RadioactiveHeat     =   ConstantRadioactiveHeat(; H_r=HUC*ρUC),     # [H] = W/m^3
        ),
        # Name              = "LithosphericMantle",
        SetMaterialParams(;
            Phase               =   2,
            Density             =   ConstantDensity(; ρ=ρM),
            HeatCapacity        =   ConstantHeatCapacity(; Cp=CpM),
            Conductivity        =   ConstantConductivity(; k=kM),
            RadioactiveHeat     =   ConstantRadioactiveHeat(; H_r=HM*ρM),       # [H] = W/m^3
        ),
    );

AddBox!(Phases,Temp,Grid, xlim=(-100,100), zlim=(-100,0),
    phase=LithosphericPhases(Layers=[20 80], Phases = [1 2], Tlab=nothing), 
    DipAngle=30.0, T=LithosphericTemp(rheology=rheology,nz=201))

@test sum(Temp[Int64(nel/2),Int64(nel/2),:]) ≈ 40513.969831615716

# 5) using flux lower boundary conditions
Temp    =   zeros(Float64, Grid.N...);
Phases  =   zeros(Int64,  Grid.N...);

AddBox!(Phases,Temp,Grid, xlim=(-100,100), zlim=(-100,0),
    phase=LithosphericPhases(Layers=[20 15 65], Phases = [1 2 3], Tlab=nothing), 
    DipAngle=30.0, T=LithosphericTemp(lbound="flux",nz=201))

@test sum(Temp[Int64(nel/2),Int64(nel/2),:]) ≈ 37359.648604722104


# Test the McKenzie thermal structure

# Create CartGrid struct
x        = LinRange(0.0,1200.0,64);
y        = LinRange(0.0,1200.0,64);
z        = LinRange(-660,50,64);
Cart     = CartData(XYZGrid(x, y, z));

# initialize phase and temperature matrix
Phase   = ones(Int32,size(Cart));
Temp    = ones(Float64,size(Cart))*1350;

# Create thermal structures
TsHC = HalfspaceCoolingTemp(Tsurface=20.0, Tmantle=1350, Age=120, Adiabat=0.4)
TsMK = McKenzie_subducting_slab(Tsurface = 20.0, Tmantle = 1350.0, v_cm_yr = 4.0, Adiabat = 0.0)

@test TsMK.v_cm_yr == 4.0
@test TsMK.it == 36

# Add a box with a McKenzie thermal structure

# horizontal 
Temp    = ones(Float64,size(Cart))*1350;
AddBox!(Phase, Temp, Cart; xlim=(0.0,600.0),ylim=(0.0,600.0), zlim=(-80.0, 0.0), phase = ConstantPhase(5), T=TsMK);
@test sum(Temp)  ≈ 3.518172093383281e8

# inclined slab
Temp    = ones(Float64,size(Cart))*1350;
AddBox!(Phase, Temp, Cart; xlim=(0.0,600.0),ylim=(0.0,600.0), zlim=(-80.0,0),StrikeAngle=0, DipAngle=45, phase = ConstantPhase(5), T=TsMK);
@test sum(Temp)  ≈ 3.5125017626287365e8



# horizontal slab, constant T
T_slab  = LinearWeightedTemperature(0,1,600.0,:X,ConstantTemp(1000), ConstantTemp(2000));
Temp    = ones(Float64,size(Cart))*1350;
AddBox!(Phase, Temp, Cart; xlim=(0.0,600.0),ylim=(0.0,600.0), zlim=(-80.0, 0.0), phase = ConstantPhase(5), T=T_slab);
@test sum(Temp)  ≈ 3.549127111111111e8

# horizontal slab, halfspace and McKenzie
T_slab = LinearWeightedTemperature(crit_dist=600, F1=TsHC, F2=TsMK);
Temp    = ones(Float64,size(Cart))*1350;
AddBox!(Phase, Temp, Cart; xlim=(0.0,600.0),ylim=(0.0,600.0), zlim=(-80.0, 0.0), phase = ConstantPhase(5), T=T_slab);
@test sum(Temp)  ≈ 3.499457641038468e8


Data_Final =   AddField(Cart,"Temp",Temp)

# Test the Bending slab geometry

# Create CartGrid struct
x        = LinRange(0.0,1200.0,128);
y        = LinRange(0.0,1200.0,128);
z        = LinRange(-660,50,128);
Cart     = CartData(XYZGrid(x, y, z));
X,Y,Z    = XYZGrid(x, y, z);

# initialize phase and temperature matrix
Phase   = ones(Int32,size(Cart));
Temp    = fill(1350.0,size(Cart));

t1 = Trench(Start = (400.0,400.0), End = (800.0,800.0),θ_max = 45, direction = 1.0, n_seg = 50, L0 = 600.0, D0 = 80.0, Lb = 500.0,d_decoupling = 100.0, type_bending =:Ribe)
@test t1.θ_max == 45.0
@test t1.D0 == 80.0
@test t1.L0 == 600.0
@test t1.Lb == 500.0

phase = LithosphericPhases(Layers=[5 7 88], Phases = [2 3 4], Tlab=nothing)
TsHC = HalfspaceCoolingTemp(Tsurface=20.0, Tmantle=1350, Age=30, Adiabat=0.4)
temp = TsHC;

addSlab!(Phase,Temp,Cart, t1, phase=phase, T = TsHC)
@test sum(Temp)  ≈ 2.7987343685251493e9

Data_Final      =   CartData(X,Y,Z,(Phase=Phase,Temp=Temp)) 
#Write_Paraview(Data_Final, "Data_Final");


Phase = ones(Int32,size(Cart));
Temp = fill(1350.0,size(Cart));
TsMK = McKenzie_subducting_slab(Tsurface = 20.0, Tmantle = 1350.0, v_cm_yr = 4.0, Adiabat = 0.0)
temp = TsMK 


Phase = ones(Int32,size(Cart));
Temp = fill(1350.0,size(Cart));
TsHC = HalfspaceCoolingTemp(Tsurface=20.0, Tmantle=1350, Age=120, Adiabat=0.4)
TsMK = McKenzie_subducting_slab(Tsurface = 20.0, Tmantle = 1350.0, v_cm_yr = 4.0, Adiabat = 0.0)
T_slab = LinearWeightedTemperature(crit_dist=600, F1=TsHC, F2=TsMK);
phase = LithosphericPhases(Layers=[5 7 88], Phases = [2 3 4], Tlab=nothing)
t1 = Trench(Start = (400.0,400.0), End = (800.0,800.0),θ_max = 90.0, direction = 1.0, n_seg = 50, L0 = 600.0, D0 = 80.0, Lb = 500.0,d_decoupling = 100.0, type_bending =:Ribe)

addSlab!(Phase,Temp,Cart, t1, phase=phase, T = TsHC)
@test sum(Temp)  ≈ 2.7836771215872355e9

Data_Final      =   CartData(X,Y,Z,(Phase=Phase,Temp=Temp)) 
