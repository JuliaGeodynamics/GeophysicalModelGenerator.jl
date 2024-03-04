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