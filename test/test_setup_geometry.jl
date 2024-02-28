# test setting geometries in the different grid types
using Test, GeophysicalModelGenerator


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


# test AboveSurface with the Grid object
Grid        =   CreateCartGrid(size=(10,20,30),x=(0.,10), y=(0.,10), z=(-10.,2.))
@test Grid.Δ[2] ≈ 0.5263157894736842

Temp        =   ones(Float64, Grid.N...)*1350;
Phases      =   zeros(Int32,  Grid.N...);


Topo_cart   =   CartData(XYZGrid(-1:.2:20,-12:.2:13,0));
ind         =   AboveSurface(Grid, Topo_cart);
@test sum(ind[1,1,:]) == 5

ind         =   BelowSurface(Grid, Topo_cart);
@test sum(ind[1,1,:]) == 25

# Create Grid & nondimensionalize it
CharDim     =   GEO_units();
Grid        =   CreateCartGrid(size=(10,20,30),x=(0.0km,10km), y=(0.0km, 10km), z=(-10.0km, 2.0km), CharDim=CharDim)
@test Grid.Δ[2] ≈ 0.0005263157894736842


# test 1D-explicit thermal solver for AddBox
Grid                =   CreateCartGrid(size=(96,96,96),x=(-100,100), y=(-100,100), z=(-400,0))
Temp                =   ones(Float64, Grid.N...)*1350;
Phases              =   zeros(Int64,  Grid.N...);

AddBox!(Phases,Temp,Grid, xlim=(-100,100), zlim=(-200,0), Origin=(0.0,0.0,0.0),
    phase=LithosphericPhases(Layers=[20 15 165], Phases = [1 2 3], Tlab=nothing), 
    DipAngle=0.0, T=LithosphericTemp(nz=301))

