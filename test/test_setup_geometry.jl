# test setting geometries in the different grid types
using Test, GeophysicalModelGenerator


# GeoData
Lon3D,Lat3D,Depth3D =   LonLatDepthGrid(1.0:1:10.0, 11.0:1:20.0, (-20:1:-10)*km);
Data                =   zeros(size(Lon3D));
Temp                =   ones(Float64, size(Data))*1350;
Phases              =   zeros(Int32,   size(Data));
Grid                =   GeoData(Lon3D,Lat3D,Depth3D,(DataFieldName=Data,))   

AddBox!(Phases,Temp,Grid, xlim=(2,4), zlim=(-15,-10), phase=ConstantPhase(3), DipAngle=10, T=LinearTemp(Tbot=1350, Ttop=200))
@test sum(Temp) ≈ 1.4254722397218528e6

AddEllipsoid!(Phases,Temp,Grid, cen=(4,15,-17), axes=(1,2,3), StrikeAngle=90, DipAngle=45, phase=ConstantPhase(2), T=ConstantTemp(600))
@test sum(Temp) ≈ 1.4037222397218528e6

# CartData 
X,Y,Z               =   XYZGrid(1.0:1:10.0, 11.0:1:20.0, -20:1:-10);
Data                =   zeros(size(X));
Temp                =   ones(Float64, size(Data))*1350;
Phases              =   zeros(Int32,   size(Data));
Grid                =   CartData(X,Y,Z,(DataFieldName=Data,))   

AddBox!(Phases,Temp,Grid, xlim=(2,4), zlim=(-15,-10), phase=ConstantPhase(3), DipAngle=10, T=LinearTemp(Tbot=1350, Ttop=200))
@test sum(Temp) ≈ 1.4254722397218528e6

AddEllipsoid!(Phases,Temp,Grid, cen=(4,15,-17), axes=(1,2,3), StrikeAngle=90, DipAngle=45, phase=ConstantPhase(2), T=ConstantTemp(600))
@test sum(Temp) ≈ 1.4037222397218528e6

# CartGrid
Grid                =   CreateCartGrid(size=(10,20,30),x=(0.,10), y=(0.,10), z=(2.,10))
Temp                =   ones(Float64, Grid.N...)*1350;
Phases              =   zeros(Int32,  Grid.N...);

AddBox!(Phases,Temp,Grid, xlim=(2,4), zlim=(4,8), phase=ConstantPhase(3), DipAngle=10, T=LinearTemp(Tbot=1350, Ttop=200))

@test maximum(Phases) == 3

# Create a CartData structure from it
Data = CartData(Grid, (T=Temp, Phases=Phases))

@test NumValue(Data.x[3,3,2]) ≈ 2.2222222222222223

# To be tested: doing the same for cross-sections or horizontal slices
Grid2D              =   CreateCartGrid(size=(10,30),x=(0.,10), z=(2.,10))