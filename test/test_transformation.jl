# tests transformations from GeoData <=> Cartesian
using Test
using GeophysicalModelGenerator

# Create 3D volume with some fake data
Lon,Lat,Depth   =   lonlatdepthGrid(5:25,20:50,(-1300:100:0)km);
Data_set3D      =   GeoData(Lon,Lat,Depth,(Depthdata=Depth*2 + Lon*km,LonData=Lon))  

proj            =   ProjectionPoint(Lon=20,Lat=35)

# Convert this 3D dataset to a Cartesian dataset (the grid will not be orthogonal)
Data_set3D_Cart =   convert2CartData(Data_set3D, proj)  
@test sum(abs.(Value(Data_set3D_Cart.x))) ≈ 5.293469089428514e6km

# Create Cartesian grid
X,Y,Z           =   xyzGrid(-400:100:400,-500:200:500,(-1300:100:0)km);
Data_Cart       =   CartData(X,Y,Z,(Z=Z,))  

# Project values of Data_set3D to the cartesian data
Data_Cart       =   projectCartData(Data_Cart, Data_set3D, proj)
@test sum(Data_Cart.fields.Depthdata) ≈ -967680.9136292854km
#@test sum(Data_Cart.fields.Depthdata) ≈ -1.416834287168597e6km
@test sum(Data_Cart.fields.LonData) ≈ 15119.086370714615


# Next, 3D surface (like topography)
Lon,Lat,Depth   =   lonlatdepthGrid(5:25,20:50,0);
Depth           =   cos.(Lon/5).*sin.(Lat)*10;
Data_surf       =   GeoData(Lon,Lat,Depth,(Z=Depth,));  
Data_surf_Cart  =   convert2CartData(Data_surf, proj); 

# Cartesian surface
X,Y,Z           =   xyzGrid(-500:10:500,-900:20:900,0);
Data_Cart       =   CartData(X,Y,Z,(Z=Z,))  

Data_Cart       =   projectCartData(Data_Cart, Data_surf, proj)
@test sum(Value(Data_Cart.z)) ≈ 1858.2487019158766km
@test sum(Data_Cart.fields.Z) ≈ 1858.2487019158766


# Cartesian surface when UTM data is used
WE,SN,depth     =   xyzGrid(420000:1000:430000, 4510000:1000:4520000, 0);

Data_surfUTM    =   UTMData(WE, SN, depth, 33, true, (Depth = WE,));
Data_Cart       =   CartData(X,Y,Z,(Z=Z,))  
Data_Cart       =   projectCartData(Data_Cart, Data_surfUTM, proj)

@test sum(Value(Data_Cart.z)) ≈ 0.0km
@test sum(Data_Cart.fields.Depth) ≈ 3.9046959539921126e9
