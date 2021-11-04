using Test
using GeophysicalModelGenerator


# Create 1D dataset with lat/lon/depth
Lat         =   1.0:10.0;
Lon         =   11.0:20.0;
Depth       =   (-20:-11)*km;
Data        =   zeros(size(Lon));
Data_set    =   GeoData(Lat,Lon,Depth,(FakeData=Data,Data2=Data.+1.))     
@test Data_set.depth[2]==-19km

Depth1      =   (-20.:-11.)*m;            # depth has units of m
Data_set1   =   GeoData(Lat,Lon,Depth1,(FakeData=Data,Data2=Data.+1.))    
@test Data_set1.depth[2]==-0.019km

Depth2      =   -20.:-11.;              # no units
Data_set2   =   GeoData(Lat,Lon,Depth2,(FakeData=Data,Data2=Data.+1.))     
@test Data_set2.depth[2]==-19km

# test that it works if we give a Data array, rather than a NamedTuple, as input (we add a default name)
Data_set3 = GeoData(Lat,Lon,Depth,Data)
@test keys(Data_set3.fields)[1] == :DataSet1

# test that it works if we give a Tuple, rather than a NamedTuple, as input (we add a default name)
Data_set4 = GeoData(Lat,Lon,Depth,(Data,))
@test keys(Data_set4.fields)[1] == :DataSet1

# Throw an error if we supply a Tuple with 2 fields (the user should really supply names in that case)
@test_throws ErrorException GeoData(Lat,Lon,Depth,(Data,Data))

# check that an error is thrown if a different size input is given for depth
Depth2      =   (-100:10:1.0)               # no units
@test_throws ErrorException GeoData(Lat,Lon,Depth2,(FakeData=Data,Data2=Data.+1.)) 

# Convert 1D vector to cartesian structure
Data_cart = convert(CartData,Data_set)

@test Data_cart.x[3] ≈ 6189.685604255086
@test Data_cart.y[3] ≈ 324.3876769792181
@test Data_cart.z[3] ≈ 1421.35608984477 

# Create Lon/Lat/Depth and X/Y/Z grids from given numbers or 1D vectors
Lon,Lat,Depth =  LonLatDepthGrid(10:20,30:40,(-10:-1)km);
X,Y,Z         =  XYZGrid(10:20,30:40,(-10:-1)km);
@test size(Lon)==(11, 11, 10)
@test Lat[2,2,2]==31.0
@test size(X)==(11, 11, 10)
@test Y[2,2,2]==31.0

Lon,Lat,Depth   =  LonLatDepthGrid(10:20,30:40,-50km);
X,Y,Z           =  XYZGrid(10:20,30:40,-50km);
@test Lon[2,2] == 11.0
@test X[2,2]   == 11.0

Lon,Lat,Depth   =  LonLatDepthGrid(10,30,(-10:-1)km); # 1D line @ given lon/lat
X,Y,Z           =  XYZGrid(10,30,(-10:-1)km);
@test size(Lon)==(1,1,10)
@test Lat[2]==30.0
@test size(X)==(1,1,10)
@test Y[2]==30.0

# throw an error if a 2D array is passed as input
@test_throws ErrorException LonLatDepthGrid(10:20,30:40,[20 30; 40 50]);
@test_throws ErrorException XYZGrid(10:20,30:40,[20 30; 40 50]);


# Create 3D arrays & convert them 
Lon,Lat,Depth   =  LonLatDepthGrid(10:20,30:40,(-10:-1)km); # 3D grid
Data            =   ustrip(Depth);

Data_set1       =   GeoData(Lon,Lat,Depth,(FakeData=Data,Data2=Data.+1.))  
@test size(Data_set1.depth)==(11, 11, 10)
@test Data_set1.depth[1,2,3]==-8.0km

# double-check that giving 3D arrays in the wrong ordering triggers a warning message
@test_throws ErrorException("It appears that the lon array has a wrong ordering")    GeoData(Lat,Lon,Depth,(FakeData=Data,Data2=Data.+1.))  

# Create 2D arrays & convert them 
Lon,Lat,Depth   =   LonLatDepthGrid(10:20,30:40,-50km);
Data            =   ustrip(Depth);
Data_set2       =   GeoData(Lon,Lat,Depth,(FakeData=Data,Data2=Data.+1.))  
@test Data_set2.depth[2,2]== -50.0km

# Convert the 2D and 3D arrays to their cartesian counterparts
Data_cart1      = convert(CartData,Data_set1)
@test size(Data_cart1.z)==(11, 11, 10)
@test Data_cart1.z[2,2,2] ≈ 3261.2581739797533

Data_cart2      = convert(CartData,Data_set2)
@test size(Data_cart2.z)==(11, 11, 1)
@test Data_cart2.z[2,2] ≈ 3240.141612908441

# Create UTM Data structure
ew          =   422123.0:100:433623.0
ns          =   4.514137e6:100:4.523637e6
depth       =   -5.4:.25:0.6
EW,NS,Depth =   XYZGrid(ew, ns, depth);
Data        =   ustrip.(Depth);
Data_set    =   UTMData(EW,NS,Depth,33, true, (FakeData=Data,Data2=Data.+1.))  

@test Data_set.EW[3,4,2]==422323.0
@test Data_set.NS[3,4,2]==4.514437e6
@test Data_set.depth[3,4,2]==-5.15km
@test Data_set.northern[1] == true
@test Data_set.zone[1] == 33

# convert from UTMData -> GeoData
Data_set1 = convert(GeoData, Data_set)
@test Data_set1.lon[20] == 14.099668158564413
@test Data_set1.lat[20] == 40.77470011887963
@test Data_set1.depth[20] == -5.4km

# Convert from GeoData -> UTMData
Data_set2 = convert(UTMData, Data_set1)
@test sum(abs.(Data_set2.EW.val-Data_set.EW.val)) < 1e-5 
