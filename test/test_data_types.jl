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

@test Data_cart.x[3] ≈ 6181.689604591349
@test Data_cart.y[3] ≈ 323.9686243936937
@test Data_cart.z[3] ≈ 1429.1140482465744

# Create Lon/Lat/Depth grids from given numbers or 1D vectors
Lon,Lat,Depth =  LonLatDepthGrid(10:20,30:40,(-10:-1)km);
@test size(Lon)==(11, 11, 10)
@test Lat[2,2,2]==31.0

Lon,Lat,Depth =  LonLatDepthGrid(10:20,30:40,-50km);
@test Lon[2,2]==11.0

Lon,Lat,Depth =  LonLatDepthGrid(10,30,(-10:-1)km); # 1D line @ given lon/lat
@test size(Lon)==(1,1,10)
@test Lat[2]==30.0

# throw an error if a 2D array is passed as input
@test_throws ErrorException LonLatDepthGrid(10:20,30:40,[20 30; 40 50]);



# Create 3D arrays & convert them 
Lon,Lat,Depth   =  LonLatDepthGrid(10:20,30:40,(-10:-1)km); # 3D grid
Data            =   ustrip(Depth);

Data_set1       =   GeoData(Lat,Lon,Depth,(FakeData=Data,Data2=Data.+1.))  
@test size(Data_set1.depth)==(11, 11, 10)
@test Data_set1.depth[1,2,3]==-8.0km

# Create 2D arrays & convert them 
Lon,Lat,Depth   =   LonLatDepthGrid(10:20,30:40,-50km);
Data            =   ustrip(Depth);
Data_set2       =   GeoData(Lat,Lon,Depth,(FakeData=Data,Data2=Data.+1.))  
@test Data_set2.depth[2,2]== -50.0km

# Convert the 2D and 3D arrays to their cartesian counterparts
Data_cart1      = convert(CartData,Data_set1)
@test size(Data_cart1.z)==(11, 11, 10)
@test Data_cart1.z[2,2,2] ≈ 1213.926828585578

Data_cart2      = convert(CartData,Data_set2)
@test size(Data_cart2.z)==(11, 11, 1)
@test Data_cart2.z[2,2] ≈ 1206.1036597751397
