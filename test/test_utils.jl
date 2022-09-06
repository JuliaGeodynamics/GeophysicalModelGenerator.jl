# Tests data_import.jl

using Test
using GeophysicalModelGenerator

# should throw an error with a 2D dataset
Lon,Lat,Depth   =   LonLatDepthGrid(10:20,30:40,-50km);
Data1           =   Depth*2;                # some data
Vx1,Vy1,Vz1     =   Data1*3,Data1*4,Data1*5
Data_set2D      =   GeoData(Lon,Lat,Depth,(Depthdata=Data1,LonData1=Lon, Velocity=(Vx1,Vy1,Vz1)))
@test_throws ErrorException CrossSection(Data_set2D, Depth_level=-10)

# Create 3D volume with some fake data
Lon,Lat,Depth   =   LonLatDepthGrid(10:20,30:40,(-300:25:0)km);
Data            =   Depth*2;                # some data
Vx,Vy,Vz        =   ustrip(Data*3),ustrip(Data*4),ustrip(Data*5);
Data_set3D      =   GeoData(Lon,Lat,Depth,(Depthdata=Data,LonData=Lon, Velocity=(Vx,Vy,Vz)))

# Create 3D cartesian dataset
Data_setCart3D  =   CartData(Lon,Lat,Depth,(Depthdata=Data,LonData=Lon, Velocity=(Vx,Vy,Vz)))

# Create 3D volume with some fake data
Lon,Lat,Depth           =   LonLatDepthGrid(10:20,30:40,(0:-25:-300)km);
Data                    =   Depth*2;                # some data
Vx,Vy,Vz                =   ustrip(Data*3),ustrip(Data*4),ustrip(Data*5);
Data_set3D_reverse      =   GeoData(Lon,Lat,Depth,(Depthdata=Data,LonData=Lon, Velocity=(Vx,Vy,Vz)))

# Create cross-sections in various directions (no interpolation which is default)
test_cross      =   CrossSection(Data_set3D, Depth_level=-100km)
@test test_cross.fields[1][41]==-200km
@test test_cross.fields[2][31]==18
@test test_cross.fields[3][1][30]==-600
@test test_cross.fields[3][2][30]==-800
@test test_cross.fields[3][3][30]==-1000

# throw error if outside bounds
@test_throws ErrorException CrossSection(Data_set3D, Depth_level=100km)

test_cross      =   CrossSection(Data_set3D, Lon_level=15)
@test test_cross.fields[1][41]==-450km
@test test_cross.fields[2][31]==15
@test test_cross.fields[3][1][30]==-1500
@test test_cross.fields[3][2][30]==-2000
@test test_cross.fields[3][3][30]==-2500

test_cross      =   CrossSection(Data_set3D, Lat_level=35)
@test test_cross.fields[1][41]==-450km
@test test_cross.fields[2][31]==18
@test test_cross.fields[3][1][30]==-1500
@test test_cross.fields[3][2][30]==-2000
@test test_cross.fields[3][3][30]==-2500

# Create cross-sections with interpolation in various directions
test_cross      =   CrossSection(Data_set3D, Depth_level=-100km, dims=(50,100), Interpolate=true)
@test size(test_cross.fields[1])    ==  (50,100,1)
@test size(test_cross.fields[3][2]) ==  (50,100,1)

test_cross      =   CrossSection(Data_set3D, Lon_level=15, dims=(50,100), Interpolate=true)
@test size(test_cross.fields[3][2])==(1,50,100)
@test Write_Paraview(test_cross, "profile_test")==nothing

test_cross      =   CrossSection(Data_set3D, Lat_level=35, dims=(50,100), Interpolate=true)
@test size(test_cross.fields[3][2])==(50,1,100)

# Diagonal cross-section
test_cross      =   CrossSection(Data_set3D, Start=(10,30), End=(20,40), dims=(50,100), Interpolate=true)
@test size(test_cross.fields[3][2])==(50,100,1)
@test Write_Paraview(test_cross, "profile_test")==nothing

#test_cross_rev  =   CrossSection(Data_set3D_reverse, Start=(10,30), End=(20,40), dims=(50,100), Interpolate=true)
#@test size(test_cross_rev.fields[3][2])==(50,100,1)
#@test Write_Paraview(test_cross_rev, "profile_test_rev")[1]=="profile_test_rev.vts"

# Cross-section with cartesian data
test_cross      =   CrossSection(Data_setCart3D, Lon_level=15, dims=(50,100), Interpolate=true)
@test size(test_cross.fields[3][2])==(1,50,100)
@test test_cross.x[1,2,3]==GeoUnit(15km)

# Extract sub-volume

# with interpolation
Data_sub_Interp = ExtractSubvolume(Data_set3D,Lon_level=(10,15), Lat_level=(30,32), Interpolate=true, dims=(51,21,32))
@test Data_sub_Interp.fields[1][11]==-600km
@test size(Data_sub_Interp.lat)==(51,21,32)

# no interpolation
Data_sub_NoInterp = ExtractSubvolume(Data_set3D,Lon_level=(10,15), Lat_level=(30,32), Interpolate=false, dims=(51,21,32))
@test Data_sub_NoInterp.fields[1][11]==-600km
@test size(Data_sub_NoInterp.lat)==(6,3,13)

# Extract subset of cross-section
test_cross      =   CrossSection(Data_set3D, Lat_level=35, dims=(50,100), Interpolate=true)
Data_sub_cross  =   ExtractSubvolume(test_cross, Depth_level=(-100km,0km), Interpolate=false)
@test Data_sub_cross.fields[1][11]==-200.00000000000003km
@test size(Data_sub_cross.lat)==(50,1,34)

# compute the mean velocity per depth in a 3D dataset and subtract the mean from the given velocities
Data_pert   =   SubtractHorizontalMean(ustrip(Data))    # 3D, no units
@test Data_pert[10] == 0.0

Data_pert   =   SubtractHorizontalMean(Data)            # 3D with units
@test Data_pert[10] == 0.0km

Data_pert   =   SubtractHorizontalMean(Data, Percentage=true)            # 3D with units
@test Data_pert[1000] == 0.0

Data2D = Data[:,1,:];
Data_pert   =   SubtractHorizontalMean(Data2D, Percentage=true)         # 2D version with units [dp the same along a vertical profile]

Data_set2D  =   GeoData(Lon,Lat,Depth,(Depthdata=Data,LonData=Lon,Pertdata=Data_pert ,Velocity=(Vx,Vy,Vz)))
@test Data_set2D.fields[3][10,8,1] == 0


# Create surface ("Moho")
Lon,Lat,Depth   =   LonLatDepthGrid(10:20,30:40,-40km);
Depth           =   Depth + Lon*km;     # some fake topography on Moho
Data_Moho       =   GeoData(Lon,Lat,Depth,(MohoDepth=Depth,LonData=Lon))


# Test intersecting a surface with 2D or 3D data sets
Above       =   AboveSurface(Data_set3D, Data_Moho);            # 3D regular ordering
@test Above[1,1,12]==true
@test Above[1,1,11]==false

Above       =   AboveSurface(Data_set3D_reverse, Data_Moho);    #  3D reverse depth ordering
@test Above[1,1,2]==true
@test Above[1,1,3]==false

Above       =   AboveSurface(Data_sub_cross, Data_Moho);        # 2D cross-section
@test Above[end]==true
@test Above[1]==false


# Test VoteMaps
Data_VoteMap = VoteMap(Data_set3D, "Depthdata<-560",dims=(10,10,10))
@test Data_VoteMap.fields[:VoteMap][101]==0
@test Data_VoteMap.fields[:VoteMap][100]==1

Data_VoteMap = VoteMap(Data_set3D_reverse, "Depthdata<-560",dims=(10,10,10))
@test Data_VoteMap.fields[:VoteMap][101]==0
@test Data_VoteMap.fields[:VoteMap][100]==1

# Combine 2 datasets
Data_VoteMap = VoteMap([Data_set3D_reverse, Data_set3D], ["Depthdata<-560","LonData>19"],dims=(10,10,10))
@test Data_VoteMap.fields[:VoteMap][10,9,1]==2
@test Data_VoteMap.fields[:VoteMap][9 ,9,1]==1
@test Data_VoteMap.fields[:VoteMap][9 ,9,2]==0

# Test rotation routines
X,Y,Z   =   LonLatDepthGrid(10:20,30:40,-50:-10);
Data_C  =   ParaviewData(X,Y,Z,(Depth=Z,))
Data_C1 =   RotateTranslateScale(Data_C, Rotate=30);
@test Data_C1.x.val[10] ≈ 20.964101615137753
@test Data_C1.y.val[10] ≈ 32.66987298107781
@test Data_C1.z.val[20] == -50

Data_C1 = RotateTranslateScale(Data_C, Scale=10, Rotate=10, Translate=(1,2,3));
@test Data_C1.x.val[10] ≈ 213.78115820908607
@test Data_C1.y.val[10] ≈ 339.4092822315127
@test Data_C1.z.val[20] == -497.0


# Test
