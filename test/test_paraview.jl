# this tests creating paraview output from the data
using Test
using GeophysicalModelGenerator



# Generate a 3D grid
Lon,Lat,Depth   =   LonLatDepthGrid(10:20,30:40,(-300:25:0)km);
Data            =   Depth*2; # some data
Data_set        =   GeoData(Lat,Lon,Depth,(Depthdata=Data,LonData=Lon))  
outfile         =   Write_Paraview(Data_set, "test_depth3D")
@test outfile[1]    ==  "test_depth3D.vts"

# Horizontal profile @ 10 km height
Lon,Lat,Depth     =   LonLatDepthGrid(10:20,30:40,10km);
Depth[2:4,2:4,1] .= 25km     # add some fake topography

Data_set2       =   GeoData(Lat,Lon,Depth,(Topography=Depth,))  
outfile2       =   Write_Paraview(Data_set2, "test2")
@test outfile2[1]    ==  "test2.vts"

# Cross sections
Lon,Lat,Depth   =   LonLatDepthGrid(10:20,35,(-300:25:0)km);
Data_set3       =   GeoData(Lat,Lon,Depth,(DataSet=Depth,))  
outfile3        =   Write_Paraview(Data_set3, "test3")
@test outfile3[1]    ==  "test3.vts"


Lon,Lat,Depth   =   LonLatDepthGrid(15,30:40,(-300:25:0)km);
Data_set4       =   GeoData(Lat,Lon,Depth,(DataSet=Depth,))  
outfile4        =   Write_Paraview(Data_set4, "test4")
@test outfile4[1]    ==  "test4.vts"

Lon,Lat,Depth   =   LonLatDepthGrid(15,35,(-300:25:0)km);
Data_set5       =   GeoData(Lat,Lon,Depth,(DataSet=Depth,))  
outfile5        =   Write_Paraview(Data_set5, "test5")
@test outfile5[1]    ==  "test5.vts"


# To be done: cases with vectors (as we will need to do vector transformation from lat/lon)

