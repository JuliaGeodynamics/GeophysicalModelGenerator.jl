# this tests creating paraview output from the data
using Test
using GeophysicalModelGenerator



# Generate a 3D grid
Lon,Lat,Depth   =   LonLatDepthGrid(10:20,30:40,(-300:25:0)km);
Data            =   Depth*2; # some data
Data_set        =   GeoData(Lat,Lon,Depth,(Depthdata=Data,LonData=Lon))  
outfile         =   Write_Paraview(Data_set, "test_depth3D")
@test outfile[1]    ==  "test_depth3D.vts"

# Horizontal profile @ 200 km
#Lon,Lat,Depth   =   LonLatDepthGrid(10:20,30:40,-200km);
#Data_set2       =   GeoData(Lat,Lon,Depth,(LonData=Lon,))  
#outfile         =   Write_Paraview(Data_set2, "test2")






