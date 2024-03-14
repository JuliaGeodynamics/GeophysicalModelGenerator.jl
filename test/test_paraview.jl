# this tests creating paraview output from the data
using Test
using GeophysicalModelGenerator

@testset "Paraview" begin

# Generate a 3D grid
Lon,Lat,Depth   =   lonlatdepth_grid(10:20,30:40,(-300:25:0)km);
Data            =   Depth*2; # some data
Data_set        =   GeoData(Lon,Lat,Depth,(Depthdata=Data,LonData=Lon))  
@test write_Paraview(Data_set, "test_depth3D") == nothing

# Horizontal profile @ 10 km height
Lon,Lat,Depth     =   lonlatdepth_grid(10:20,30:40,10km);
Depth[2:4,2:4,1] .= 25km     # add some fake topography

Data_set2       =   GeoData(Lon,Lat,Depth,(Topography=Depth,))  
@test write_Paraview(Data_set2, "test2")    ==  nothing

# Cross sections
Lon,Lat,Depth   =   lonlatdepth_grid(10:20,35,(-300:25:0)km);
Data_set3       =   GeoData(Lon,Lat,Depth,(DataSet=Depth,))  
@test write_Paraview(Data_set3, "test3")   ==  nothing


Lon,Lat,Depth   =   lonlatdepth_grid(15,30:40,(-300:25:0)km);
Data_set4       =   GeoData(Lon,Lat,Depth,(DataSet=Depth,))  
@test write_Paraview(Data_set4, "test4")    ==  nothing

Lon,Lat,Depth   =   lonlatdepth_grid(15,35,(-300:25:0)km);
Data_set5       =   GeoData(Lon,Lat,Depth,(DataSet=Depth,))  
@test write_Paraview(Data_set5, "test5")   ==  nothing


# Test saving vectors 
Lon,Lat,Depth           =   lonlatdepth_grid(10:20,30:40,50km);
Ve                      =   zeros(size(Depth)) .+ 1.0;
Vn                      =   zeros(size(Depth));
Vz                      =   zeros(size(Depth));
Velocity                =   (copy(Ve),copy(Vn),copy(Vz))              # tuple with 3 values, which 
Data_set_vel            =   GeoData(Lon,Lat,Depth,(Velocity=Velocity, Veast=Velocity[1]*cm/yr, Vnorth=Velocity[2]*cm/yr, Vup=Velocity[3]*cm/yr))  
@test  write_Paraview(Data_set_vel, "test_Vel") == nothing

# Test saving colors
red                     = zeros(size(Lon)); 
green                   = zeros(size(Lon)); 
blue                    = zeros(size(Lon)); 
Data_set_color          =   GeoData(Lon, Lat, Depth, (Velocity=Velocity,colors=(red,green,blue),color2=(red,green,blue)))
@test write_Paraview(Data_set_color, "test_Color") == nothing

# Manually test the in-place conversion from spherical -> cartesian (done automatically when converting GeoData->ParaviewData  )
Vel_Cart                =   (copy(Ve),copy(Vn),copy(Vz)) 
velocity_spherical_to_cartesian!(Data_set_vel, Vel_Cart);
@test Vel_Cart[2][15]   ≈   0.9743700647852352
@test Vel_Cart[1][15]   ≈   -0.224951054343865
@test Vel_Cart[3][15]   ≈   0.0

# Test saving unstructured point data (EQ, or GPS points)
Data_set_VelPoints          =       GeoData(Lon[:],Lat[:],ustrip.(Depth[:]),(Velocity=(copy(Ve[:]),copy(Vn[:]),copy(Vz[:])), Veast=Ve[:]*mm/yr, Vnorth=Vn[:]*cm/yr, Vup=Vz[:]*cm/yr))  
@test write_Paraview(Data_set_VelPoints, "test_Vel_points", PointsData=true)    ==      nothing

end