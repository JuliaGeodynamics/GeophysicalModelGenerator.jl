# tests transformations from GeoData <=> Cartesian
using Test
using GeophysicalModelGenerator

# Create 3D volume with some fake data
Lon,Lat,Depth   =   LonLatDepthGrid(10:20,30:40,(-300:25:0)km);
Data_set3D      =   GeoData(Lon,Lat,Depth,(Depthdata=Depth*2,LonData=Lon))  



# Transform the grid to cartesian & back, without flattening (so no distortion)
data_Cart, refPoint     = GeoData_To_Cartesian(Data_set3D, referencePoint="CenterBottom", Flatten=false);
Data_set3D_transformed  = Cartesian_To_GeoData(data_Cart, refPoint, Flatten=false) # transform back

@test Data_set3D_transformed.lon[1] ≈ Data_set3D.lon[1]
@test Data_set3D_transformed.lat[1] ≈ Data_set3D.lat[1]
@test Data_set3D_transformed.depth[1] ≈ Data_set3D.depth[1]

@test Data_set3D_transformed.lon[600] ≈ Data_set3D.lon[600]
@test Data_set3D_transformed.lat[600] ≈ Data_set3D.lat[600]
@test Data_set3D_transformed.depth[600] ≈ Data_set3D.depth[600]

# Next, we flatten the setup. That 


Write_Paraview(Data_set3D_transformed, "data_3d_transformed")