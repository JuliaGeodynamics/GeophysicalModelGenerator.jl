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
#=
# Take topography of the Alps, flatten it, shift it & transfer to km
using GeophysicalModelGenerator, GMT

# Topo in lon/lat/depth
Topo            =   ImportTopo([4,20,37,49], file="@earth_relief_01m.grd")


Topo_UTM        =   Convert2UTMzone(Topo,33,true);
centerUTM       =   (mean(Topo_UTM.EW.val), mean(Topo_UTM.NS.val))
Topo_Cart       =   Convert2CartData(Topo_UTM,center=centerUTM)

# New grid on which we want to interpolate the data
X,Y,Z           =   LonLatDepthGrid(-500:10:500, -600:10:600,0);
Data_surf       =   CartData(X,Y,Z,(Data=Z,))



# Convert the cartesian grid to lon/lat (in fact we can do that with a grid w/out data sets):
Data_UTM        = Convert2UTMzone(Data_surf, 33, true, center=centerUTM);
Data_lonlat     = convert(GeoData,Data_UTM)   

fields_new      = InterpolateDataFields2D(Topo,Data_lonlat.lon.val, Data_lonlat.lat.val)
Data_surf       = CartData(Data_surf.x.val, Data_surf.y.val, fields_new.Topography, (Topography=fields_new.Topography,))





Data_surf1       =   CartData(Data_surf.x.val,Data_surf.y.val,Z2d,(Data=Z2d,))
Write_Paraview(Data_surf1,"Data_surf1")

=#

Write_Paraview(Data_set3D_transformed, "data_3d_transformed")