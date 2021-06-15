# This shows how to read in moho topography, plot the data as points and fit a surface through the data

# Spada, M., Bianchi, I., Kissling, E., Agostinetti, N.P., Wiemer, S., 2013. *Combining controlled-source seismology and receiver function information to derive 3-D Moho topography for Italy.* Geophysical Journal International 194, 1050–1068. [doi:10.1093/gji/ggt148](https://doi.org/10.1093/gji/ggt148)
#
# We have uploaded it here: https://seafile.rlp.net/d/a50881f45aa34cdeb3c0/
#
# Please download the files `Moho_Map_Data-WesternAlps-SpadaETAL2013_Moho1.txt`, `Moho_Map_Data-WesternAlps-SpadaETAL2013_Moho2.txt` and `Moho_Map_Data-WesternAlps-SpadaETAL2013_Moho3.txt`

using DelimitedFiles, GeophysicalModelGenerator

# read European Moho
data                =   readdlm("Moho_Map_Data-WesternAlps-SpadaETAL2013_Moho1.txt",' ',Float64,'\n', skipstart=38,header=false)
lon, lat, depth     =   data[:,1], data[:,2], -data[:,3];
data_Moho1          =   GeoData(lon,lat,depth,(MohoDepth=depth*km,))

# Write as data points:
Write_Paraview(data_Moho1, "Spada_Moho_Europe", PointsData=true) 

# Do the same with the other Moho's:
data =readdlm("Moho_Map_Data-WesternAlps-SpadaETAL2013_Moho2.txt",' ',Float64,'\n', skipstart=38,header=false);
lon, lat, depth        = data[:,1], data[:,2], -data[:,3];
Tutorial_MohoSpada_LonLat_Paraview = GeoData(lon,lat,depth,(MohoDepth=depth*km,))
Write_Paraview(data_Moho2, "Spada_Moho_Adria", PointsData=true) 
data =readdlm("Moho_Map_Data-WesternAlps-SpadaETAL2013_Moho3.txt",' ',Float64,'\n', skipstart=38,header=false);
lon, lat, depth        = data[:,1], data[:,2], -data[:,3];
data_Moho3 = GeoData(lon,lat,depth,(MohoDepth=depth*km,))
Write_Paraview(data_Moho3, "Spada_Moho_Tyrrhenia", PointsData=true) 

# Combine all data points into one data set:
lon   = [data_Moho1.lon.val;   data_Moho2.lon.val;   data_Moho3.lon.val];
lat   = [data_Moho1.lat.val;   data_Moho2.lat.val;   data_Moho3.lat.val];
depth = [data_Moho1.depth.val; data_Moho2.depth.val; data_Moho3.depth.val];
data_Moho_combined = GeoData(lon, lat, depth, (MohoDepth=depth*km,))

# Define regular grid
Lon,Lat,Depth     =   LonLatDepthGrid(4.1:0.1:11.9,42.5:.1:49,-30km);

# Interpolate using NN:
using NearestNeighbors
kdtree = KDTree([lon'; lat'])
idxs, dists = knn(kdtree, [Lon[:]'; Lat[:]'], 1, true)
Depth = zeros(size(Lon))*km;
for i=1:length(idxs)
    Depth[i] = depth[idxs[i]][1]
end
data_Moho = GeoData(Lon, Lat, Depth, (MohoDepth=Depth,))
Write_Paraview(data_Moho, "Spada_Moho_combined") 

