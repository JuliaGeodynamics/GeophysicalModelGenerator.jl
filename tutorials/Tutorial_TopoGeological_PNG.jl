"""
This is Tutorial_TopoGeological_PNG.jl. It contains all the necessary commands to import the netCDF file from ETOPO1 and 
to convert it to the GeoData format. The, a png file with the geological map of the Alps. The following steps are performed:
1. Read data from this file.
2. Select a longitude and latitude range.
3. Put the data in a GeoData format (this is the format that is used internally in the GMG).
4. Export that data to a format readable by Paraview.

ETOPO1 data files used in this example can be downloaded here:  
[https://ngdc.noaa.gov/mgg/global/global.html](https://ngdc.noaa.gov/mgg/global/global.html)
"""

using GMT, NearestNeighbors, GeoParams, GeophysicalModelGenerator

# 1. define where the file is located on your computer
filename_topo = "./ETOPO1/ETOPO1_Ice_g_gmt4.grd" # this only works if you are in the directory where the data is located
filename_geo  = "./tectonic_maps_4dmb_2020_09_17/tectonicmap_SPP.png" #tectonicmap_schmid_2004.png"

# 2. Select latitude and longitude range for topo
lat_min = 37.0
lat_max = 49.0
lon_min = 4.0
lon_max = 20.0

# 2.1 load desired topography data using GMT.jl and adapt the data to make it compatible with GeoData
G = gmtread(filename_topo, limits=[lon_min,lon_max,lat_min,lat_max]);
Lon,Lat,Depth    =   LonLatDepthGrid(G.x[1:end],G.y[1:end],0);
numel_topo       =   prod(size(Lon));
Depth[:,:,1]     =   1e-3*G.z';
DataTopo         =   GeophysicalModelGenerator.GeoData(Lon, Lat, Depth, (Topography=Depth*km,))

# 3. load data for the geological map
Corner_LowerLeft    =   ( 5.0, 43.6 , 0.0)
Corner_UpperRight   =   (17.0, 48.5 , 0.0)

Corner_LowerLeft    =   ( lon_min, lat_min , 0.0)
Corner_UpperRight   =   (lon_max, lat_max , 0.0)

DataPNG = Screenshot_To_GeoData(filename_geo, Corner_LowerLeft, Corner_UpperRight)

# 4. interpolate geological map data colors onto topo grid using nearest neighbor interpolation

# 4.1 set up the tree and determine nearest neighbors 
coord = [vec(DataPNG.lon.val)';vec(DataPNG.lat.val)'];
kdtree = KDTree(coord; leafsize = 10);
points = [vec(Lon)';vec(Lat)'];
idx,dist = nn(kdtree, points);

# 4.2 assign color data to topo grid
red   = zeros(size(Depth));
green = zeros(size(Depth));
blue  = zeros(size(Depth));

tmp                = DataPNG.fields.colors[1];
red[1:numel_topo] = tmp[idx];
tmp                = DataPNG.fields.colors[2];
green[1:numel_topo] = tmp[idx];
tmp                = DataPNG.fields.colors[3];
blue[1:numel_topo] = tmp[idx];

# 4.3 set colors outside the box to white
ind_tmp = Lat .<  Corner_LowerLeft[2];
red[ind_tmp] .= 1;
green[ind_tmp] .= 1;
blue[ind_tmp] .= 1;

ind_tmp = Lat .> Corner_UpperRight[2];
red[ind_tmp] .= 1;
green[ind_tmp] .= 1;
blue[ind_tmp] .= 1;

ind_tmp = Lon .<  Corner_LowerLeft[1];
red[ind_tmp] .= 1;
green[ind_tmp] .= 1;
blue[ind_tmp] .= 1;

ind_tmp = Lon .> Corner_UpperRight[1];
red[ind_tmp] .= 1;
green[ind_tmp] .= 1;
blue[ind_tmp] .= 1;

# Set up the Data structure
Data_set        =   GeoData(Lon, Lat, Depth, (Topography=Depth*km,colors=(red,green,blue)))

# Export the data structure to Paraview format
Write_Paraview(Data_set, "test_GeoMap")
