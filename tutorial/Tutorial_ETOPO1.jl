"""
This is Tutorial_ETOPO1.jl. It contains all the necessary commands to import the netCDF file from ETOPO1, 
to convert it to the GeoData format and to export it to a Paraview format. The following steps are performed:
1. Read data from this file.
2. Select a longitude and latitude range.
3. Put the data in a GeoData format (this is the format that is used internally in the GMG).
4. Export that data to a format readable by Paraview.

ETOPO1 data files used in this example can be downloaded here:  
[https://ngdc.noaa.gov/mgg/global/global.html](https://ngdc.noaa.gov/mgg/global/global.html)

doi:10.7289/V5C8276M
"""

using NetCDF, GeophysicalModelGenerator

# 1. define where the file is located on your computer
filename = "./ETOPO1_Ice_g_gmt4.grd" # this only works if you are in the directory where the data is located

# 2. load desired data
lat = ncread(filename,"y")
lon = ncread(filename,"x")
topo = ncread(filename,"depth")

# 3. Select latitude and longitude range
lat_min = 
lat_max
lon_min
lon_max


# For netCDF data, 3D coordinates of a regular grid are only given as 1D vectors. As we need to compute Cartesian coordinates for
# Paraview, we need the full matrix of grid coordinates
Lon3D,Lat3D,Depth3D = LonLatDepthGrid(lon, lat, 0km);

# Set up the Data structure
Data_set1       =   GeoData(Lat3D,Lon3D,Depth3D,(VS=vs,))

# Export the data structure to Paraview format
Write_Paraview(Data_set, "test_netcdf_3D")
