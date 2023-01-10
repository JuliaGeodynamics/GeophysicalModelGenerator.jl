"""
This is Tutorial_netCDF.jl. It contains all the necessary commands to import a netCDF file from a seismic tomography,
to convert it to the GeoData format and to export it to a Paraview format. The following steps are performed:
1. Read data from this file.
2. Put the data in a GeoData format (this is the format that is used internally in the GMG).
3. Export that data to a format readable by Paraview.

The data file used in this example can be downloaded here:
[https://ds.iris.edu/files/products/emc/emc-files/El-Sharkawy-etal-G3.2020-MeRE2020-Mediterranean-0.0.nc](https://ds.iris.edu/files/products/emc/emc-files/El-Sharkawy-etal-G3.2020-MeRE2020-Mediterranean-0.0.nc)
"""

using NetCDF, GeophysicalModelGenerator

# 1. define where the file is located on your computer
filename = "./El-Sharkawy-etal-G3.2020-MeRE2020-Mediterranean-0.0.nc" # this only works if you are in the directory where the data is located

# 2. load desired data

# Now check with ncinfo(filename), what the variables are called exactly and what the contents of your netCDF file are

lat = ncread(filename,"latitude");
lon = ncread(filename,"longitude");
depth = ncread(filename,"depth");
vs    = ncread(filename,"Vs");
depth = -1 .* depth # CAREFUL: MAKE SURE DEPTH IS NEGATIVE, AS THIS IS THE ASSUMPTION IN GeoData

# For netCDF data, 3D coordinates of a regular grid are only given as 1D vectors. As we need to compute Cartesian coordinates for
# Paraview, we need the full matrix of grid coordinates
Lon3D,Lat3D,Depth3D = LonLatDepthGrid(lon, lat, depth);

# Set up the Data structure
Data_set       =   GeoData(Lon3D,Lat3D,Depth3D,(VS=vs,))

# Export the data structure to Paraview format
Write_Paraview(Data_set, "test_netcdf_3D")
