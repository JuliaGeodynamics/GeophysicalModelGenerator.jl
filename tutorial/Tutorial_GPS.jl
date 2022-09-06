# This reads the GPS data of Sanchez et al. (2018) https://essd.copernicus.org/articles/10/1503/2018/#section7
#
# They can be downloaded from: https://doi.pangaea.de/10.1594/PANGAEA.886889
#
using GeophysicalModelGenerator
using DataFrames, CSV

# Read in coordinates of the grids (not stations as they are given in a different reference frame)
#
# Note: the Vz velocity is given on a fully regular grid; yet Ve/Vn only on on-land stations which makes this
# a bit tricky.
#
# The approach we take here is to first read in the Vz data points & reshape it to 2D matrixes
# Next, we read the Ve/Vn data and add them to the Vz grid
data_file               =   CSV.File("ALPS2017_DEF_VT.GRD",datarow=18,header=false,delim=' ')

num_columns             =   4;
data                    =   ParseColumns_CSV_File(data_file, num_columns);     # Read numerical data from the file

# Reshape data to 2D (and 3D) matrixes
lon_Vz, lat_Vz, Vz_vec  =   data[:,1], data[:,2], data[:,3];
lon, lat, Vz            =   zeros(41,31,1),zeros(41,31,1),zeros(41,31,1)
lon[:,:,1]              =   reshape(lon_Vz,(41,31))
lat[:,:,1]              =   reshape(lat_Vz,(41,31))
Vz[:,:,1]               =   reshape(Vz_vec,(41,31))

Vz                      =   Vz*1000;                #in mm/year (original data in m/yr)
Ve                      =   ones(size(Vz))*NaN;
Vn                      =   ones(size(Vz))*NaN;

# Now read the Ve/Vn (horizontal) velocities:
data_file                       =   CSV.File("ALPS2017_DEF_HZ.GRD",datarow=18,header=false,delim=' ')
data                            =   ParseColumns_CSV_File(data_file, 10)
lon_Hz, lat_Hz, Ve_Hz, Vn_Hz    =   data[:,1], data[:,2], data[:,3],  data[:,4];
Ve_Hz                   =   Ve_Hz*1000; #in mm/year
Vn_Hz                   =   Vn_Hz*1000; #in mm/year

for i in eachindex(lon_Hz)
    ind = intersect(findall(x->x==lon_Hz[i], lon), findall(x->x==lat_Hz[i], lat))
    Ve[ind] .= Ve_Hz[i];
    Vn[ind] .= Vn_Hz[i];
end

Vmagnitude          =   sqrt.(Ve.^2 + Vn.^2 + Vz.^2);  # velocity magnitude in mm/yr


# Finally, it would be nice to put the data points on the topography.
# The elevation of the points is not given in the GPS dataset, so we use GMT to extract a topographic grid
# and interpolate the elevation on the GPS grid locations
using GMT, Interpolations
Elevation               =   gmtread("@earth_relief_01m.grd", limits=[3,17,42,50]);
Lon_Topo,Lat_Topo,dummy =   LonLatDepthGrid(Elevation.x[1:end-1],Elevation.y[1:end-1],0);


Lon_vec     =  Lon_Topo[:,1,1];
Lat_vec     =  Lat_Topo[1,:,1];
interpol    =   LinearInterpolation((Lon_vec, Lat_vec), Elevation.z');      # create interpolation object
height      =   interpol.(lon,lat)/1e3;


# At this stage we have lon/lat/height of all points as well as velocity components

GPS_Sanchez_grid        =   GeoData(lon,lat,height,(Velocity_mm_year=(Ve,Vn,Vz),V_north=Vn*mm/yr, V_east=Ve*mm/yr, V_vertical=Vz*mm/yr, Vmagnitude = Vmagnitude*mm/yr, Topography = height*km))

# Save paraview
Write_Paraview(GPS_Sanchez_grid, "GPSAlps_Sanchez_2017_grid")
