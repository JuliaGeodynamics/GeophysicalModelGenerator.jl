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

using GMT, GeoParams, GeophysicalModelGenerator, Shapefile, Luxor

# 1. define where the file is located on your computer
filename_topo = "./ETOPO1/ETOPO1_Ice_g_gmt4.grd" # this only works if you are in the directory where the data is located
filename_geo  = "./tectonic_maps_4dmb_2020_09_17/shape_files/tect_units_alcapadi.shp"
filename_faults = "./tectonic_maps_4dmb_2020_09_17/shape_files/faults_alcapadi.shp"
filename_gmt = "./tectonic_maps_4dmb_2020_09_17/GMT_example/alcapadi_polygons.gmt"
plot("alcapadi_polygons.gmt", show=true)
# 3. Select latitude and longitude range
lat_min = 37.0
lat_max = 49.0
lon_min = 4.0
lon_max = 20.0
region = [lon_min lon_max lat_min lat_max]

# 2. load desired topography data using GMT.jl and adapt the data to make it compatible with GeoData
G = gmtread(filename_topo, limits=[lon_min,lon_max,lat_min,lat_max]);
Lon,Lat,Depth    =   LonLatDepthGrid(G.x[1:end],G.y[1:end],0);
Depth[:,:,1]     =   1e-3*G.z';


# 3. load data for the geological map
TectUnitData = zeros(size(Lon));

# assign "tectonic units" as phases to surface grid using shapefile import
table = Shapefile.Table(filename_geo)
tectunits = unique(table.tect_unit); # get tectonic units

    for ishape = 1:length(table)
        phase = findall(==(0),cmp.(table.tect_unit[ishape], tectunits)); # assign a certain phase to the respective tectonic unit
        
        #xv = zeros(length(table.geometry[ishape].points),1); # preallocate the vector to hold the different points
        #yv = zeros(length(table.geometry[ishape].points),1); # preallocate the vector to hold the different points
        xv = Vector{Float64}(undef,length(table.geometry[ishape].points));
        yv = Vector{Float64}(undef,length(table.geometry[ishape].points));

        for ipoint = 1:length(table.geometry[ishape].points)
            xv[ipoint] = table.geometry[ishape].points[ipoint].x;
            yv[ipoint] = table.geometry[ishape].points[ipoint].y;
        end
    
        xa = vec(Lon);
        ya = vec(Lat);

        polygon = Point.(xv,yv);
        points  =  Point.(xa,ya);

        inside = [isinside(p, polygon; allowonedge=true) for p in points]
        
        inside = inside .& (vec(TectUnitData).==0);

        TectUnitData[inside] .= phase;
    end
    
    
# Set up the Data structure
Data_set        =   GeoData(Lon, Lat, Depth, (Topography=Depth*km,tect_unit=TectUnitData))

# Export the data structure to Paraview format
Write_Paraview(Data_set, "test_GeoMap")

















# Set up the Data structure
data_Topo        =   GeoData(Lon, Lat, Depth, (Topography=Depth*km,))

# Export the data structure to Paraview format
Write_Paraview(Data_set, "test_netcdf_ETOPO1")
