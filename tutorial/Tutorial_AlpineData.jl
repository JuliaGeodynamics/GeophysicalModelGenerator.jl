# # Alpine Data Visualization
# This is a tutorial to:
# 1. download datasets from known sources
# 2. process and unify these datasets with GeophysicalModelGenerator
# 3. save the resulting dataset
# 4. export the datasets to Paraview
#
# This is a rather lengthy tutorial that combines different tutorials that you can find here (ADD LINK), but it will guide you through all the steps necessary to obtain a somewhat comprehensive view of the European Alps and their subsurface from a geodynamical point of view.

#  ## 1. Surface Topography
# In many cases, we want to add topographic data to our visualization. Here we use [GMT.jl](https://github.com/GenericMappingTools/GMT.jl) to download data from a certain region, and transfer that to GMG. 
# To add the GMT package, simply add it with the julia package manager:
# ```julia
# julia> ]
# (@v1.8) pkg> add GMT
# ```
# and load both GMG and GMT with:

using GeophysicalModelGenerator, GMT

# When loading both packages, several GMT routines within GMG will be loaded. One of these routines is the function ImportTopo, where one simply has to provide the region for which to download the topographic data and the data source.

Topo = ImportTopo([4,20,37,49], file="@earth_relief_01m.grd")
# The data is available in different resolutions; see [here](http://gmt.soest.hawaii.edu/doc/latest/grdimage.html) for an overview. Generally, it is advisable to not use the largest resolution if you have a large area. 

# If you have issues with loading the topography with GMT, there is also the alternative to download the data yourself and import it using Rasters.jl. This procedure is explained here. (ADD LINK TO TUTORIAL)

# We can now export this data to a VTK format so that we can visualize it with Paraview. To do so, GMG provides the function Write_Paraview:
Write_Paraview(Topo, "Topography_Alps") 

# Also, if you want to save this data for later use in julia, you can save it as JLD2 file using the function save_GMG:
save_GMG("Topography_Alps",Topo)


# ##  2. Moho topography
# When looking at data concerning the Alpine subsurface, we are often interested in the depth of the [Moho](https://en.wikipedia.org/wiki/Mohorovi%C4%8Di%C4%87_discontinuity). 
# ### 2.1 Download and import of the data
# Here, we will use the dataset from Mroczek et al. (2023). This dataset is publicly available and can be downloaded from [here](https://datapub.gfz-potsdam.de/download/10.5880.GFZ.2.4.2021.009NUEfb/2021-009_Mroczek-et-al_SWATHD_moho_jul22.csv). 
# To allow for downloadind such data, we use the julia Package Downloads.jl, which is a dependenvy of GMG and thus should be automatically available (is that so???). To download the Moho data in the current directory, simply type:
download_data("https://datapub.gfz-potsdam.de/download/10.5880.GFZ.2.4.2021.009NUEfb/2021-009_Mroczek-et-al_SWATHD_moho_jul22.csv","MohoMroczek2023.csv")

# Here the downloaded file gets the name MohoMroczek2023.dat and will be saved in the current directory. A quick look at the file shows that it contains a header consisting of 11 lines, one line with the description of the different columns and then the actual data. Have a look at the file yourself.
# To import the CSV file, we will use the package [DelimitedFiles.jl](https://github.com/JuliaData/DelimitedFiles.jl). Have a look at its documentation to see the different import options.
using DelimitedFiles
data_mroczek = readdlm("MohoMroczek2023.csv",',',header=false,skipstart=11)
# Note that we skipped the first 11 lines of the file as they contain the file header. The result of this operation will look like this:
#
# ```julia
# julia> data_mroczek = readdlm("MohoMroczek2023.csv",',',header=false,skipstart=11)
# 40786×10 Matrix{Any}:
#      "X"      "Y"      "Z"    "lat"    "lon"    "depth"   "tPs"   "k"   "interp"  "tag"
#  4185.26  1005.66  4640.51  47.152   13.5113  41.54      4.82    1.62  0          "PA"
#  4186.9   1005.95  4639.04  47.1319  13.5099  41.4893    4.81    1.62  0          "PA"
#  4188.54  1006.24  4637.57  47.1118  13.5085  41.4281    4.8     1.62  0          "PA"
#  4190.2   1006.53  4636.11  47.0917  13.5072  41.3568    4.8     1.63  0          "PA"
#  4191.86  1006.82  4634.66  47.0716  13.5058  41.2761    4.79    1.63  0          "PA"
#  4193.51  1007.12  4633.22  47.0516  13.5045  41.1865    4.78    1.63  0          "PA"
#  4195.18  1007.41  4631.78  47.0315  13.5031  41.0887    4.77    1.63  0          "PA"
#  4196.85  1007.71  4630.34  47.0114  13.5018  40.9835    4.76    1.63  0          "PA"
#  4198.53  1008.01  4628.91  46.9913  13.5004  40.8725    4.75    1.63  0          "PA"
#  4183.59  1007.45  4642.47  47.1721  13.5396  40.92      4.75    1.62  0          "PA"
#  4185.22  1007.73  4640.99  47.152   13.5382  40.8866    4.75    1.63  0          "PA"
#  4186.85  1008.02  4639.51  47.1319  13.5368  40.8435    4.75    1.63  0          "PA"
#  4188.49  1008.31  4638.04  47.1118  13.5355  40.7913    4.74    1.63  0          "PA"
#  4190.14  1008.6   4636.57  47.0917  13.5341  40.7305    4.73    1.63  0          "PA"
#     ⋮                                          ⋮                                  
#  4123.28  1111.52  4669.18  47.5537  15.0867  43.4301    5       1.67  0          "EU"
#  4124.02  1111.57  4666.69  47.5336  15.0848  44.7763    5.13    1.67  0          "EU"
#  4124.67  1111.59  4664.11  47.5136  15.0828  46.2535    5.28    1.67  0          "EU"
#  4125.23  1111.59  4661.42  47.4935  15.0808  47.868     5.44    1.67  0          "EU"
#  4125.71  1111.56  4658.63  47.4734  15.0788  49.6226    5.61    1.67  0          "EU"
#  4126.09  1111.51  4655.74  47.4533  15.0768  51.5103    5.79    1.66  0          "EU"
#  4126.41  1111.45  4652.77  47.4332  15.0749  53.496     6       1.66  0          "EU"
#  4126.7   1111.38  4649.78  47.4131  15.0729  55.5207    6.21    1.66  0          "EU"
#  4126.93  1111.28  4646.73  47.393   15.0709  57.6306    6.45    1.66  0          "EU"
#  4127.05  1111.16  4643.54  47.3729  15.0689  59.9233    6.7     1.66  0          "EU"
#  4127.0   1111.0   4640.2   47.3529  15.067   62.4358    6.98    1.66  1          "EU"
#  4126.2   1113.34  4649.82  47.4131  15.1     55.4728    6.21    1.66  0          "EU"
#  4126.41  1113.24  4646.74  47.393   15.098   57.6211    6.45    1.66  0          "EU"
#  4126.5   1113.11  4643.52  47.3729  15.096   59.9569    6.71    1.66  0          "EU"
# ```

# We are now only interested in the depth of the Moho at a given longitude/latitude. To obtain these values, we now have to extract columns 4-6. In addition, we also extract the 10th column, as it contains an identifier for the tectonic unit the respective point belongs to. 
lon        = zeros(size(data_mroczek,1)-1);lon        .= data_mroczek[2:end,5];
lat        = zeros(size(data_mroczek,1)-1);lat        .= data_mroczek[2:end,4];
depth      = zeros(size(data_mroczek,1)-1);depth      .= -1.0*data_mroczek[2:end,6]; # multiplied with -1, as we consider depth to be negative
tag        = string.(data_mroczek[2:end,10]); # get unit identifiers und convert them to strings

# As a next step, we then determine how many different tectonic units there are:
units = unique(tag) # get different units
# We will use these units later to save the Moho data separately for each tectonic unit.

# ### 2.2 Converting the data to a GMG dataset
# To convert this data to a GMG dataset, we now have to interpolate it to a regular grid. You can generate the respective grid with the GMG function  LonLatDepthGrid

Lon,Lat,Depth     =   LonLatDepthGrid(9.9:0.02:15.1,45.0:.02:49.0,0km);

# To interpolate the Moho data of the different units to this grid, we have here decided to employ a simple Nearest Neighbor interpolation for simplicity.
using NearestNeighbors
# !!! note Interpolating data is tricky and may result in unnecessary smoothing of the data. There are different ways to interpolate data on a regular grid. Have a look at our data interpolation tutorial to see the different possibilities.

# Now that we have generated the grid, we can loop over our different tectonic units, extract the relevant data points and interpolate them to the regular grid:
for iunit = 1:length(units)
    Dist              = zeros(size(Lon))
    # get all points belonging to the unit
    ind_unit    = findall( x -> occursin(units[iunit], x), tag) # index of the points belonging to that unit
    lon_tmp     = lon[ind_unit]
    lat_tmp     = lat[ind_unit]
    depth_tmp   = depth[ind_unit]
    
    # for later checking, we can save the original point data as a VTK file: 
    data_Moho = GeophysicalModelGenerator.GeoData(lon_tmp,lat_tmp,depth_tmp,(MohoDepth=depth_tmp*km,))
    filename = "Mroczek_Moho_" * units[iunit]
    Write_Paraview(data_Moho, filename, PointsData=true)
    
    # Now we create a KDTree for an effective nearest neighbor determination;
    kdtree = KDTree([lon_tmp';lat_tmp']; leafsize = 10)
    points = [Lon[:]';Lat[:]']
    idxs, dists = knn(kdtree, points, 1, true) # get the distance to the nearest data point
    dists = reduce(vcat,dists)
    idxs  = reduce(vcat,idxs)
    idxs  = reduce(vcat,idxs)
    
    # Having determined the nearest neighbor for each point in our regular grid, we can now directly assign the respective depth. Whenever the nearest neighbor is further than a certain distance away, we assume that there is no Moho at this point and do not assign a depth to that point.
    for i=1:length(idxs)
        if dists[i]<0.02 
            Depth[i] = depth_tmp[idxs[i]]*km
        else
            Depth[i] = NaN*km
        end
        Dist[i]  = dists[i]    
    end
    
    # As we will be using the data later, we would also like to provide some Metadata so that we know where it is coming from:
    Data_attribs   = Dict(
        "author"=>  "Mroczek et al.",
        "year"=> "2023",
        "doi"=>"https://doi.org/10.5880/GFZ.2.4.2021.009",
        "url"=>"https://nextcloud.gfz-potsdam.de/s/zB5dPNby6X2Kjnj",
    )
    
    # Finally, we can now export that data to VTK and save a jld2 file 
    Data_Moho = GeophysicalModelGenerator.GeoData(Lon, Lat, Depth, (MohoDepth=Depth,PointDist=Dist),Data_attribs)
    filename = "Mrozek_Moho_Grid_" * units[iunit]
    Write_Paraview(Data_Moho, filename)
    save_GMG(filename,Topo)

end

# ## 3. Seismicity
# Earthquakes are always interesting, so we will now import the seismicity data from ISC. 
# ### 3.1 Download and import
# ISC provides a method to download parts of it's catalogue via a web interface. See the description of the interface [here](http://www.isc.ac.uk/iscbulletin/search/webservices/bulletin/).
# We will now download all reviewed earthquake data between 1990 and 2000 in the same region as the extracted topography. We will only consider earthquakes with a magnitude larger than 3. The resulting dataset is quite large, so consider to either limit the time range or the magnitude range.
download_data("http://www.isc.ac.uk/cgi-bin/web-db-run?out_format=QuakeML&request=REVIEWED&searchshape=RECT&bot_lat=37&top_lat=49&left_lon=4&right_lon=20&start_year=1990&start_month=01&start_day=01&start_time=00:00:00&end_year=2000&end_month=12&end_day=31&end_time=15:00:00&min_mag=3.0&req_mag_type=Any&prime_only=on","ISCData.xml")

# ### 3.2 Read data into Julia
The main data-file, `ISC1.dat`, has 23 lines of comments (indicated with `#`), after which the data starts. We can use the julia package [https://github.com/JuliaData/CSV.jl](CSV.jl) to read in the data, and tell it that the data is separated by `,`.
```julia-repl
julia> using CSV, GeophysicalModelGenerator
julia> data_file        =   CSV.File("ISC1.dat",datarow=24,header=false,delim=',')
```
As this data contains a lot of information that we are not interested in at the moment and which is given in non-numeric formats (e.g. date, time etc.), we will use our helper function *ParseColumns_CSV_File* to only extract columns with numeric data.
```julia-repl
julia> data        =   ParseColumns_CSV_File(data_file, 14)
julia> lon         = data[:,2];
julia> lat         = data[:,1];
julia> depth       = -1* data[:,3];
julia> magnitude   = data[:,4];
```
Converting this data to a GeoStruct data and to export is to Paraview is then straightforward.
```julia-repl
julia> EQ_Data = GeoData(lon,lat,depth,(Magnitude=magnitude,Depth=depth));
julia> Write_Paraview(EQ_Data, "EQ_ISC", PointsData=true)
```
The result the looks like this (plotted here together with the topography):

![Tutorial_ISC](../assets/img/Tutorial_ISC.png)

In case you are interested: we are employing the `oleron` scientific colormap from [Fabio Crameri's scientific colormap package](https://www.fabiocrameri.ch/colourmaps/) here.

# ## 4. GPS data
# Besides data on the structure of the subsurface, it is also nice to see the dynamics of a region. Dynamic processes can be nicely seen in the surface velocities given by GPS data. 
# As GPS data consists of three-dimensional vectors, we have to treat it differently than the seismicity data in the previous section. The example is based on a paper by Sanchez et al. (2018) [https://essd.copernicus.org/articles/10/1503/2018/#section7](https://essd.copernicus.org/articles/10/1503/2018/#section7).
#
# ### 3.1. Download and import GPS data: 
# The data related to the paper can be downloaded from: [here](https://doi.pangaea.de/10.1594/PANGAEA.886889). There you will find links to several data sets. 
# Some are the data on the actual stations and some are interpolated data on a grid. Here, we will use the gridded data as an example (which interpolates the ), 
# and will therefore download the following data sets:
#
# - ALPS2017_DEF_HZ	Surface deformation model of the Alpine Region	[https://store.pangaea.de/Publications/Sanchez-etal_2018/ALPS2017_DEF_HZ.GRD](https://store.pangaea.de/Publications/Sanchez-etal_2018/ALPS2017_DEF_HZ.GRD)
# - ALPS2017_DEF_VT	Vertical deformation model of the Alpine Region	[https://store.pangaea.de/Publications/Sanchez-etal_2018/ALPS2017_DEF_VT.GRD](https://store.pangaea.de/Publications/Sanchez-etal_2018/ALPS2017_DEF_VT.GRD)

download_data("https://store.pangaea.de/Publications/Sanchez-etal_2018/ALPS2017_DEF_HZ.GRD","ALPS2017_DEF_HZ.GRD")
download_data("https://store.pangaea.de/Publications/Sanchez-etal_2018/ALPS2017_DEF_VT.GRD","ALPS2017_DEF_VT.GRD")

# Next, we have a look at the data. We will use the package `CSV.jl` to load the comma-separated data.
# Let's have a look at the file `ALPS2017_DEF_VT.GRD`. If we open it with a text editor, we see that the data starts at line 18, and has the following format:
#
# Column 1: Longitude [degrees]
# Column 2: Latitude [degrees]
# Column 3: Velocity in the height direction [m/a]
# Column 4: Uncertainty of the height component [m/a]
#
# So we have 4 columns with data values, and the data is separated by spaces.
# We can load that in julia as:
data_vhz = readdlm("ALPS2017_DEF_HZ.GRD",header=false,skipstart=18);
data_vz = readdlm("ALPS2017_DEF_VT.GRD",header=false,skipstart=17);
  
lon_vz =   data_vz[:,1]
lat_vz =   data_vz[:,2]
vz     =   data_vz[:,3];


#=

```

#### 2. Check & reshape vertical velocity

Let's have a look at the data, by plotting it:
```julia
julia> using Plots
julia> Plots.scatter(lon_Vz,lat_Vz)
```
![Tutorial_GPS_1](../assets/img/Tutorial_GPS_1.png)

So clearly, this is a fully regular grid.
We can determine the size of the grid with 
```julia
julia> unique(lon_Vz)
41-element Vector{Float64}:
  4.0
  4.3
  4.6
  4.9
  5.2
  5.5
  5.8
  ⋮
 14.5
 14.8
 15.1
 15.4
 15.7
 16.0
julia> unique(lat_Vz)
31-element Vector{Float64}:
 43.0
 43.2
 43.4
 43.6
 43.8
 44.0
 44.2
  ⋮
 48.0
 48.2
 48.4
 48.6
 48.8
 49.0
```
So we have a `41` by `31` grid. GMG requires 3D matrixes for the data (as we want to plot the results in paraview in 3D). That is why we first initialize 3D matrixes for `lon,lat,Vz`:
```julia
julia> lon, lat, Vz            =   zeros(41,31,1),zeros(41,31,1),zeros(41,31,1)
```
And we can reshape the vectors accordingly:
```julia
julia> lon[:,:,1]              =   reshape(lon_Vz,(41,31))
julia> lat[:,:,1]              =   reshape(lat_Vz,(41,31))
julia> Vz[:,:,1]               =   reshape(Vz_vec,(41,31))
```


#### 3. Load horizontal velocities
Next, we load the horizontal velocities from the file `ALPS2017_DEF_HZ.GRD`

```julia
julia> data_file                       =   CSV.File("ALPS2017_DEF_HZ.GRD",datarow=18,header=false,delim=' ');
julia> data                            =   ParseColumns_CSV_File(data_file, 6);
julia> lon_Hz, lat_Hz, Ve_Hz, Vn_Hz    =   data[:,1], data[:,2], data[:,3],  data[:,4];
```

Let's plot the data as well:
```julia
julia> Plots.scatter(lon_Hz,lat_Hz)
```
![Tutorial_GPS_2](../assets/img/Tutorial_GPS_2.png)

So it appears that the horizontal velocities are given on the same regular grid as well, but not in the water. 
This thus requires a bit more work. The strategy we take is to first define 2D matrixes with horizontal velocities with the same size as Vz which are initialized with `NaN` (not a number), which is treated specially by Paraview.

```julia
julia> Ve = ones(size(Vz))*NaN;
julia> Vn = ones(size(Vz))*NaN;
```

Next we loop over all points in `lon_Hz,lat_Hz` and place them into the 2D matrixes:
```julia
julia>  for i in eachindex(lon_Hz)
            ind = intersect(findall(x->x==lon_Hz[i], lon), findall(x->x==lat_Hz[i], lat))
            Ve[ind] .= Ve_Hz[i];
            Vn[ind] .= Vn_Hz[i];
        end
```

At this stage, we have horizontal and vertical velocities in units of `m/yr`. Yet, given the small velocities in the Alps, it makes more sense to have them in units of `mm/yr`:
```julia
julia> Vz = Vz*1000;
julia> Ve = Ve*1000;
julia> Vn = Vn*1000;
```
And the magnitude is:
```julia
julia> Vmagnitude  =   sqrt.(Ve.^2 + Vn.^2 + Vz.^2);  
```

#### 4. Interpolate topography on grid
At this stage we have the 3D velocity components on a grid. Yet, we don't have information yet about the elevation of the stations (as the provided data set did not give this). 
We could ignore that and set the elevation to zero, which would allow saving the data directly to paraview.

Yet, a better way is to load the topographic map of the area and interpolate the elevation to the velocity grid. We are using the `GMT.jl` to load the topographic data:
```julia
julia> using GMT
julia> Elevation =   gmtread("@earth_relief_01m.grd", limits=[3,17,42,50]);
```

Next, we use the `Interpolations.jl` package to interpolate the topography:
```julia
julia> using Interpolations
julia> interpol = LinearInterpolation((Elevation.x[1:end-1], Elevation.y[1:end-1]), Elevation.z');    
julia> height   = interpol.(lon,lat)/1e3;
```

#### 5. Saving and plotting in Paraview
At this stage, we have all we need. As the velocity is a vector field, we need to save it as a data structure with 3 components. When saving to paraview, GMG internally does a vector transformation. As this transformation does not retain the `east/north` components of the velocity field, it is a good idea to save them as separate fields so we can color the vectors accordingly in Paraview. Also note that we do not attach units to the vector fields, but we do have them for the scalar fields:

```julia
julia> GPS_Sanchez_grid = GeoData(lon,lat,height,(Velocity_mm_year=(Ve,Vn,Vz),V_north=Vn*mm/yr, V_east=Ve*mm/yr, V_vertical=Vz*mm/yr, Vmagnitude = Vmagnitude*mm/yr, Topography = height*km))
GeoData 
  size  : (41, 31, 1)
  lon   ϵ [ 4.0 : 16.0]
  lat   ϵ [ 43.0 : 49.0]
  depth ϵ [ -2.6545 km : 3.426 km]
  fields: (:Velocity_mm_year, :V_north, :V_east, :V_vertical, :Vmagnitude, :Topography)
```
Saving this to paraview is as always:
```julia
julia> Write_Paraview(GPS_Sanchez_grid, "GPSAlps_Sanchez_2017_grid")
```

Opening and plotting the vertical field gives:
![Tutorial_GPS_3](../assets/img/Tutorial_GPS_3.png)

In order to plot the velocities as arrows, you need to select the `Glyph` tool (red circle). Also specify `Velocity_mm_year ()` as both Orientation and Scale Array, and add `50` as scale factor. Once you push `Apply` it should look like:
![Tutorial_GPS_4](../assets/img/Tutorial_GPS_4.png)

The arrows can now be colored by the individual velocity components or its magnitude.


=#
