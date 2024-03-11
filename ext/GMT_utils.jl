# NOTE: these are useful routines that are only made available when the GMT package is already loaded in the REPL
module GMT_utils

import GeophysicalModelGenerator: ImportTopo, ImportGeoTIFF

# We do not check `isdefined(Base, :get_extension)` as recommended since
# Julia v1.9.0 does not load package extensions when their dependency is
# loaded from the main environment.
if VERSION >= v"1.9.1"
  using GMT
else
  using ..GMT
end

using GeophysicalModelGenerator: LonLatDepthGrid, GeoData, UTMData, km, remove_NaN_Surface!

println("Loading GMT routines within GMG")


"""
    Topo = ImportTopo(limits; file::String="@earth_relief_01m.grd", maxattempts=5) 
    
Uses `GMT` to download the topography of a certain region, specified with limits=[lon_min, lon_max, lat_min, lat_max].
Sometimes download fails because of the internet connection. We do `maxattempts` to download it.

Note: 
====
- latitude values in the southern hemisphere should have a minus sign (e.g., -2.8)
- longitude values that are "west" should *either* come with a minus sign *or* are defined by values >180

| Dataset                 |   Resolution |   Description                                               |
|:----------------        | ------------ | ----------------------------------------------------------- |
| "@earth\\_relief\\_01s" |	1 arc sec 	 | SRTM tiles (14297 tiles, land only, 60S-60N) [NASA/USGS]    |
| "@earth\\_relief\\_03s"	|   3 arc sec	 | SRTM tiles (14297 tiles, land only, 60S-60N) [NASA/USGS]    |
| "@earth\\_relief\\_15s"	|  15 arc sec	 | SRTM15+ [David Sandwell, SIO/UCSD]                          |
| "@earth\\_relief\\_30s"	|  30 arc sec	 | SRTM30+ [Becker et al., 2009, SIO/UCSD]                     |
| "@earth\\_relief\\_01m"	|   1 arc min	 | ETOPO1 Ice surface [NEIC/NOAA]                              |
| "@earth\\_relief\\_02m"	|   2 arc min	 | ETOPO2v2 Ice surface [NEIC/NOAA]                            |
| "@earth\\_relief\\_03m"	|   3 arc min	 | ETOPO1 after Gaussian spherical filtering (5.6 km fullwidth)|
| "@earth\\_relief\\_04m"	|   4 arc min	 | ETOPO1 after Gaussian spherical filtering (7.5 km fullwidth)|
| "@earth\\_relief\\_05m"	|   5 arc min	 | ETOPO1 after Gaussian spherical filtering (9 km fullwidth)  |
| "@earth\\_relief\\_06m"	|   6 arc min	 | ETOPO1 after Gaussia30n spherical filtering (10 km fullwidth) |
| "@earth\\_relief\\_10m"	|  10 arc min	 | ETOPO1 after Gaussian spherical filtering (18 km fullwidth) |
| "@earth\\_relief\\_15m"	|  20 arc min	 | ETOPO1 after Gaussian spherical filtering (28 km fullwidth) |
| "@earth\\_relief\\_20m"	|  20 arc min	 | ETOPO1 after Gaussian spherical filtering (37 km fullwidth) |
| "@earth\\_relief\\_30m"	|  30 arc min	 | ETOPO1 after Gaussian spherical filtering (55 km fullwidth) |
| "@earth\\_relief\\_60m"	|  60 arc min	 | ETOPO1 after Gaussian spherical filtering (111 km fullwidth)|

*Note*: this routine is only available once the GMT.jl package is loaded in the REPL

# Example 
```julia-repl
julia> Topo = ImportTopo([4,20,37,49]);
GeoData 
  size  : (960, 720, 1)
  lon   ϵ [ 4.0 : 19.983333333333334]
  lat   ϵ [ 37.0 : 48.983333333333334]
  depth ϵ [ -3.8725 km : 4.2495 km]
  fields: (:Topography,)
```
And you can save this to Paraview with
```julia
julia> Write_Paraview(Topo,"Topo_Alps")
1-element Vector{String}:
 "Topo_Alps.vts"
```
"""
function ImportTopo(limits; file::String="@earth_relief_01m.grd", maxattempts=5)

    # Correct if negative values are given (longitude coordinates that are west)
    ind = findall(limits[1:2] .< 0); 
    
    if (limits[1] < 0) && (limits[2] < 0)
      limits[ind] .= 360 .+ limits[ind]; 
      limits[1:2] = sort(limits[1:2])   
    end

    # Download topo file  - add a few attempts to do so
    G = [];
    attempt = 0
    while attempt<maxattempts
      try
        G       =   gmtread(file, limits=limits, grid=true);
        break
      catch
        @warn "Failed downloading GMT topography on attempt $attempt/$maxattempts"
        sleep(5)  # wait a few sec
      end
      attempt += 1
    end
    if isempty(G)
      error("Could not download GMT topography data")
    end

    # Transfer to GeoData
    nx,ny           =   size(G.z,2), size(G.z,1)
    Lon,Lat,Depth   =   LonLatDepthGrid(G.x[1:nx],G.y[1:ny],0);
    Depth[:,:,1]    =   1e-3*G.z';
    Topo            =   GeoData(Lon, Lat, Depth, (Topography=Depth*km,))
    
    return Topo

end

"""
  ImportTopo(; lat::Vector{2}, lon::Vector{2}, file::String="@earth_relief_01m.grd", maxattempts=5)

Imports topography (using GMT), by specifying keywords for latitude and longitude ranges

# Example
=========
```julia
julia> Topo = ImportTopo(lat=[30,40], lon=[30, 50] )
```
The values can also be defined as tuples:
```julia
julia> Topo = ImportTopo(lon=(-50, -40), lat=(-10,-5), file="@earth_relief_30s.grd")
```

"""
ImportTopo(; lat=[37,49], lon=[4,20], file::String="@earth_relief_01m.grd", maxattempts=5) = ImportTopo([lon[1],lon[2], lat[1], lat[2]], file=file, maxattempts=maxattempts)


"""
  data_GMT = ImportGeoTIFF(fname::String; fieldname=:layer1, negative=false, iskm=true, NorthernHemisphere=true, constantDepth=false, removeNaN_z=false, removeNaN_field=false)

This imports a GeoTIFF dataset (usually containing a surface of some sort) using GMT.
The file should either have `UTM` coordinates of `longlat` coordinates. If it doesn't, you can 
use QGIS to convert it to `longlat` coordinates.

Optional keywords:
- `fieldname` : name of the field (default=:layer1)
- `negative`  : if true, the depth is multiplied by -1 (default=false)
- `iskm`      : if true, the depth is multiplied by 1e-3 (default=true)
- `NorthernHemisphere`: if true, the UTM zone is set to be in the northern hemisphere (default=true); only relevant if the data uses UTM projection
- `constantDepth`: if true we will not warp the surface by z-values, but use a constant value instead
- `removeNaN_z`  : if true, we will remove NaN values from the z-dataset
"""
function ImportGeoTIFF(fname::String; fieldname=:layer1, negative=false, iskm=true, NorthernHemisphere=true, constantDepth=false, removeNaN_z=false, removeNaN_field=false)
  G = gmtread(fname);

  # Transfer to GeoData
  nx,ny = length(G.x)-1, length(G.y)-1
  Lon,Lat,Depth   =   LonLatDepthGrid(G.x[1:nx],G.y[1:ny],0);
  if  hasfield(typeof(G),:z) 
    Depth[:,:,1]    =   G.z';
    if negative
      Depth[:,:,1]  =   -G.z';
    end
    if iskm
      Depth    *=   1e-3*km;
    end
  end

  # Create GeoData structure
  data = zero(Lon)
  if hasfield(typeof(G),:z)
    data = Depth
  
  elseif hasfield(typeof(G),:image)
    if length(size(G.image)) == 3
      data = permutedims(G.image,[2, 1, 3]);
    elseif length(size(G.image)) == 2
      data[:,:,1] = G.image'
    end

  end
  
  if removeNaN_z
    remove_NaN_Surface!(Depth, Lon, Lat)
  end
  if removeNaN_field
    remove_NaN_Surface!(data, Lon, Lat)
  end
  data_field  = NamedTuple{(fieldname,)}((data,));

  if constantDepth
    Depth = zero(Lon)
  end

  if contains(G.proj4,"utm")
    zone = parse(Int64,split.(split(G.proj4,"zone=")[2]," ")[1]); # retrieve UTM zone
    data_GMT    = UTMData(Lon, Lat, Depth, zone, NorthernHemisphere, data_field)
  
  elseif contains(G.proj4,"longlat") 
    data_GMT    = GeoData(Lon, Lat, Depth, data_field)

  else
    error("I'm sorry, I don't know how to handle this projection yet: $(G.proj4)\n
           We recommend that you transfer your GeoTIFF to longlat by using QGIS \n
           Open the GeoTIFF there and Export -> Save As , while selecting \"EPSG:4326 - WGS 84\" projection.")
  end

  return data_GMT
end


end
