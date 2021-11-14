# NOTE: these are useful routines that are only loaded when the GMT package is already loaded in the REPL

using .GMT

export ImportTopo

"""
    Topo = ImportTopo(limits; file::String="@earth_relief_01m.grd") 
    
Uses GMT to download the topography of a certain region, specified with limits=[lon_min, lon_max, lat_min, lat_max] 

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

*Note*: this routine is only available once the GMT.jl package is loaded on the REPL

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
function ImportTopo(limits; file::String="@earth_relief_01m.grd")

    # Correct if negative values are given (longitude coordinates that are west)
    ind = findall(limits[1:2] .< 0); 
    
    if (limits[1] < 0) && (limits[2] < 0)
      limits[ind] .= 360 .+ limits[ind]; 
      limits[1:2] = sort(limits[1:2])   
    end

    # Download topo file  
    G               =   gmtread(file, limits=limits, grid=true);

    # Transfer to GeoData
    nx,ny           =   size(G.z,2), size(G.z,1)
    Lon,Lat,Depth   =   LonLatDepthGrid(G.x[1:nx],G.y[1:ny],0);
    Depth[:,:,1]    =   1e-3*G.z';
    Topo            =   GeoData(Lon, Lat, Depth, (Topography=Depth*km,))
    
    return Topo

end

"""
  ImportTopo(; lat::Vector{2}, lon::Vector{2}, file::String="@earth_relief_01m.grd")

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
ImportTopo(; lat=[37,49], lon=[4,20], file::String="@earth_relief_01m.grd") = ImportTopo([lon[1],lon[2], lat[1], lat[2]], file=file)
