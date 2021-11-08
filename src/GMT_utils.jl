# NOTE: these are useful routines that are only loaded when the GMT package is already loaded in the REPL

using .GMT

export ImportTopo

"""
    Topo = ImportTopo(limits; file::String="@earth_relief_01m.grd") 
    
Uses GMT to download the topography of a certain region, specified with limits=[lon_min, lon_max, lat_min, lat_max] 

| Dataset           |   Resolution   |   Description                                               |
|:----------------  | -------------- | ----------------------------------------------------------- |
| @earth_relief_01s |	1 arc sec 	 | SRTM tiles (14297 tiles, land only, 60S-60N) [NASA/USGS]    |
| @earth_relief_03s	|   3 arc sec	 | SRTM tiles (14297 tiles, land only, 60S-60N) [NASA/USGS]    |
| @earth_relief_15s	|  15 arc sec	 | SRTM15+ [David Sandwell, SIO/UCSD]                          |
| @earth_relief_30s	|  30 arc sec	 | SRTM30+ [Becker et al., 2009, SIO/UCSD]                     |
| @earth_relief_01m	|   1 arc min	 | ETOPO1 Ice surface [NEIC/NOAA]                              |
| @earth_relief_02m	|   2 arc min	 | ETOPO2v2 Ice surface [NEIC/NOAA]                            |
| @earth_relief_03m	|   3 arc min	 | ETOPO1 after Gaussian spherical filtering (5.6 km fullwidth)|
| @earth_relief_04m	|   4 arc min	 | ETOPO1 after Gaussian spherical filtering (7.5 km fullwidth)|
| @earth_relief_05m	|   5 arc min	 | ETOPO1 after Gaussian spherical filtering (9 km fullwidth)  |
| @earth_relief_06m	|   6 arc min	 | ETOPO1 after Gaussian spherical filtering (10 km fullwidth) |
| @earth_relief_10m	|  10 arc min	 | ETOPO1 after Gaussian spherical filtering (18 km fullwidth) |
| @earth_relief_15m	|  20 arc min	 | ETOPO1 after Gaussian spherical filtering (28 km fullwidth) |
| @earth_relief_20m	|  20 arc min	 | ETOPO1 after Gaussian spherical filtering (37 km fullwidth) |
| @earth_relief_30m	|  30 arc min	 | ETOPO1 after Gaussian spherical filtering (55 km fullwidth) |
| @earth_relief_60m	|  60 arc min	 | ETOPO1 after Gaussian spherical filtering (111 km fullwidth)|

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

    G               =   gmtread(file, limits=limits, grid=true);
    Lon,Lat,Depth   =   LonLatDepthGrid(G.x[1:end-1],G.y[1:end-1],0);
    Depth[:,:,1]    =   1e-3*G.z';
    Topo            =   GeoData(Lon, Lat, Depth, (Topography=Depth*km,))
    
    return Topo

end