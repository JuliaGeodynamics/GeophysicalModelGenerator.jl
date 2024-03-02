# Extract topographic data from GMT.jl 

## Goal

In many cases, we want to add topographic data as well to our visualization. This tutorial shows how to use [GMT.jl](https://github.com/GenericMappingTools/GMT.jl) to download data from a certain region, and transfer that.

!!! note
    It used to be tricky to get `GMT.jl` installed and working correctly on your system but that has improved since version 1.0 which now comes with precompiled binaries. So as long as you make sure that your `GMT` version is >1.0, it should work.

## Steps

#### 1. Download topographic data of the Alpine region

The nice thing about GMT is that it automatically downloads data for you for a certain region and with a certain resolution. As this is a routine that you may use often in your daily workflow, we added the function `ImportTopo` that simplifies this. Note that this function only is available once `GMT` is loaded. 

```julia
julia> using GeophysicalModelGenerator, GMT
julia> Topo = ImportTopo([4,20,37,49], file="@earth_relief_01m.grd")
GeoData 
  size  : (960, 720, 1)
  lon   ϵ [ 4.0 : 19.983333333333334]
  lat   ϵ [ 37.0 : 48.983333333333334]
  depth ϵ [ -3.8725 km : 4.2495 km]
  fields: (:Topography,)
```
The data is available in different resolutions; see [here](http://gmt.soest.hawaii.edu/doc/latest/grdimage.html) for an overview. Generally, it is advisable to not use the largest resolution if you have a large area. 

#### 2. Save
Transforming this to Paraview is a piece of cake:

```julia
julia> Write_Paraview(Topo, "Topography_Alps") 
```
The result is shown here, together with Moho data

![Tutorial_GMT_topography](../assets/img/Tutorial_GMT_topography.png)

In case you are interested: we are employing the `oleron` scientific colormap [here](https://www.fabiocrameri.ch/colourmaps/).