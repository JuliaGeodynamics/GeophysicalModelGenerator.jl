# Extract topographic data from GMT.jl

## Goal

In many cases, we want to add topographic data as well to our visualization. This tutorial shows how to use [GMT.jl](https://github.com/GenericMappingTools/GMT.jl) to download data from a certain region, and transfer that.

!!! note
    It may be tricky to get GMT.jl installed and working correctly on your system (at least until someone provides a BinaryBuilder package for julia, that is). You first need to have a working version of GMT on your system and only after that, you can install `GMT.jl`. See the installation instructions on their webpage for details.
    On a MacBook Pro, a tested procedure to install GMT and to make it work with julia is to directly install the binaries for Julia, GMT (and possibly Ghostscript) and not use any package manager (such as spack or homebrew).


## Steps

#### 1. Download topographic data of the Alpine region

The nice thing about GMT is that it automatically downloads data for you, from a certain region:

```julia
julia> using GMT
julia> G = gmtread("@earth_relief_01m.grd", limits=[4,20,37,49], grid=true);
```
The data is available in different resolutions; see [here](http://gmt.soest.hawaii.edu/doc/latest/grdimage.html) for an overview. Generally, it is advisable to not use the largest

#### 2. Save

Transforming this to Paraview is piece of cake:

```julia
julia> Lon,Lat,Depth    =   LonLatDepthGrid(G.x[1:end-1],G.y[1:end-1],0);
julia> Depth[:,:,1]     =   1e-3*G.z';
julia> data_Topo        =   GeoData(Lon, Lat, Depth, (Topography=Depth*km,))
julia> Write_Paraview(data_Topo, "Topography_Alps")
```
The result is shown here, together with Moho data

![Tutorial_GMT_topography](../assets/img/Tutorial_GMT_topography.png)

In case you are interested: we are employing the `oleron` scientific colormap here.
