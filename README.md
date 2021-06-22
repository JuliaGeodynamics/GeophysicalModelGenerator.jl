# Geophysical Model Generator

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliageodynamics.github.io/GeophysicalModelGenerator.jl/dev)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliageodynamics.github.io/GeophysicalModelGenerator.jl/dev/)
[![Build Status](https://github.com/JuliaGeodynamics/GeophysicalModelGenerator.jl/workflows/CI/badge.svg)](https://github.com/JuliaGeodynamics/GeophysicalModelGenerator.jl/actions)

Creating consistent 3D images of geophysical and geological datasets is often challenging. The aim of this package is to help with this, by providing a number of routines to easily import data and create a consistent 3D visualisation from it in the VTK-toolkit format, which can for example be viewed with [Paraview](https://www.paraview.org).

## Main features
Some of the key features are:
- Create 3D volumes of seismic tomography models.
- Handle 2D data (e.g., along a cross-section), including surfaces such as the Moho depth.
- Plot data along lines (e.g., drillholes) or at points (e.g., earthquake locations, GPS velocities).
- Handle both scalar and vector data sets.
- Grab screenshots of cross-sections or maps in published papers and view them in 3D (together with other data).
- Create a consistent overview that includes all available data of a certain region.

All data is transformed ibto a `GeoData` structure which contains info about the longitude, latitude and depth  along with an arbitrary number of scalar/vector datasets.
  
## Usage 
The best way to learn how to use this is to install the package (see below) and look at the tutorials in the [manual](https://juliageodynamics.github.io/GeophysicalModelGenerator.jl/dev/).


## Development roadmap
In the pipeline: 
- More tutorials
- Add more import tools.
- Compute gravity anomalies for lon/lat datasets, rather than just x/y/z.
- Provide an interface to [geomIO](https://bitbucket.org/geomio/geomio/wiki/Home) (currently being translated from MATLAB to python) in order to allow creating 3D geometric model setups by drawing in Inkscape. 
- Provide tools to create and export 3D geodynamic model setups. 
 
  
## Installation 
First, you need to install julia on your machine. We recommend to use the binaries from [https://julialang.org](https://julialang.org).
Next, start julia and install this and required dependencies within the julia package manager.
Note that we first install the [GeoParams.jl](https://github.com/JuliaGeodynamics/GeoParams.jl) package after which we install GeophysicalModelGenerator.
```julia
julia> ]
(@v1.6) pkg> add https://github.com/JuliaGeodynamics/GeoParams.jl
(@v1.6) pkg> add https://github.com/JuliaGeodynamics/GeophysicalModelGenerator
```
You can test whether it works on your system with
```julia
julia> ]
(@v1.6) pkg> test GeophysicalModelGenerator
```
and use it with
```julia
julia> using GeophysicalModelGenerator
```

## Dependencies
We rely on a number of additional packages. All of them are automatically installed, except `GeoParams.jl`, which you currenty have to add yourself
- [GeoParams.jl](https://github.com/JuliaGeodynamics/GeoParams.jl) Defines dimensional units, and makes it easy to convert for km/s to m/s, etc.
- [WriteVTK.jl](https://github.com/jipolanco/WriteVTK.jl) writes VTK files (to be opened with Paraview).
- [ImageIO.jl](https://github.com/JuliaIO/ImageIO.jl), [FileIO.jl](https://github.com/JuliaIO/FileIO.jl), [Colors.jl](https://github.com/JuliaGraphics/Colors.jl) to import screenshots from papers.
- [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl) for interpolations (for example related to importing screenshots).


## Visualising alpine data
We have used this package to interpret various data sets of the Alps (mostly openly available, sometimes derived from published papers). You can download the resulting paraview files here (using the `*.vts` format), where we also included the julia scripts to do the work (some of which are also described in more detail in the tutorials). Just unzip the files and open the corresponding `*.vts` in Paraview. 

[https://seafile.rlp.net/d/22b0fb85550240758552/](https://seafile.rlp.net/d/22b0fb85550240758552/)

If you want your data be included here as well, give us an email (or even better: send the files with julia scripts).
## Contributing
You are very welcome to request new features and point out bugs by opening an issue. You can also help by adding features and creating a pull request.

## Funding
Development of this software package was funded by the German Research Foundation (DFG grants TH2076/7-1 and KA3367/10-1), which are part of the [SPP 2017 4DMB project](http://www.spp-mountainbuilding.de) project as well as by the European Research Council under grant ERC CoG #771143 - [MAGMA](https://magma.uni-mainz.de).
