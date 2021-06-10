# Geophysical Model Generator

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://github.com/JuliaGeodynamics/GeophysicalModelGenerator.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://github.com/JuliaGeodynamics/GeophysicalModelGenerator.jl/dev)
[![Coverage](https://codecov.io/gh/JuliaGeodynamics/GeophysicalModelGenerator.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaGeodynamics/GeophysicalModelGenerator.jl)
[![Build Status](https://github.com/JuliaGeodynamics/GeophysicalModelGenerator.jl/workflows/CI/badge.svg)](https://github.com/JuliaGeodynamics/GeophysicalModelGenerator.jl/actions)


Creating consistent 3D images of various geophysical datasets is often challenging. The aim of this package is to help with this, by providing a number of routines to import and export data in various formats, for example to visualize them in [Paraview](https://www.paraview.org).

## Development roadmap
Currently, this package is under development. We plan to add the following features in approximately the following order:
-   Create VT* files from 3D volumes with lon/lat/depth and data (e.g., seismic tomography)
-   Create VT* files with horizontal or vertical cross-sections from lon/lat/depth and data
-   Create 1D lines with seismic reflection data below a station (amplitude)
-   Optionally interface with `GMT.jl` to allow importing  
  

## Installation 
You can install this and required dependencies within the julia package manager.
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
- [WriteVTK.jl](https://github.com/jipolanco/WriteVTK.jl) writes VTK files (to be openened with Paraview)

## Contributing
You are very welcome to request new features and point out bugs by opening an issue. You can also help by adding features and creating a pull request

## Funding
Develoment of this software package was funded by the German Research Foundation (DFG grant TH2076/7-1 KA3367/10-1), which is part of of the [SPP 2017 4DMB project](http://www.spp-mountainbuilding.de) project as well as by the European Research Council (ERC CoG #771143 - [MAGMA](https://magma.uni-mainz.de)).
