# Geophysical Model Generator

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliageodynamics.github.io/GeophysicalModelGenerator.jl/dev)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliageodynamics.github.io/GeophysicalModelGenerator.jl/dev/)
[![Build Status](https://github.com/JuliaGeodynamics/GeophysicalModelGenerator.jl/workflows/CI/badge.svg)](https://github.com/JuliaGeodynamics/GeophysicalModelGenerator.jl/actions)

Creating consistent 3D images of geophysical and geological datasets and turning that into an input model for geodynamic simulations is often challenging. The aim of this package is to help with this, by providing a number of routines to easily import data and create a consistent 3D visualisation from it in the VTK-toolkit format, which can for example be viewed with [Paraview](https://www.paraview.org). In addition, we provide a range of tools that helps to generate input models to perform geodynamic simulations and import the results of such simulations back into julia. 

![README_img](./docs/src/assets/img/Readme_pic.png)
### Contents
* [Main features](#main-features) 
* [Usage](#usage)
* [Installation](#installation)
* [Dependencies](#dependencies)
* [Alpine Data](#visualising-alpine-data)
* [Contributing](#contributing)
* [Development roadmap](#development-roadmap)
* [Funding](#funding)

## Main features
Some of the key features are:
- Create 3D volumes of seismic tomography models.
- Handle 2D data (e.g., along a cross-section), including surfaces such as the Moho depth.
- Plot data along lines (e.g., drillholes) or at points (e.g., earthquake locations, GPS velocities).
- Handle both scalar and vector data sets.
- Grab screenshots of cross-sections or maps in published papers and view them in 3D (together with other data).
- Create a consistent overview that includes all available data of a certain region.
- Create initial model setups for the 3D geodynamic code [LaMEM](https://bitbucket.org/bkaus/lamem/src/master/).
- Import LaMEM timesteps. 

All data is transformed into either a `GeoData` or a `UTMData`  structure which contains info about `longitude/latitude/depth`, `ew/ns/depth` coordinates along with an arbitrary number of scalar/vector datasets, respectively. All data can be exported to Paraview with the `Write_Paraview` routine, which transfers the data to a `ParaviewData` structure (that contains Cartesian Earth-Centered-Earth-Fixed (ECEF) `x/y/z` coordinates, used for plotting)
  
## Usage 
The best way to learn how to use this is to install the package (see below) and look at the tutorials in the [manual](https://juliageodynamics.github.io/GeophysicalModelGenerator.jl/dev/).

## Installation 
First, you need to install julia on your machine. We recommend to use the binaries from [https://julialang.org](https://julialang.org).
Next, start julia and switch to the julia package manager using `]`, after which you can add the package.
```julia
julia> ]
(@v1.6) pkg> add GeophysicalModelGenerator
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
We rely on a number of additional packages, which are all automatically installed.
- [GeoParams.jl](https://github.com/JuliaGeodynamics/GeoParams.jl) Defines dimensional units, and makes it easy to convert for km/s to m/s, etc.
- [WriteVTK.jl](https://github.com/jipolanco/WriteVTK.jl) writes VTK files (to be opened with Paraview).
- [ImageIO.jl](https://github.com/JuliaIO/ImageIO.jl), [FileIO.jl](https://github.com/JuliaIO/FileIO.jl), [Colors.jl](https://github.com/JuliaGraphics/Colors.jl) to import screenshots from papers.
- [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl) for interpolations (for example related to importing screenshots).


## Visualising Alpine data
We have used this package to interpret various data sets of the Alps (mostly openly available, sometimes derived from published papers). You can download the resulting paraview files here (using the `*.vts` format), where we also included the julia scripts to do the work (some of which are also described in more detail in the tutorials). Just unzip the files and open the corresponding `*.vts` in Paraview. 

[https://seafile.rlp.net/d/22b0fb85550240758552/](https://seafile.rlp.net/d/22b0fb85550240758552/)

If you want your data be included here as well, give us an email (or even better: send the files with julia scripts).
## Contributing
You are very welcome to request new features and point out bugs by opening an issue. You can also help by adding features and creating a pull request.

## Development roadmap
In the pipeline: 
- More tutorials
- Add more import tools.
- Compute gravity anomalies for lon/lat datasets, rather than just x/y/z.
- Provide an interface to [geomIO](https://bitbucket.org/geomio/geomio/wiki/Home) (currently being translated from MATLAB to python) in order to allow creating 3D geometric model setups by drawing in Inkscape. 
- Provide tools to create and export 3D geodynamic model setups. 
 
## Funding
Development of this software package was funded by the German Research Foundation (DFG grants TH2076/7-1 and KA3367/10-1), which are part of the [SPP 2017 4DMB project](http://www.spp-mountainbuilding.de) project as well as by the European Research Council under grant ERC CoG #771143 - [MAGMA](https://magma.uni-mainz.de).
