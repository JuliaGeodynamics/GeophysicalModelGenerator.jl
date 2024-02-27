---
title: 'GeophysicalModelGenerator.jl: A Julia package to visualize geoscientific data and generate model '
tags:
  - Julia
  - geosciences
  - geodynamics
  - tectonics
  - geophysics
  - computational geosciences
authors:
  - name: Boris J.P. Kaus
    orcid: 0000-0002-0247-8660
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Marcel Thielmann
    orcid: 0000-0003-1185-3730
    affiliation: 2
  - name: Arne Spang
    orcid: 0000-0002-6107-0403
    affiliation: 2
  - name: Luca de Siena
    orcid: 0000-0002-3615-5923
    affiliation: 3
  - name: Pascal Aellig
    orcid: 0009-0008-9039-5646
    affiliation: 1
  - name: Jacob Frasukiewicz
    #orcid: 0000-0003-1185-3730
    affiliation: 1
affiliations:
 - name: Johannes Gutenberg University Mainz, Germany
   index: 1
 - name: University of Bayreuth, Germany
   index: 2
 - name: Bologna University, Italy
   index: 3

date: 27 February 2024
bibliography: paper.bib
---

# Summary

Geoscientific data exists in wide different variety of formats. Yet, to make a consistent interpretation of a certain region, it is often helpful to jointly visualize all available data using the same coordinates, and compare, for example seismic tomography, surface geology, Moho depth, earthquakes with GPS velocities. If one wishes to create mechanical or thermomechanical geodynamic simulations of the region, creating an input model that honors the given constraints is helpful. And since most numerical codes work in cartesian boxes, we need tools to project the data from geographic to cartesian coordinates.

A significant challenge in doing this is that there is no standard data format for geoscientific data sets. Seismic tomography, for example, may come as ASCII data with lon/lat/depth, or as NetCDF files, with the ordering of the data typically differing from one dataset to the other. In ideal cases, geological surfaces may be provided as GeoTIFFs which can be readily imported. Yet, in many cases, the underlying data are not even available in digital format and we only have access to figures in the publications. Yet, it is still helpful to image these in 3D, along with more recent, digitally available, datasets.

The aim of the GeophysicalModelGenerate package is twofold:

1) Simplify collecting and jointly visualising a wide variety of geoscientific data that is provided as point (e.g. earthquake locations), surface (e.g., topography) or volumetric data (seismic tomography)  
2) Generate input setups for numerical models (and read back the output simulations)

# Statement of need

`GeophysicalModelGenerator.jl` is a Julia package that helps collecting and visualizing a wide variety of geophysical and geoscientific data in a coherent manner. It also simplifies the process of generating a 2D or 3D models that can, for example, be used as input models in geodynamic simulations. It provides functions that transfer data from one format to the other, or project them from geographic (Longitude/Latitude/Depth) or (UTM) coordinates to cartesian coordinates (in kilometers). It allows performing tasks such as creating cross-sections though volumetric data, import screenshots from published papers, download digital elevation data and save the resulting data in VTK format, that can be visualized with [Paraview](www.paraview.org).

Many geoscientists likely have their own python/matlab/bash scripts to visualize their own data and thus perform part of this job already. Yet, we are not aware of an open-source package that performs all these tasks in one place while also making a straightforward link.


# Related software packages
Perhaps the most widely used package in geophysics to create figures or maps is the Generic Mapping Tools ([`GMT`]()), which also provides a Julia interface [GMT.jl](https://github.com/GenericMappingTools/GMT.jl). It mostly focusses on generating maps and postscript/pdf images and is therefore not ideally suited for interactive 3D data visualisation or to generate input models for geodynamic codes.
 
The [`Geodynamic World Builder`](https://github.com/GeodynamicWorldBuilder/WorldBuilder) is a C++ library to create model setups [`@se-10-1785-2019`]. The focus is on generating input models for geodynamic simulations, such as an initial subduction zone and related thermal structure. It has C and Fortran wrappers and can thus be embedded in geodynamic codes. Users of the `Geodynamic World Builder` have to generate JSON files to define the model geometry, which is less interactive than by using the Julia `REPL`. In the currently available version, there is no straighforward way to integrate existing geophysical/geological data in the workflow and compare model results with them.   

[GemPy]() is a is Python-based, open-source geomodeling library that can construct 3D geological models of folded structures, fault networks and unconformities, while taking uncertainties into account. Its focus is on creating geometric models rather than on integrating a wide variety of geoscientific datasets.

There are also a number of commercial software solutions as well: 
- [Petrel subsurface software](https://www.software.slb.com/products/petrel) (by Schlumberger), which is mostly used by the hydrocarbon industry and is particularly powerful in integrating seismic reflection and well-data, 
- [GOCAD Mining Suite](https://www.mirageoscience.com/mining-industry-software/gocad-mining-suite/) (by MiraGeoscience) helps generate geometric models of the sub surface in the vicinity of mines, based on sparse geological measurements and drillhole data.
- [GeoModeller](https://www.intrepid-geophysics.com/products/geomodeller/) (by Intrepid Geophysics) creates surface-near geometric geological models by implicit modelling of surface measurements while taking geophysical constraints into account.    

In all cases, the commercial license fees are far beyond what most researchers can afford, even if reduced license fees are often available for academia. The closed-source nature of the software packages makes them non-extendable by the community.

The `GeophysicalModelGenerator.jl` package is already used to generate input models for the geodynamic codes [LaMEM](https://github.com/UniMainzGeo/LaMEM), [JustRelax.jl](https://github.com/PTsolvers/JustRelax.jl), and [MagmaThermokinematics.jl](https://github.com/boriskaus/MagmaThermoKinematics.jl). It is also already used in a number of shortcourses and lectures at the university of Mainz.

# Basic usage

The core of the package consists of the  `GeoData`, `UTMData`, `ParaviewData` and `CartData` structures which essentially hold the 3D data. We can illustrate this 

```julia
julia> using GeophysicalModelGenerator
julia> Write_Paraview()
```


# Examples of usage
A full set of tutorials is given on the [GitHub]() page. In the following we present a number of examples that illustrate 
In the following

### Visualize data of the Alps

### Jura mountains

### La Palma volcanic eruption


### Slab model setup 
To illustrate 



# 


# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge funding from ERC Consolidator Grant 771143, from the CHEESE-2p Center of Excellence (funded by EuroHPC-JU) and by the BMBF project DEGREE.

# References
