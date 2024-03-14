---
title: 'GeophysicalModelGenerator.jl: A Julia package to visualise geoscientific data and create numerical model setups'
tags:
  - julia
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
  - name: Pascal Aellig
    orcid: 0009-0008-9039-5646
    affiliation: 1
  - name: Albert de Montserrat
    orcid: 0000-0003-1694-3735 
    affiliation: 3
  - name: Luca de Siena
    orcid: 0000-0002-3615-5923
    affiliation: 4
  - name: Jacob Frasukiewicz
    orcid: 0009-0002-5049-4259
    affiliation: 1
  - name: Lukas Fuchs
    orcid: 0000-0002-9165-6384
    affiliation: 5
  - name: Andrea Piccolo
    orcid: 0000-0003-3074-6041
    affiliation: 2
  - name: Hendrik Ranocha
    orcid: 0000-0002-3456-2277
    affiliation: 1
  - name: Nicolas Riel
    orcid: 0000-0002-5037-5519
    affiliation: 1
  - name: Christian Schuler
    orcid: 0009-0004-9873-9774
    affiliation: 1
  - name: Arne Spang
    orcid: 0000-0002-6107-0403
    affiliation: 2
  - name: Tatjana Weiler
    #orcid: 0000-0002-6107-0403
    affiliation: 2

affiliations:
 - name: Johannes Gutenberg University Mainz, Germany
   index: 1
 - name: University of Bayreuth, Germany
   index: 2
 - name: ETH Zürich, Switzerland
   index: 3
 - name: Bologna University, Italy
   index: 4
 - name: Goethe University Frankfurt, Germany
   index: 5

date: 14 March 2024
bibliography: paper.bib
---

# Summary

Geoscientific data exists in wide variety of different formats. Yet, to make a consistent interpretation of a certain region, it is often helpful to jointly visualise all this data using the same coordinates, and compare, for example seismic tomography, surface geology, Moho depth, Earthquake locations, and GPS surface velocities. If one wishes to create mechanical or thermo-mechanical numerical models of the region, creating an input model that honors these  constraints is helpful. And since most numerical codes work in Cartesian boxes, it is helpful to have tools to project the data from geographic to Cartesian coordinates.

A significant challenge in doing this is that there is no standard format for geoscientific data. Seismic tomography, for example, may come in the form of ASCII data with `lon/lat/depth` axes, or as NetCDF files, with the ordering of the data typically differing from one dataset to the other. In ideal cases, geological surfaces may be provided as GeoTIFF images. In many cases, however, the underlying data discussed in publications are not available in digital format and only shown as figures in the paper. It is nevertheless still helpful to image these in 3D in the correct coordinates, along with more recent, digitally available, datasets.

The aim of the `GeophysicalModelGenerator.jl` package is therefore two-fold:

1) Simplify collecting and visualising a wide variety of geoscientific data that is provided as point (e.g., earthquake locations), surface (e.g., topography) or volumetric data (e.g., seismic tomography).  
2) Create input setups for 2D or 3D numerical models.

# Statement of need

`GeophysicalModelGenerator.jl` is a Julia [@bezanson2017julia] package that helps collecting and visualising a wide variety of geophysical and geoscientific data in a coherent manner. It also simplifies the process of generating 2D or 3D models that can, for example, be used as input models in geodynamic simulations. It provides functions that transfer data from one format to the other, or project them from geographic `Longitude/Latitude/Depth` or `UTM` coordinates to Cartesian coordinates (`kilometers`). It allows performing tasks such as creating cross-sections though volumetric data, importing screenshots from published papers, downloading digital elevation data and saving the resulting data in `VTK` format, that can, for instance, be visualised with open source tools such as [Paraview](www.paraview.org).

Most geoscientists tend to have their own python/matlab/bash visualisation and therefore perform part of this job already. Yet, having all functionality in one place in an easy-to-use package, makes this more extendable and will likely facilitate sharing data along with their interpretations. 

# Related software packages
Perhaps the most widely used package in geophysics to create figures or maps is the Generic Mapping Tools package ([`GMT`](https://www.generic-mapping-tools.org)), which also provides a Julia interface [GMT.jl](https://github.com/GenericMappingTools/GMT.jl) [@Wessel_Luis_Uieda_Scharroo_Wobbe_Smith_Tian_2019]. It mostly focuses on generating maps and postscript/pdf images and is therefore not ideally suited for interactive 3D data visualisation or to generate input models for numerical codes.
 
The [`Geodynamic World Builder`](https://github.com/GeodynamicWorldBuilder/WorldBuilder) is a C++ library to create model setups [@se-10-1785-2019]. The focus is on generating input models for geodynamic simulations, such as subduction zones and related thermal structures. It has C and Fortran wrappers and can thus be embedded in geodynamic codes. Users of the `Geodynamic World Builder` have to generate JSON files to define the model geometry, which is less interactive than by using the Julia `REPL`. There is no straightforward way to integrate existing geophysical/geological data in the workflow and compare model results with them.   

[GemPy](https://www.gempy.org) is a Python-based, open-source geomodeling library that can construct 3D geological models of folded structures, fault networks and unconformities, while taking uncertainties into account [@DeLaVarga_Schaaf_Wellmann_2019]. It focuses on creating geometric models with uncertainties rather than on integrating a wide variety of geoscientific datasets.

There are also a number of commercial software solutions: 

- [Petrel subsurface software](https://www.software.slb.com/products/petrel) (by Schlumberger), which is mostly used by the hydrocarbon industry and is particularly powerful in integrating seismic reflection and well-data.

- [GOCAD Mining Suite](https://www.mirageoscience.com/mining-industry-software/gocad-mining-suite/) (by MiraGeoscience) helps generating geometric models of the sub surface in the vicinity of mines, based on sparse geological measurements and drillhole data.

- [GeoModeller](https://www.intrepid-geophysics.com/products/geomodeller/) (by Intrepid Geophysics) creates surface-near geometric geological models by implicit modelling of surface measurements while taking geophysical constraints into account.    

In all cases, the commercial license fees are far beyond what most researchers can afford, even if reduced license fees are often available for academia. The closed-source nature of the software packages makes them non-extendable by the community.

The `GeophysicalModelGenerator.jl` julia package is already used to generate input models for the geodynamic codes [LaMEM](https://github.com/UniMainzGeo/LaMEM) [@Kaus_Popov_Baumann_Pusok_Bauville_Fernandez_Collignon_2016], [JustRelax.jl](https://github.com/PTsolvers/JustRelax.jl), and [MagmaThermokinematics.jl](https://github.com/boriskaus/MagmaThermoKinematics.jl). It is also already used in a number of shortcourses and lectures at the universities of Mainz, Heidelberg, and Bologna.

# Basic usage

The core of the package consists of the  `GeoData`, `UTMData`, `ParaviewData`, and `CartData` structures which hold the 3D data along with coordinates (and potentially metadata) information. `GeophysicalModelGenerator.jl` can be installed using the build-in Julia package manager:
```julia
julia> ]
(@v1.10) pkg> add GeophysicalModelGenerator
```
which comes with a test-suite:
```julia
(@v1.10) pkg> test GeophysicalModelGenerator
```
and can be loaded with:
```julia
julia> using GeophysicalModelGenerator
```

As a first example, we will download a 3D seismic tomography dataset for the Alpine region (from @Paffrath_Friederich_Schmid_Handy_2021):
```julia
julia> Tomo_Alps_full = load_GMG(
  "https://zenodo.org/records/10738510/files/Paffrath_2021_SE_Pwave.jld2?download=1");
```
We can download the topography for the Alpine region with:
```julia
julia> Topo_Alps = load_GMG(
  "https://zenodo.org/records/10738510/files/AlpsTopo.jld2?download=1");
```

The seismic data covers a much wider region than the Alps itself, but in much of that region there is poor data coverage. We can therefore extract a part of the data that has coverage:
```julia
julia> Tomo_Alps = extract_subvolume(Tomo_Alps_full, Lon_level=(4,20), 
                        Lat_level=(36,50), Depth_level=(-600,-10));
```

At this stage, we can save the data to `VTK`  format:
```julia
julia> write_paraview(Tomo_Alps,"Tomo_Alps");
julia> write_paraview(Topo_Alps,"Topo_Alps");
```
and open it with Paraview (see \autoref{fig:basic}a).
We can create vertical and horizontal cross-sections through the data with:
```julia
julia> Cross_200km = cross_section(Tomo_Alps, Depth_level=-200, Interpolate=true);
julia> Cross_vert  = cross_section(Tomo_Alps, Start=(5,47), End=(15,44));
julia> write_paraview(Cross_vert, "Cross_vert");
julia> write_paraview(Cross_200km,"Cross_200km");
```
and visualise them along with the volumetric data (\autoref{fig:basic}a).

![Example of visualising 3D seismic data of the Alps, using a) geographic coordinates (`GeoData`) or b) Cartesian coordinates (`CartData`) projected from geographic coordinates. Shown are topography as well as several slices through the 3D seismic tomography P-wave model of [@Paffrath_Friederich_Schmid_Handy_2021].  \label{fig:basic} ](Basic_Tutorial.png){ width=100% }

One complication with geographic data is that Paraview does not have native support for geographic coordinates, and accordingly it is not always straightforward to use the build-in tools, for example, to create slices through the data. 
In addition, many numerical models work in (orthogonal) Cartesian rather than in spherical coordinates, which appears to be a good first-order approximation for many geodynamic applications [@Macherel_Räss_Schmalholz_2024].

`GeophysicalModelGenerator.jl` includes tools to transfer the data from geographic to Cartesian coordinates, which requires defining a projection point, along which the projection is performed:
```julia
julia> proj = ProjectionPoint(Lon=12.0,Lat =43)
ProjectionPoint(43.0, 12.0, 255466.98055255096, 4.765182932801006e6, 33, true)
```
We can now project the topography with:
```julia
julia> Topo_cart = convert2CartData(Topo_Alps, proj);
```
which returns a `CartData` (Cartesian data) structure. The disadvantage of doing this projection is that the resulting Cartesian grid is no longer strictly orthogonal which is a problem for some Cartesian numerical models (e.g., using finite difference discretisations).
We can project the data on an orthogonal grid as well, by first creating appropriately sized orthogonal grids for the tomography and topography:
```julia
julia> Tomo_rect = CartData(xyz_grid(-550.0:10:600, -500.0:10:700, -600.0:5:-17));
julia> Topo_rect = CartData(xyz_grid(-550.0:1:600, -500.0:1:700, 0)); 
```
Next, we can project the data to the orthogonal grids with:
```julia
julia> Topo_rect = project_CartData(Topo_rect, Topo_Alps, proj);
julia> Tomo_rect = project_CartData(Tomo_rect, Tomo_Alps, proj);
julia> write_paraview(Tomo_rect,"Tomo_rect");
julia> write_paraview(Topo_rect,"Topo_rect");
```
We can now use the build-in tools of Paraview to visualise the data (see \autoref{fig:basic} b), and use this as inspiration to create an initial numerical model setup. It is also possible to interpolate other seismic tomography datasets to the same grid and subsequently compute a "votemap" to count in how many tomographic models a specific seismic anomaly is present.


# Examples of usage
`GeophysicalModelGenerator.jl` comes with build-in (CI/CD) tests and [tutorials](https://juliageodynamics.github.io/GeophysicalModelGenerator.jl/stable) that explain the most important use cases, from importing data to generating input model setups for numerical simulations. In the following, we present a number of examples that illustrate various aspects of the package.

### Visualise data of the Alps
The European Alps are among the best studied mountain belts on the planet, and have therefore been the focus of numerous geological and geophysical studies. Different seismic tomography model have been published (using different parameterisations and datasets), and those do not necessarily agree with each other. 

In `Tutorial_AlpineData.jl`, users learn how to load the topography of the region, import Moho data, load and visualise GPS vectors, import and plot earthquake locations, along with cross-sections through the model (\autoref{fig:alps}).

![Example of combined data of the Alps, which shows the GPS surface velocity (arrows), topography, earthquake locations (colored dots) and cross-sections through a recent anisotropic P-wave tomography model by [@Rappisi_VanderBeek_Faccenda_Morelli_Molinari_2022]. \label{fig:alps}](../src/assets/img/GMG_AlpineData.png){ width=90% }

### La Palma volcanic eruption
The 2019 Cumbre Viejo eruption in La Palma, Canary Islands, was accompanied by seismic activity. In `Tutorial_LaPalma.jl`, users learn to generate a Cartesian block model of the island, import seismicity and use that to generate a 3D volumetric seismic activity map (\autoref{fig:lapalma}). 

![Example of a model of La Palma which shows seismicity during the 2019 Cumbre Viejo eruption. \label{fig:lapalma}](../src/assets/img/Tutorial_LaPalma.png){ width=100% }


### Jura mountains
The Jura mountains are a small-scale fold and thrust belt located in the Switzerland and France. Thanks to seismic cross-sections and boreholes, quite a bit of information is available about its structure at depth, which was used to generate extensive 3D models of the subsurface including thickness maps of various geological units, generate a new geological map of the region, and create balanced reconstructions [@Schori_2021].  

In `Tutorial_Jura.jl` users learn how to drape the geological map over the topography, import surfaces from GeoTIFF images (such as basement topography), and include screenshots from geological cross-sections. The data is rotated and transferred to Cartesian coordinates such that we obtain a 3D block model that is perpendicular to the strike of the mountain range (\autoref{fig:jura}).

![Example of creating a 3D Cartesian block model that runs perpendicular to the Jura mountains, combining surface geology, with screenshots from interpreted cross-sections (gray drawing), and digital data of the the basement topography [using data of @Schori_2021]. \label{fig:jura}](../src/assets/img/Jura_2.png){ width=100% }


### Slab model setup 
In `Tutorial_NumericalModel_3D.jl`, users learn how to generate a 3D geodynamic model setup with subducting slabs, a mid oceanic ridge and an overriding cratonic lithosphere. The thermal structure of the oceanic slab increases away from the ridge until the trench following a halfspace cooling analytical solution, whereas the thermal structure of the subducted part of the slab is based on an analytical solution that takes heating from the surrounding, hot, mantle into account [after @McKenzie_1969]. Between the mantle and the trench the slab uses a mixture between these two thermal models, and a weak-zone is added above the slab (to facilitate subduction in numerical models). A sedimentary wedge is situated at the continental margin and a grid-like pattern is put on top of the oceanic slab to simplify tracking deformation throughout the simulation (\autoref{fig:slab3d}).

![Example of a geodynamic setup of a subducting oceanic plate beneath a continental lithosphere with a smoothly bending slab. \label{fig:slab3d}](../src/assets/img/Mechanical3D_Tutorial_2.png){ width=100% }

# Acknowledgements
We acknowledge funding from ERC Consolidator Grant 771143 (MAGMA), by the German Ministry of Science and Education (BMBF) as part of project DEGREE, by the CHEESE-2p Center of Excellence (co-funded by both EuroHPC-JU and the BMBF), by the German Research Foundation (DFG grants TH2076/7-1 and KA3367/10-1) as  part of the SPP 2017 4DMB project project, and by DFG Emmy Noether grant TH2076/8-1.

# References
