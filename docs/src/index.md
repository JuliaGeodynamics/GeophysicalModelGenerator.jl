```@meta
CurrentModule = GeophysicalModelGenerator
```

# GeophysicalModelGenerator

Documentation for [GeophysicalModelGenerator](https://github.com/JuliaGeodynamics/GeophysicalModelGenerator.jl).

The main purpose of this package is to simplify the process of going from 1D/2D/3D geophysical data to a 3D consistent model of the region. By simplifying the process of plotting the data, it becomes easier to compare different data sets, and generate a 3D models that can be used for other computations such as geodynamic simulations, or forward modelling of gravity anomalies.

For this, GeophysicalModelGenerator provides the following functionality:
- A consistent GeoData structure, that holds the data along with lon/lat/depth information. 
- Routines to generate VTK files from the GeoData structure in order to visualize results in Paraview.
- The ability to deal with points, 2D profiles and 3D volumes, for both scalar and vector values.
- Rapidly import screenshots of published papers and compare them with other data sets in 3D using paraview.
- Create movies for representation of earthquake or wave propagation.
- Create geodynamic input models (for LaMEM)

The best way to get started is to look at the tutorials.
