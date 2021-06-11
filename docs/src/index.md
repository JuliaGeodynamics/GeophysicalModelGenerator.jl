```@meta
CurrentModule = GeophysicalModelGenerator
```

# GeophysicalModelGenerator

Documentation for [GeophysicalModelGenerator](https://github.com/JuliaGeodynamics/GeophysicalModelGenerator.jl).

The main purpose of this package is to simplify the process of going from 1D/2D/3D geophysical data to a 3D consistent models. By simplifying the process of plotting the data, it becomes easier to 

For this we provide the following functionality:
- A consistent GeoData data structure, that holds the data sets along with lat/lon/depth information. 
- Routines to generate VTK files from the GeoData structure in order to visualie results in Paraview.
- The ability to deal with points, 3D volumes, 2D profiles for both scalar and vector values.
