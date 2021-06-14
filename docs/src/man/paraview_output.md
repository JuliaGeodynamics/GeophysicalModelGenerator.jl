# Paraview output

We have one main routine to generate Paraview output for data that is either stored in a `GeoData` structure (that has lat/lon info), or `CartData` (Cartesian).
If `GeoData` is supplied it is internally automatically converted to the right format. Vectors, such as velocity, are also converted accordingly.

```@docs
GeophysicalModelGenerator.Write_Paraview
```
