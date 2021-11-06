# Data structures

The main data structure used in GeophysicalModelGenerator.jl is `GeoData`, which contains info about the `longitude`,`latitude`, and `depth` of a data set, as well as several data sets itself.
We also provide a `UTMData`, which is essentially the same but with UTM coordinates, and a `CartData` structure, which has Cartesian coordinates in kilometers (as used in many geodynamic codes).
For plotting, we transfer this into the `ParaviewData` structure, which has cartesian coordinates around the center of the Earth. We wmploy the `wgs84` reference ellipsoid as provided by the [Geodesy.jl](https://github.com/JuliaGeo/Geodesy.jl) package to perform this transformation.

```@docs
GeophysicalModelGenerator.GeoData
GeophysicalModelGenerator.UTMData
GeophysicalModelGenerator.ParaviewData
GeophysicalModelGenerator.CartData
GeophysicalModelGenerator.LonLatDepthGrid
GeophysicalModelGenerator.XYZGrid
```
