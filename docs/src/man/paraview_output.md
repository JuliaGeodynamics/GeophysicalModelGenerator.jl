# Paraview output

We have one main routine to generate Paraview output for data that is either stored in a `GeoData` structure (that has lat/lon info), or `ParaviewData`, `CartData` (Cartesian).
If `GeoData` is supplied it is internally automatically converted to the right format. Vectors, such as velocity, are also converted accordingly.
You can also visualize time-dependent data, or combine existing paraview files into a `*.pvd` paraview collection (that can be used to show a movie)
```@docs
write_paraview
movie_paraview
make_paraview_collection
```
