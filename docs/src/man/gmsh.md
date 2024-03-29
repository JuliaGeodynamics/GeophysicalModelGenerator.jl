# Gmsh

[Gmsh](https://gmsh.info) is a widely used 3D unstructured mesh generator that produces tetrahedral meshes which can be tagged by region. It is possible to import such meshes and use that to set region info in a `GMG` data set, or use it to generate pTatin input.

```@docs
GeophysicalModelGenerator.import_Gmsh
```