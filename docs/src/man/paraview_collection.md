# Paraview collection

We have one main routine to generate `*.pvd` files from existing `vtk` files. This is useful if no `*.pvd` file was generated during the simulation, or if you want to generate a `*.pvd` file from a collection of `vtk` files that were generated in different simulations. The `*.pvd` file can be used to animate temporal data in paraview. You can either provide a Vector of the files or the specific time step or the function reads the directory and assigns a pseudo time step to the `*.pvd` file.
```@docs
GeophysicalModelGenerator.make_paraview_collection
```
