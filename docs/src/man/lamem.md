# LaMEM

In order to generate geodynamic simulations from setups created with `GeophysicalModelGenerator.jl`, we provide a few routines that directly create marker input files for the 3D geodynamic modelling software [LaMEM](https://github.com/UniMainzGeo/LaMEM), which is an open-source cartesian code to perform crustal and lithospheric-scale simulations. 
If you want to learn how to run LaMEM simulations, the easiest way to get started is by looking at [LaMEM.jl](https://github.com/JuliaGeodynamics/LaMEM.jl) which is integrated with `GMG`

The routines provided here have the following functionality:
- Read LaMEM *.dat files (to get the size of the domain)
- Read LaMEM processor partitioning file
- Save LaMEM marker files in serial or in parallel
- Read a LaMEM timestep

```@docs
GeophysicalModelGenerator.ReadLaMEM_InputFile
GeophysicalModelGenerator.GetProcessorPartitioning
GeophysicalModelGenerator.Save_LaMEMTopography
GeophysicalModelGenerator.Save_LaMEMMarkersParallel
GeophysicalModelGenerator.ReadData_PVTR
GeophysicalModelGenerator.LaMEM_grid
GeophysicalModelGenerator.CreatePartitioningFile
```
