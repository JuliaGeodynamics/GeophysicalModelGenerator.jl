# LaMEM

In order to generate geodynamic simulations from setups created with `GeophysicalModelGenerator.jl`, we provide a few routines that directly create marker input files for the 3D geodynamic modelling software [LaMEM](https://bitbucket.org/bkaus/lamem), which is an open-source cartesian code that is well-suited to perform crustal and lithospheric-scale simulations. 
If you want to learn how to run LaMEM simulations, please have a look at the [wiki page](https://bitbucket.org/bkaus/lamem/wiki/Home). 

The routines provided here have the following functionality:
- Read LaMEM *.dat files (to get the size of the domain)
- Read LaMEM processor partitioning file
- Add lithospheric boxes to a setup, that may have a layered structure and various thermal structures
- Save LaMEM marker files in serial or in parallel
- Read a LaMEM timestep

```@docs
GeophysicalModelGenerator.ReadLaMEM_InputFile
GeophysicalModelGenerator.GetProcessorPartitioning
GeophysicalModelGenerator.Save_LaMEMTopography
GeophysicalModelGenerator.Save_LaMEMMarkersParallel
GeophysicalModelGenerator.ReadData_PVTR
GeophysicalModelGenerator.AddBox!
GeophysicalModelGenerator.ConstantTemp
GeophysicalModelGenerator.LinearTemp
GeophysicalModelGenerator.HalfspaceCoolingTemp
GeophysicalModelGenerator.SpreadingRateTemp
GeophysicalModelGenerator.ConstantPhase
GeophysicalModelGenerator.LithosphericPhases
GeophysicalModelGenerator.LaMEM_grid
```
