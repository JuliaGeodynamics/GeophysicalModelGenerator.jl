# LaMEM I/O

In order to generate geodynamic simulations from setups created with `GeophysicalModelGenerator.jl`, we provide a few routines that directly create marker input files for the 3D geodynamic modelling software [LaMEM](https://bitbucket.org/bkaus/lamem), which is an open-source cartesian code that is well-suited to perform crustal and lithospheric-scale simulations. 
If you want to learn how to run LaMEM simulations, please have a look at the [wiki page](https://bitbucket.org/bkaus/lamem/wiki/Home). 

```@docs
GeophysicalModelGenerator.ReadLaMEM_InputFile
GeophysicalModelGenerator.Save_LaMEMMarkersParallel
GeophysicalModelGenerator.GetProcessorPartitioning
GeophysicalModelGenerator.LaMEM_grid
GeophysicalModelGenerator.ReadData_PVTR
```
