# Geodynamic model setups

In order to generate geodynamic simulations from setups created with `GeophysicalModelGenerator.jl`, we provide a few routines that directly create setups 

The routines provided here have the following functionality:
- Add lithospheric boxes to a setup, that may have a layered structure and various thermal structures

```@docs
GeophysicalModelGenerator.AddBox!
GeophysicalModelGenerator.AddLayer!
GeophysicalModelGenerator.AddSphere!
GeophysicalModelGenerator.AddEllipsoid!
GeophysicalModelGenerator.AddCylinder!
GeophysicalModelGenerator.makeVolcTopo
GeophysicalModelGenerator.ConstantTemp
GeophysicalModelGenerator.LinearTemp
GeophysicalModelGenerator.HalfspaceCoolingTemp
GeophysicalModelGenerator.SpreadingRateTemp
GeophysicalModelGenerator.LithosphericTemp
GeophysicalModelGenerator.ConstantPhase
GeophysicalModelGenerator.Compute_Phase
GeophysicalModelGenerator.LithosphericPhases
```
