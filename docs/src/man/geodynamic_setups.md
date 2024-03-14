# Numerical model setups

In order to generate numerical simulations from setups created with `GeophysicalModelGenerator.jl`, we provide a few routines that directly create setups. 

The routines provided here have the following functionality:
- Add lithospheric boxes to a setup, that may have a layered structure and various thermal structures
- Add various geometries (spheres, cuylinders, ellipsoids)
- Add lithospheric structure
- Add various 1D thermal structures (and possibilities to combine them)

```@docs
GeophysicalModelGenerator.add_box!
GeophysicalModelGenerator.add_layer!
GeophysicalModelGenerator.add_sphere!
GeophysicalModelGenerator.add_ellipsoid!
GeophysicalModelGenerator.add_cylinder!
GeophysicalModelGenerator.add_stripes!
GeophysicalModelGenerator.add_slab!
GeophysicalModelGenerator.make_volc_topo
GeophysicalModelGenerator.ConstantTemp
GeophysicalModelGenerator.LinearTemp
GeophysicalModelGenerator.HalfspaceCoolingTemp
GeophysicalModelGenerator.SpreadingRateTemp
GeophysicalModelGenerator.LithosphericTemp
GeophysicalModelGenerator.ConstantPhase
GeophysicalModelGenerator.compute_phase
GeophysicalModelGenerator.LithosphericPhases
GeophysicalModelGenerator.McKenzie_subducting_slab
GeophysicalModelGenerator.LinearWeightedTemperature
```
