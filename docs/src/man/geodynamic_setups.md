# Numerical model setups

In order to generate numerical simulations from setups created with `GeophysicalModelGenerator.jl`, we provide a few routines that directly create setups. 

The routines provided here have the following functionality:
- Add lithospheric boxes to a setup, that may have a layered structure and various thermal structures
- Add various geometries (spheres, cuylinders, ellipsoids)
- Add lithospheric structure
- Add various 1D thermal structures (and possibilities to combine them)

Map-view outlines of plates, weak zones, faults and trenches can be drawn in any vector graphics program, such as
Inkscape, Affinity or Illustrator, exported as SVG and read with `read_outlines`, which returns them in drawing order:
closed ones as polygons for `add_plate!`, open ones as `Trace` of `add_fault!` and `Trench`. Name each shape (or its
layer) after the outline it represents, and draw a reference line whose end points have known model coordinates, or let the page span the model
domain. `write_outlines` stores outlines as a CSV file.

Layers with curved or perturbed interfaces, listric faults, faults with a curved trace and Gaussian thermal
anomalies are shown in the tutorial [Curved layers, curved faults and thermal anomalies](@ref).

## Functions

```@docs
GeophysicalModelGenerator.add_box!
GeophysicalModelGenerator.add_layer!
GeophysicalModelGenerator.add_sphere!
GeophysicalModelGenerator.add_ellipsoid!
GeophysicalModelGenerator.add_cylinder!
GeophysicalModelGenerator.add_stripes!
GeophysicalModelGenerator.add_slab!
GeophysicalModelGenerator.add_fault!
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
GeophysicalModelGenerator.GaussianTemp
GeophysicalModelGenerator.Trench
GeophysicalModelGenerator.Tear
GeophysicalModelGenerator.BezierSlab
GeophysicalModelGenerator.add_plate!
GeophysicalModelGenerator.add_polygon!
GeophysicalModelGenerator.add_volcano!
GeophysicalModelGenerator.compute_slab_surface
GeophysicalModelGenerator.compute_thermal_structure
GeophysicalModelGenerator.bezier_curve
GeophysicalModelGenerator.read_outlines
GeophysicalModelGenerator.write_outlines
```
