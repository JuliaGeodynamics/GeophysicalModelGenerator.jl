# Surfaces

We have a number of functions to deal with horizontal surfaces, which are defined as `GeoData` or `CartData` structures that have 3 dimensional `flat` coordinate array. Flat implies that the resolution is `n` by `m` by `1`, which is required to correctly display them in paraview.
The most straightforward surface is the topography (which you can obtain with `ImportTopo`).
Surfaces can be added and subtracted.

```@docs
drape_on_topo
fit_surface_to_points
aboveSurface
belowSurface
interpolateDataOnSurface
is_surface
remove_NaN_Surface!
```
