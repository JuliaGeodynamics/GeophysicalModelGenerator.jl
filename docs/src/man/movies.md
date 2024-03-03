# Create movies

We have two routines to help create movies from the results.
The first one takes `*.png` images generated with the `Save Animation` option in paraview and puts them together in compressed movies (either `mp4` or `mov`):
```@docs
movie_from_images
```

The other one creates `*.pvd` files that can be saved with the `pvd=...` optional option in `Write_Paraview`, such that you can animate temporal data in paraview (yif you're happy you can save the result as images and use `movies_from_pics`). See the corresponding tutorial on how to generate `*.pvd` files.
```@docs
Write_Paraview
```
