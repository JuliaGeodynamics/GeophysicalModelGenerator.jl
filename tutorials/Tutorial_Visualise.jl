#=
# Quick look with GLMakie

### Aim
Writing Paraview files and opening them is the most flexible way to look at data, but often we only want a quick
look without leaving Julia. When `GLMakie` is loaded together with `GeophysicalModelGenerator`, two tools become
available for this:

- `visualise`, an interactive 3D widget that shows slices through a 3D data set, optionally with topography on top;
- `heatmap` and `heatmap!` methods that plot a horizontal 2D data set, such as topography or a map-view cross-section.

This tutorial shows both for seismic tomography of the Alps and for a numerical model setup.
It needs a computer with a display, as `GLMakie` opens a window.

### 1. Loading the data
We load both packages:
=#
using GeophysicalModelGenerator, GLMakie

#=
We use the P-wave tomography of [Paffrath et al. 2021](https://se.copernicus.org/articles/12/2671/2021/se-12-2671-2021.html)
and topography of the Alps, which are also used in the Getting started tutorial. Both are downloaded from
[Zenodo](https://zenodo.org/records/10738510). We keep the region around the Alps:
=#
Tomo_Alps_full = load_GMG("https://zenodo.org/records/10738510/files/Paffrath_2021_SE_Pwave.jld2?download=1")
Tomo_Alps = extract_subvolume(Tomo_Alps_full, Lon_level = (4, 20), Lat_level = (36, 50), Depth_level = (-600, -10))
Topo_Alps = load_GMG("https://zenodo.org/records/10738510/files/AlpsTopo.jld2?download=1")

#=
### 2. Exploring a 3D GeoData set
`visualise` opens the widget for a 3D `GeoData` set. It returns the figure, which you can save with `save`:
=#
fig_tomo = visualise(Tomo_Alps)

# ![Tutorial_Visualise_GeoData](../assets/img/Tutorial_Visualise_GeoData.png)

#=
The first call takes a while because Julia compiles the plotting code; later calls are fast.

The widget has the following controls:

- **Dataset**: the field that is shown. All fields of the data set are listed.
- **Colormap**: the colormap of the slices.
- **Color axis limits**: the lower and upper limit of the colormap. Type a number and press `Enter`.
  A lower limit that is not below the upper limit is rejected, and the box turns red.
- **Isovalue**: the toggle shows or hides an isosurface of the field, the text box sets its value
  (press `Enter`), and the `α` slider sets its opacity.
- **lon**, **lat**, **depth**: the sliders move the three slices through the volume. The toggle next to each slider
  shows or hides that slice.
- **Transparency topo**: the opacity of the topography, if given. The toggle next to it shows or hides the topography.

In the 3D view, drag with the left mouse button to rotate, drag with the right mouse button to move, and scroll to zoom.

A topography is drawn on top of the data set with the `Topography` keyword:
=#
fig_tomo_topo = visualise(Tomo_Alps, Topography = Topo_Alps)

# ![Tutorial_Visualise_GeoData_Topo](../assets/img/Tutorial_Visualise_GeoData_Topo.png)

#=
By default, the topography colormap is centered around zero and spans the highest point. Use `Topo_range`
to set the limits yourself, e.g. `visualise(Tomo_Alps, Topography = Topo_Alps, Topo_range = (-4, 4))`.

### 3. Exploring Cartesian data
`GeoData` is shown in lon/lat/depth, which distorts distances. For a view in kilometers, we project the data onto
a Cartesian grid with `project_CartData`. `visualise` needs an orthogonal grid, so we create one that lies inside
the projected data (see the Getting started tutorial):
=#
proj = ProjectionPoint(Lon = 12.0, Lat = 43)
Tomo_rect = CartData(xyz_grid(-550.0:10:600, -500.0:10:700, -600.0:5:-17))
Tomo_rect = project_CartData(Tomo_rect, Tomo_Alps, proj)
Topo_rect = CartData(xyz_grid(-550.0:5:600, -500.0:5:700, 0))
Topo_rect = project_CartData(Topo_rect, Topo_Alps, proj)

#=
For `CartData`, the axes have equal scaling, so the model has its true shape:
=#
fig_cart = visualise(Tomo_rect, Topography = Topo_rect)

# ![Tutorial_Visualise_CartData](../assets/img/Tutorial_Visualise_CartData.png)

#=
### 4. Quick 2D plots with heatmap
For 2D horizontal data sets, which have one point in the vertical direction, `heatmap` plots a field in a new figure.
The default field is `:Topography`:
=#
heatmap(Topo_Alps)

#=
Use the `field` keyword for other fields, e.g. for a map-view cross-section through the tomography at 200 km depth:
=#
Tomo_200km = cross_section(Tomo_Alps, Depth_level = -200)
fig_200km, ax_200km, hm_200km = heatmap(Tomo_200km, field = :dVp_Percentage, colormap = :roma, colorrange = (-3, 3))
Colorbar(fig_200km[1, 2], hm_200km, label = "dVp [%]")

# ![Tutorial_Visualise_heatmap](../assets/img/Tutorial_Visualise_heatmap.png)

#=
`heatmap!` adds the plot to the current axis. For `CartData`, `colorbar = true` adds a colorbar:
=#
fig_topo = Figure()
Axis(fig_topo[1, 1], aspect = DataAspect(), xlabel = "x [km]", ylabel = "y [km]")
heatmap!(Topo_rect, colormap = :oleron, colorrange = (-4, 4), colorbar = true)

# ![Tutorial_Visualise_heatmap_cart](../assets/img/Tutorial_Visualise_heatmap_cart.png)

#=
`heatmap` and `heatmap!` also accept an array instead of a field, e.g. `heatmap(Topo_Alps, Topo_Alps.depth.val)`.

### 5. Checking a numerical model setup
`visualise` is also a quick way to check a model setup before writing Paraview or LaMEM files. We build a simple
3D subduction setup with an oceanic plate and an inclined slab (see the 3D model setups tutorial for details):
=#
nx, ny, nz = 128, 128, 64
Grid = CartData(xyz_grid(range(-1000, 1000, nx), range(-1000, 1000, ny), range(-660, 0, nz)))
Phases = fill(2, nx, ny, nz)
Temp = fill(1350.0, nx, ny, nz)
lith = LithosphericPhases(Layers = [15 45 10], Phases = [0 1 2], Tlab = 1250)
add_box!(
    Phases, Temp, Grid; xlim = (-800.0, 0.0), ylim = (-400, 400.0), zlim = (-80.0, 0.0), phase = lith,
    T = SpreadingRateTemp(SpreadingVel = 3, MORside = "right"), StrikeAngle = 30
)
add_box!(
    Phases, Temp, Grid; xlim = (0.0, 300.0), ylim = (-400.0, 400.0), zlim = (-80.0, 0.0), phase = lith,
    T = HalfspaceCoolingTemp(Age = 800.0e3 / (3 / 100) / 1.0e6), DipAngle = 30, StrikeAngle = 30
)
Grid = addfield(Grid, (; Phases, Temp))

#=
The widget opens with the first field, `Z`, which `xyz_grid` adds. Select `Phases` in the **Dataset** menu and move
the `X [km]` slider to check that the slab has the intended position and dip; select `Temp` to check the thermal
structure. The screenshot shows `Phases` with the slice at `x = 150` km:
=#
fig_model = visualise(Grid)

# ![Tutorial_Visualise_model](../assets/img/Tutorial_Visualise_model.png)

#=
Once the setup looks right, write it to disk, for example with `write_paraview(Grid, "Subduction3D")` or with
`save_LaMEM_markers_parallel`.
=#
