# # Getting started
# ### Aim:
# The aim of this tutorial is to give you a feeling for the `GeophysicalModelGenerator` package.
#
# ### 1. Loading data
# We start with loading the package, and assume that you have installed it:
using GeophysicalModelGenerator;

# Lets load a first 3D seismic tomography data set of the Alps. We uploaded an example file on Zenodo [here](https://zenodo.org/records/10738510), which comes from the paper of [Paffrath et al. 2021](https://se.copernicus.org/articles/12/2671/2021/se-12-2671-2021.html).
# We can load this in `GMG` with:
Tomo_Alps_full = load_GMG("https://zenodo.org/records/10738510/files/Paffrath_2021_SE_Pwave.jld2?download=1")

# This is a so-called `GeoData` object, which is a 3D grid of seismic velocities as a function of longitude, latitude and depth, which can include various fields (here we only have a single field: `:dVp_Percentage`)
# We can save this in `VTK` format, which is a widely used format that can for exampke be read by the 3D open-source visualization tool [Paraview](https://www.paraview.org/):
Write_Paraview(Tomo_Alps_full,"Tomo_Alps_full")

# We also uploaded a dataset with the topography of the Alpine region which can be downloaded with:
Topo_Alps = load_GMG("https://zenodo.org/records/10738510/files/AlpsTopo.jld2?download=1")

# Different than the 3D tomographic model, the topography has size 1 for the last index which indicates that this is a 3D surface. As you can see, the depth varies, which is because it is a warped surface.
# We can write this to disk as well
Write_Paraview(Topo_Alps,"Topo_Alps")

# If we open both datasets in Paraview, we see this (after giving some color to the topography):
# ![Basic_Tutorial_1](../assets/img/Basic_Tutorial_1.png)
# Note that I use the `Oleron` scientific colormap for the tomography which you can download [here](https://www.fabiocrameri.ch/colourmaps/)
#
# ### 2. Extract subset of data
# As you can see the tomographic data covers a much larger area than the Alps itself, and in most of that area there is no data.
# It is thus advantageous to cut out a piece of the dataset that we are interested in which can be done with `ExtractSubVolume`:
Tomo_Alps = ExtractSubvolume(Tomo_Alps_full,Lon_level=(4,20),Lat_level=(36,50), Depth_level=(-600,-10))

Write_Paraview(Tomo_Alps,"Tomo_Alps");

# Which looks like:
# ![Basic_Tutorial_2](../assets/img/Basic_Tutorial_2.png)

# ### 3. Create cross sections
# Paraview has the option to `Slice` through the data but it is not very intuitive to do this in 3D. Another limitation of Paraview is that it does not have native support for spherical coordinates, and therefore the data is translated to cartesian (`x`,`y`,`z`) coordinates (with the center of the Earth at `(0,0,0)`).
# That makes this a bit cumbersome to make a cross-section at a particular location.
# If you are interested in this you can use the `CrossSection` function:
data_200km = CrossSection(Tomo_Alps, Depth_level=-200)
#
# As you see, this is not exactly at 200 km depth, but at the closest `z`-level in the data sets. If you want to be exactly at 200 km, use the `Interpolate` option:
data_200km_exact = CrossSection(Tomo_Alps, Depth_level=-200, Interpolate=true)

# In general, you can get help info for all functions with `?`:
# ```julia
# help?> CrossSection
# search: CrossSection CrossSectionVolume CrossSectionPoints CrossSectionSurface FlattenCrossSection
# 
#   CrossSection(DataSet::AbstractGeneralGrid; dims=(100,100), Interpolate=false, Depth_level=nothing, Lat_level=nothing, Lon_level=nothing, Start=nothing, End=nothing, Depth_extent=nothing, section_width=50km)
# 
#   Creates a cross-section through a GeoData object.
# 
#     •  Cross-sections can be horizontal (map view at a given depth), if Depth_level is specified
# 
#     •  They can also be vertical, either by specifying Lon_level or Lat_level (for a fixed lon/lat), or by defining both Start=(lon,lat) & End=(lon,lat) points.
# 
#     •  Depending on the type of input data (volume, surface or point data), cross sections will be created in a different manner:
# 
#     1. Volume data: data will be interpolated or directly extracted from the data set.
# 
#     2. Surface data: surface data will be interpolated or directly extracted from the data set
# 
#     3. Point data: data will be projected to the chosen profile. Only data within a chosen distance (default is 50 km) will be used
# 
#     •  Interpolate indicates whether we want to simply extract the data from the data set (default) or whether we want to linearly interpolate it on a new grid, which has dimensions as specified in dims NOTE: THIS ONLY APPLIES TO VOLUMETRIC AND SURFACE DATA
#        SETS
# 
#     •  'section_width' indicates the maximal distance within which point data will be projected to the profile
# 
#   Example:
#   ≡≡≡≡≡≡≡≡
# 
#   julia> Lon,Lat,Depth   =   LonLatDepthGrid(10:20,30:40,(-300:25:0)km);
#   julia> Data            =   Depth*2;                # some data
#   julia> Vx,Vy,Vz        =   ustrip(Data*3),ustrip(Data*4),ustrip(Data*5);
#   julia> Data_set3D      =   GeoData(Lon,Lat,Depth,(Depthdata=Data,LonData=Lon, Velocity=(Vx,Vy,Vz))); 
#   julia> Data_cross      =   CrossSection(Data_set3D, Depth_level=-100km)  
#   GeoData 
#     size  : (11, 11, 1)
#     lon   ϵ [ 10.0 : 20.0]
#     lat   ϵ [ 30.0 : 40.0]
#     depth ϵ [ -100.0 km : -100.0 km]
#     fields: (:Depthdata, :LonData, :Velocity)
# ```
#
# Let's use this to make a vertical cross-section as well:
Cross_vert = CrossSection(Tomo_Alps, Start=(5,47), End=(15,44))

# And write them to paraview:
Write_Paraview(Cross_vert,"Cross_vert");
Write_Paraview(data_200km,"data_200km");

# ![Basic_Tutorial_3](../assets/img/Basic_Tutorial_3.png)
# In creating this image, I used the `Clip` tool of Paraview to only show topography above sealevel and made it 50% transparent.
#
# ### 4. Cartesian data
# As you can see, the curvature or the Earth is taken into account here. Yet, for many applications it is more convenient to work in Cartesian coordinates (kilometers) rather then in geographic coordinates.
# `GeophysicalModelGenerator` has a number of tools for this.
# First we need do define a `ProjectionPoint`  around which we project the data
proj = ProjectionPoint(Lon=12.0,Lat =43)

Topo_cart = Convert2CartData(Topo_Alps, proj)

# And do the same with the tomography:
Tomo_cart = Convert2CartData(Tomo_Alps, proj)

# Save:
Write_Paraview(Tomo_cart,"Tomo_cart");
Write_Paraview(Topo_cart,"Topo_cart");

# ![Basic_Tutorial_4](../assets/img/Basic_Tutorial_4.png)
# As the coordinates are now aligned with the `x`,`y`,`z` coordinate axes in Paraview it is now straightforward to use the build-in tools to explore the data.

# ### 5. Rectilinear data
#=
Yet, because of the curvature of the Earth, the resulting 3D model is not strictly rectilinear, which is often a requiment for cartesian numerical models.
This can be achieved in a relatively straightforward manner, by creating a new 3D dataset that is slightly within the curved boundaries of the projected data set:
=#

Tomo_rect = CartData(XYZGrid(-550.0:10:600, -500.0:10:700, -600.0:5:-17)); 
#
# the routine `ProjectCartData` will then project the data from the geographic coordinates to the new rectilinear grid:
Tomo_rect = ProjectCartData(Tomo_rect, Tomo_Alps, proj)

# we can do the same with topography:
Topo_rect = CartData(XYZGrid(-550.0:1:600, -500.0:1:700, 0)) 
Topo_rect = ProjectCartData(Topo_rect, Topo_Alps, proj)

# Save it:
Write_Paraview(Tomo_rect,"Tomo_rect");
Write_Paraview(Topo_rect,"Topo_rect");
# ![Basic_Tutorial_5](../assets/img/Basic_Tutorial_5.png)

# At this stage, the data can directly be used to generate cartesian numerical model setups, as explained in the other tutorials. 


#src Note: The markdown page is generated using:
#src Literate.markdown("tutorials/Tutorial_Basic.jl","docs/src/man",keepcomments=true, execute=true, codefence = "```julia" => "```")

