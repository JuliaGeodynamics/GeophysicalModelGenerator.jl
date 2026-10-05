# Create profiles for the AdriaArrayGeometryPicker

## Goal
The [AdriaArrayGeometryPicker](https://github.com/JuliaGeodynamics/AdriaArrayGeometryPicker.jl) is a graphical user interface to compare different geophysical datasets along profiles and to pick structures (e.g., slabs or the Moho) in them. It does not read the original datasets itself, but works with *profiles* that have been created beforehand with GeophysicalModelGenerator: every dataset is projected onto the profile, and the result is saved as a `.pgmg` file that is then opened in the picker.

This tutorial shows how to create these `.pgmg` files for a whole set of vertical and horizontal profiles at once. You need two text files:

1. a **profile file** that lists the profiles (e.g. `profiles_mixed.txt`),
2. a **dataset file** that lists the datasets to be projected onto the profiles (e.g. `Datasets_Adria.txt`).

Both are described in detail below. Some general background on profile processing is given in the [Profile Processing](profile_processing.md) section.

!!! note
    The AdriaArrayGeometryPicker currently relies on the `profile_processing_mt` branch of GeophysicalModelGenerator. If you install the picker, this branch is installed automatically.

## 1. The profile file
The profile file is a comma-separated text file. The first line is a header and is skipped. Every following line defines one profile, which can be either horizontal or vertical:

| Profile type | Columns | Meaning |
|:--|:--|:--|
| horizontal | `Num, Depth` | profile number and depth of the horizontal slice in km (negative below the surface) |
| vertical | `Num, LonStart, LatStart, LonEnd, LatEnd` | profile number and the longitude/latitude (in degrees) of the start and end point |

Horizontal and vertical profiles can be mixed in one file. This is the file `profiles_mixed.txt`, which contains two horizontal slices (at 200 and 300 km depth) and five vertical profiles across the Adriatic region:

```
Num,Depth
1,-200
2,-300
3,9.777223439354351,42.936916884179524,15.465423125301268,47.79656290565586
4,10.160585214657562,42.71827704599679,15.860636558746865,47.55888190360081
5,10.541291181660043,42.49821887191029,16.252277335405516,47.31986389797744
6,10.919360031561707,42.27677067866074,16.64038623036512,47.07954242435165
7,11.294811014000167,42.05396063315099,17.0250043232569,46.837950655777426
```

A few remarks:
- The header is only skipped, its content does not matter. For a file with only vertical profiles you can, e.g., use `Num,LonStart,LatStart,LonEnd,LatEnd`.
- A line with one value after `Num` is interpreted as a horizontal profile, a line with four values as a vertical profile. Any other number of values results in an error.
- The profiles are processed in the order in which they appear in the file; `Num` itself is not used.
- A vertical profile starts at (`LonStart`,`LatStart`). The distance along the profile (`x_profile`, in km) is measured from this point.

The file is read with `read_picked_profiles`, which returns a vector of `ProfileData` objects:
```julia
julia> using GeophysicalModelGenerator
julia> profiles = read_picked_profiles("profiles_mixed.txt");
julia> profiles[1]
Horizontal ProfileData
  depth      : -200.0
julia> profiles[3]
Vertical ProfileData
  lon/lat    : (9.777223439354351, 42.936916884179524)-(15.465423125301268, 47.79656290565586)
```

## 2. The dataset file
The dataset file is also a comma-separated text file, with one line per dataset. Again, the first line is a header and is skipped. Each line has the columns

| Column | Meaning |
|:--|:--|
| `Name` | Name of the dataset. It appears in the picker and is used as prefix of the field names in the profile (e.g. `Giacomuzzi_dVp_perc`). Avoid spaces and special characters. |
| `Location` | Location of the dataset: a local `*.jld2` file (absolute path, or relative to the directory in which you run julia), or a url starting with `http` from which the file is downloaded. The `.jld2` file must contain a GMG data structure, such as one saved with `save_GMG`. |
| `Type` | One of `Volume`, `Surface`, `Point`, `Topography` or `Screenshot` (see below). |
| `Active` | Optional, `true` or `false`. Inactive datasets are not loaded. If the column is missing or empty, `true` is assumed. |

The meaning of the different types is:
- `Volume`: 3D data, such as seismic tomography models. A cross-section is computed through every volume dataset.
- `Surface`: surfaces such as the Moho depth. For vertical profiles, the intersection of the surface with the profile is computed, which is shown as a line.
- `Point`: point data, such as earthquake hypocenters. All points within a band of width `section_width` around the profile are projected onto it.
- `Topography`: the topography of the region. It is shown above vertical profiles and in the map overview of the picker. The extent of the (first) topography dataset also determines the resolution of the profiles in the script below, so it should cover the region of interest.
- `Screenshot`: screenshots of figures from publications (see [Import screenshots](tutorial_Screenshot_To_Paraview.md)). These are not used in this tutorial.

This is the file `Datasets_Adria.txt` used in this tutorial (with the paths shortened to `DataSet/...`, a directory relative to where julia is started):
```
Name,Location,Type, [Active]
Giacomuzzi,DataSet/Giacomuzzi.jld2,Volume,true
Timko2023,DataSet/CaPaREA2023_Mantle_Timko.jld2,Volume,true
Friederich2025,DataSet/Friederich2025_Pwave_Alps.jld2,Volume,true
Koulakov,DataSet/Koulakov_Europe.jld2,Volume,true
ElSharkawy2020,DataSet/MeRE2020_El-Sharkawy.jld2,Volume,true
MIT08_Li,DataSet/MIT08_Pwave.jld2,Volume,true
Piromallo2003,DataSet/Piromallo2003.jld2,Volume,true
REVEAL,DataSet/REVEAL_rel_avg_AK135.jld2,Volume,true
SAVANI,DataSet/SAVANI_Auer2014.jld2,Volume,true
UU07_Amaru,DataSet/UU07_Pwave.jld2,Volume,true
Zhu2015,DataSet/Zhu2015.jld2,Volume,true
ETOPO1,DataSet/etopo1.jld2,Topography,true
Grad2009,DataSet/Grad2009_EU_Amr.jld2,Surface,true
ISC,DataSet/isc24_MedS.jld2,Point,true
```
It contains eleven tomographic models, the ETOPO1 topography, the Moho depth of Grad et al. (2009) and earthquake locations from the ISC catalogue. To exclude a dataset temporarily (for example, a large one while you test your profiles), set its last column to `false`.

!!! warning
    Since the file is comma-separated, the paths and urls must not contain commas.

The dataset file is read with `load_dataset_file`, and `load_GMG` then loads all active datasets and sorts them by type:
```julia
julia> datasets = load_dataset_file("Datasets_Adria.txt");
julia> data = load_GMG(datasets);
julia> keys(data)
(:Volume, :Surface, :Point, :Screenshot, :Topography)
julia> keys(data.Volume)
(:Giacomuzzi, :Timko2023, :Friederich2025, :Koulakov, :ElSharkawy2020, :MIT08_Li, :Piromallo2003, :REVEAL, :SAVANI, :UU07_Amaru, :Zhu2015)
```
Every entry of `data` is a `NamedTuple` with the datasets of that type, named as in the `Name` column.

## 3. Create the profiles
The following script projects all datasets onto all profiles and saves every profile as a `.pgmg` file. Besides GeophysicalModelGenerator, it needs the [JLD2](https://github.com/JuliaIO/JLD2.jl) package (`] add JLD2`).

```julia
using GeophysicalModelGenerator, JLD2

file_profiles = "profiles_mixed.txt"
file_datasets = "Datasets_Adria.txt"

# read the profiles and load all active datasets
profiles = read_picked_profiles(file_profiles)
datasets = load_dataset_file(file_datasets)
data     = load_GMG(datasets)
```

Next, we determine the resolution of the profiles. We take the lon/lat extent of the topography and the finest grid spacing of all volume datasets, so that no information of the highest-resolution model is lost:
```julia
lon_range = extrema(data.Topography[1].lon.val)
lat_range = extrema(data.Topography[1].lat.val)

res_lon = minimum(minimum(diff(vol.lon.val, dims = 1)) for vol in data.Volume)
res_lat = minimum(minimum(diff(vol.lat.val, dims = 2)) for vol in data.Volume)

nlon = round(Int, (lon_range[2] - lon_range[1]) / res_lon)
nlat = round(Int, (lat_range[2] - lat_range[1]) / res_lat)
```

Now we loop over the profiles, project the data onto every profile with `extract_ProfileData!`, and save the result:
```julia
for (i, profile) in enumerate(profiles)
    if profile.vertical
        # vertical profile: (points along the profile, points in depth)
        dims_vol = (nlon, 200)
        prefix   = "Profile_vertical"
    else
        # horizontal profile: (points in lon, points in lat)
        dims_vol = (nlon, nlat)
        prefix   = "Profile_horizontal"
    end

    extract_ProfileData!(profile, data.Volume, data.Surface, data.Point;
                         TopoData      = data.Topography,
                         DimsVolCross  = dims_vol,
                         Depth_extent  = (-400, 0),
                         DimsSurfCross = (nlon,),
                         section_width = 50km)

    # surfaces and topography are only intersected with vertical profiles,
    # so for horizontal profiles we add the full datasets for the map view
    if !profile.vertical
        profile.SurfData = data.Surface
        profile.TopoData = data.Topography
    end

    # save the profile for the picker
    jldsave("$(prefix)$(i).pgmg"; profile)
end
```

The options of `extract_ProfileData!` are:
- `DimsVolCross`: number of grid points of the cross-sections through the volume data (see above).
- `Depth_extent`: minimum and maximum depth (in km) of vertical profiles. Here, we use the uppermost 400 km.
- `DimsSurfCross`: number of points along the profile at which surfaces are intersected. The topography is sampled with 5 times as many points.
- `section_width`: width of the band around the profile from which point data (earthquakes) are projected onto the profile.

You can also call `extract_ProfileData!` with `nothing` (or an empty `NamedTuple()`) for any dataset type that you do not have, e.g. `extract_ProfileData!(profile, data.Volume, nothing, data.Point)`.

After running the script, you have one file per profile: `Profile_horizontal1.pgmg`, `Profile_horizontal2.pgmg` and `Profile_vertical3.pgmg` to `Profile_vertical7.pgmg`. A processed vertical profile contains:
```julia
julia> profiles[3]
Vertical ProfileData
  lon/lat    : (9.777223439354351, 42.936916884179524)-(15.465423125301268, 47.79656290565586)
    VolData  : (:x_profile, :Giacomuzzi_dVp_perc, :Giacomuzzi_dVs_perc, ..., :Piromallo2003_Vp, :Piromallo2003_dVp_perc, ...)
    SurfData : (:Grad2009,)
    TopoData : (:ETOPO1,)
        PointData: (:ISC,)
```
All volume datasets are combined into a single `GeoData` structure, in which the name of each field is prefixed with the name of its dataset. `x_profile` is the distance along the profile in km.

!!! tip
    Processing many large tomographic models can take a while and needs a fair amount of memory, as all datasets are loaded at the same time. Test your profiles first with only a few active datasets.

## 4. The `.pgmg` file format
A `.pgmg` file is a [JLD2](https://github.com/JuliaIO/JLD2.jl) file that contains a single `ProfileData` object; the name under which it is stored does not matter. The picker uses the following fields:

| Field | Content |
|:--|:--|
| `vertical` | `true` for a vertical profile, `false` for a horizontal slice |
| `start_lonlat`, `end_lonlat` | start and end point of a vertical profile |
| `depth` | depth of a horizontal slice |
| `VolData` | the cross-sections through all volume datasets, combined in one `GeoData` |
| `SurfData` | `NamedTuple` of `GeoData`, one per surface dataset |
| `PointData` | `NamedTuple` of `GeoData`, one per point dataset |
| `TopoData` | `NamedTuple` of `GeoData`, one per topography dataset |

You can load a profile back into julia with
```julia
julia> using JLD2
julia> profile = load("Profile_vertical3.pgmg", "profile");
```

## 5. Open the profiles in the AdriaArrayGeometryPicker
Install the picker as described in its [README](https://github.com/JuliaGeodynamics/AdriaArrayGeometryPicker.jl), and open a profile with
```julia
julia> using AdriaArrayGeometryPicker
julia> geometry_picker("Profile_vertical3.pgmg")
```
Alternatively, start `geometry_picker()` without arguments and use **Load Profile** in the menu. How to pick structures and save the picks is described in the documentation of the picker.
