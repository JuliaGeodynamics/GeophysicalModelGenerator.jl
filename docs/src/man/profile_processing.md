# Profile processing

We have number of routines that makes it easier to read various datasets into `GMG` & create profiles through that data
The routines provided here have the following functionality:
- Read datasets (remote or local) that contains volumetric, surface, point, topograpy or screenshot data
- Define a profile (horizontal, vertical) with space for (projected) data
- Project earthquake (point) data onto the profile or intersect surfaces with a vertical profile (e.g., Moho data)

```@docs
load_GMG
save_GMG
cross_section
ProfileData
extract_ProfileData
create_ProfileData
GMG_Dataset
load_dataset_file
combine_vol_data
extract_ProfileData!
read_picked_profiles
```