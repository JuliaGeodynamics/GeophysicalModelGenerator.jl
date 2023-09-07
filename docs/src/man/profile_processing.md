# Profile processing

We have number of routines that makes it easier to read various datasets into GMG & create profiles through that data
The routines provided here have the following functionality:
- Read datasets (remote or local) that contains volumetric, surface, point, topograpy or screenshot data
- Define a profile (horizontal, vertical) with space for (projected) data
- Project earthquake (point) data onto the profile or intersect surfaces with a vertical profile (e.g., Moho data)

```@docs
GeophysicalModelGenerator.load_GMG
GeophysicalModelGenerator.save_GMG
GeophysicalModelGenerator.ProfileData
GeophysicalModelGenerator.ExtractProfileData
GeophysicalModelGenerator.CreateProfileData
GeophysicalModelGenerator.GMG_Dataset
GeophysicalModelGenerator.Load_Dataset_file
GeophysicalModelGenerator.combine_VolData
GeophysicalModelGenerator.ExtractProfileData!
GeophysicalModelGenerator.ReadPickedProfiles
```