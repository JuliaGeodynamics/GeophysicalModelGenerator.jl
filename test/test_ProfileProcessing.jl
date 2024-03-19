using Test
using GeophysicalModelGenerator

pkg_dir = pkgdir(GeophysicalModelGenerator)
cd(joinpath(pkg_dir,"test"))

# Test profile processing dataset routines
data_Surf = GMG_Dataset("Mrozek_Moho_Grid_EU","Surface","https://seafile.rlp.net/f/483d9c7c808a4087ba9e/?dl=1", true)
@test data_Surf.DirName == "https://seafile.rlp.net/f/483d9c7c808a4087ba9e/?dl=1"
@test data_Surf.Type == "Surface"
@test data_Surf.active == true
@test data_Surf.Name == "Mrozek_Moho_Grid_EU"

# Specify a few more profiles
data_EQ     = GMG_Dataset("AlpArraySeis","Point","https://seafile.rlp.net/f/87d565882eda40689666/?dl=1", true)
data_SS     = GMG_Dataset("Handy_etal_SE_Profile1","Screenshot","https://seafile.rlp.net/f/5ffe580e765e4bd1bafe/?dl=1", true)

# Note: the volumetric datasets are chosen as they are smaller in size (less download)
data_Vol1   = GMG_Dataset("Hua2017","Volume","https://seafile.rlp.net/f/1fb68b74e5d742d39e62/?dl=1", true)
data_Vol2   = GMG_Dataset("Plomerova2022","Volume","https://seafile.rlp.net/f/abccb8d3302b4ef5af17/?dl=1", true)
#data_Vol1   = GMG_Dataset("Paffrath2021","Volume","https://seafile.rlp.net/f/5c8c851af6764b5db20d/?dl=1", true)
#data_Vol2   = GMG_Dataset("Zhao2016","Volume","https://seafile.rlp.net/f/e81a6d075f6746609973/?dl=1", true)

# Now load these datasets into NamedTuples

SurfData        =   load_GMG(data_Surf)
PointData       =   load_GMG(data_EQ)
ScreenshotData  =   load_GMG(data_SS)
VolData         =   load_GMG(data_Vol1)
VolData         =   merge(VolData, load_GMG(data_Vol2))

# Combine all Datasets into one file
Datasets        =   [data_Vol1,data_Vol2, data_Surf, data_EQ, data_SS]

# Some tests with the loaded datasets
@test SurfData.Mrozek_Moho_Grid_EU.fields.MohoDepth[100,100] ≈ -58.6889km
@test keys(VolData) == (:Hua2017, :Plomerova2022)

# read datasets from file
Datasets_temp = load_dataset_file("test_files/AlpineData.txt")
@test Datasets_temp[2].DirName == GMG_Dataset("INGV","Point","./Seismicity/CLASS/class_seis_alps.jld2", true).DirName

# Load data of all Datasets & split them in type of data
Data = load_GMG(Datasets)
@test keys(Data.Volume) == (:Hua2017, :Plomerova2022)

# Combine volumetric datasets into one
VolData_combined1 = combine_vol_data(Data.Volume)
@test keys(VolData_combined1.fields) == (:Hua2017_Vp, :Hua2017_dVp_perc, :Plomerova2022_Vp, :Plomerova2022_dVp)

VolData_combined2 = combine_vol_data(Data.Volume, dims=(50,51,52))
@test VolData_combined2.fields.Hua2017_Vp[1000] ≈ 10.6904

VolData_combined3 = combine_vol_data(Data.Volume, lon=(1,22), lat=(40,52), dims=(50,51,52))
@test isnan(VolData_combined3.fields.Hua2017_Vp[1000])

# Define horizontal & vertical profiles
prof1 = ProfileData(start_lonlat=(5,45), end_lonlat=(15,49))
prof2 = ProfileData(depth = -100)
prof3 = ProfileData(start_lonlat=(5,45), end_lonlat=(5,49))
prof4 = ProfileData(depth = -20)

# test internal routines to intersect profile with volumetric data:
GeophysicalModelGenerator.create_profile_volume!(prof1, VolData_combined1)
@test prof1.VolData.fields.Hua2017_Vp[30,40] ≈ 9.141520976523731

GeophysicalModelGenerator.create_profile_volume!(prof2, VolData_combined1)
@test prof2.VolData.fields.Hua2017_Vp[30,40] ≈ 8.177263544536272

GeophysicalModelGenerator.create_profile_volume!(prof1, VolData_combined1,  Depth_extent=(-300, -100))
@test extrema(prof1.VolData.depth.val) == (-300.0, -100.0)

# Intersect surface data:
GeophysicalModelGenerator.create_profile_surface!(prof1,Data.Surface)
@test prof1.SurfData[1].fields.MohoDepth[80] ≈ -37.58791461075397km

# ditto with EQ data:
GeophysicalModelGenerator.create_profile_point!(prof1,Data.Point, section_width=5km)
GeophysicalModelGenerator.create_profile_point!(prof4,Data.Point, section_width=10km)
@test  length(prof1.PointData[1].lon) == 13
@test  length(prof4.PointData[1].lon) == 445


# Test the main profile extraction routines:
extract_ProfileData!(prof1, VolData_combined1, Data.Surface, Data.Point)
extract_ProfileData!(prof2, VolData_combined1, Data.Surface, Data.Point)
extract_ProfileData!(prof3, VolData_combined1, Data.Surface, Data.Point)
extract_ProfileData!(prof4, VolData_combined1, Data.Surface, Data.Point)

extract_ProfileData!(prof1, VolData_combined2, Data.Surface, Data.Point)
extract_ProfileData!(prof2, VolData_combined2, Data.Surface, Data.Point)
extract_ProfileData!(prof3, VolData_combined2, Data.Surface, Data.Point)
extract_ProfileData!(prof4, VolData_combined2, Data.Surface, Data.Point)

extract_ProfileData!(prof1, VolData_combined3, Data.Surface, Data.Point)
extract_ProfileData!(prof2, VolData_combined3, Data.Surface, Data.Point)
extract_ProfileData!(prof3, VolData_combined3, Data.Surface, Data.Point)
extract_ProfileData!(prof4, VolData_combined3, Data.Surface, Data.Point)


# Test that it works if only EQ's are provided:
prof4 = ProfileData(depth = -20)
extract_ProfileData!(prof4, nothing, NamedTuple(), Data.Point)
@test isnothing(prof4.VolData)
@test isempty(prof4.SurfData)
@test length(prof4.PointData[1].depth) == 3280

@test prof1.SurfData[1].fields[1][80] ≈ -37.58791461075397km
@test isempty(prof2.SurfData)
@test isnan(prof3.SurfData[1].fields[1][80])
@test isempty(prof4.SurfData)

# Read profiles from file
profile_list = read_picked_profiles("test_files/PickedProfiles.txt")
@test profile_list[5].start_lonlat == ProfileData(start_lonlat=(9.40627872242647, 45.5128223429144), end_lonlat=(7.85480813419117, 47.8635353553922)).start_lonlat

# Try the convenience function
DepthVol=nothing
DimsVolCross=(100,100)
Depth_extent=nothing
DimsSurfCross=(100,)
section_width=50km

profile_backwards_compat = extract_ProfileData("test_files/PickedProfiles.txt",1,"test_files/AlpineData_remote.txt",DimsVolCross=DimsVolCross,DepthVol=Depth_extent,DimsSurfCross=DimsSurfCross,WidthPointProfile=section_width)

@test length(profile_backwards_compat.PointData[1].lon) == 440
