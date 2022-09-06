# Tests data_import.jl

using Test
using GeophysicalModelGenerator


# test the import of a csv file with depth given as positive values and units in km

# test loaded structure
# NOTE: currently deactivated these tests -----------------------------
#LoadedData   =   ReadCSV_LatLon("TestData.csv", "positive");
#@test LoadedData.lat.name  ==  "lat"
#@test LoadedData.lat.unit  ==  "deg"
#@test LoadedData.lat.values  ==  [36.424, 36.194, 36.144, 35.486, 36.106, 35.767, 36.24, 36.301, 36.401, 36.321]

#@test LoadedData.lon.name  ==  "lon"
#@test LoadedData.lon.unit  ==  "deg"
#@test LoadedData.lon.values  ==  [68.74, 71.111, 70.757, 70.38, 71.051, 70.55, 68.95, 71.303, 71.567, 70.829]

#@test LoadedData.depth.name  ==  "depth"
#@test LoadedData.depth.unit  ==  "km"
#@test LoadedData.depth.values  ==  [-52.04, -60.51, -60.74, -61.0, -75.5, -81.52, -82.55, -84.36, -90.9, -110.41]

#@test LoadedData.values.varnames == ["local_magnitude", "absError(km)"]
#@test LoadedData.values.vals == [5.5   9.7;5.6  10.0;5.1  10.7;5.5  20.4;5.1   6.5;5.9  10.1;5.2   5.5;5.1   3.2;5.3   3.6;5.3   3.6]
# ---------------------------------------------------------------------


# test loading images (profiles & mapviews)

# Extract & save profile in GeoData format
filename            =   "test.png";             # fake png
Corner_LowerLeft    =   (18.0, 51.0, -590.0)
Corner_UpperRight   =   (9.0, 42.0,    0.0)
data_Image          =   Screenshot_To_GeoData(filename,Corner_LowerLeft, Corner_UpperRight)
@test data_Image.lon[1000] ≈ 17.592964824120603
@test data_Image.lat[1000] ≈ 50.59296482412061
@test Value(data_Image.depth[1000])==-590km
@test Write_Paraview(data_Image, "Profile_1")==nothing

# Test in CartData
data_Image          =   Screenshot_To_GeoData(filename,Corner_LowerLeft, Corner_UpperRight, Cartesian=true)
@test Value(data_Image.x[22]) == 18.0km
@test Value(data_Image.y[22]) == 51.0km
@test Value(data_Image.z[22]) ≈ -125.15151515151516km


# Test in UTM zone [note that depth should be in m]
data_Image          =   Screenshot_To_GeoData(filename,Corner_LowerLeft, Corner_UpperRight, UTM=true, UTMzone=33, isnorth=true)
@test data_Image.EW.val[22] == 18.0
@test data_Image.NS.val[22] == 51.0
@test Value(data_Image.depth[22]) ≈ -125.15151515151516m

# Mapview (distorted) in GeoData format
filename            =   "test.png";             # fake png
Corner_LowerLeft    =   (2.0,  40.0, -15.0)
Corner_UpperRight   =   (22.0, 51.0, -15.0)
Corner_LowerRight   =   (20.0, 40.0, -15.0)
Corner_UpperLeft    =   (0.0,  51.0, -15.0)
data_Image          =   Screenshot_To_GeoData(filename,Corner_LowerLeft, Corner_UpperRight, Corner_LowerRight=Corner_LowerRight, Corner_UpperLeft=Corner_UpperLeft)
@test data_Image.lon[1000]  ≈  2.814070351758794
@test data_Image.lat[1000]  ≈ 40.00000000000001
@test Value(data_Image.depth[1000])==-15km
@test Write_Paraview(data_Image, "MapView_1") == nothing

# MapView in CartData
data_Image          =   Screenshot_To_CartData(filename,Corner_LowerLeft, Corner_UpperRight, Corner_LowerRight=Corner_LowerRight, Corner_UpperLeft=Corner_UpperLeft)
@test Value(data_Image.x[22]) ≈ 0.42424242424242425km
@test Value(data_Image.y[22]) ≈ 48.666666666666664km
@test Value(data_Image.z[22]) ≈ -15km

# MapView in UTMData
data_Image          =   Screenshot_To_UTMData(filename,Corner_LowerLeft, Corner_UpperRight, Corner_LowerRight=Corner_LowerRight, Corner_UpperLeft=Corner_UpperLeft, UTMzone=33, isnorth=true)
@test data_Image.EW.val[22] ≈ 0.42424242424242425
@test data_Image.NS.val[22] ≈ 48.666666666666664
@test Value(data_Image.depth[22]) ≈ -15.0m
