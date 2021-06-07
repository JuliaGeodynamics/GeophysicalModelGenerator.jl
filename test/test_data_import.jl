# Tests data_import.jl

using Test
using GeophysicalModelGenerator

# test the import of a csv file with depth given as positive values and units in km

# test loaded structure
LoadedData   =   ReadCSV_LatLon("TestData.csv", "positive");
@test LoadedData.lat.name  ==  "lat"
@test LoadedData.lat.unit  ==  "deg"
@test LoadedData.lat.values  ==  [36.424, 36.194, 36.144, 35.486, 36.106, 35.767, 36.24, 36.301, 36.401, 36.321]

@test LoadedData.lon.name  ==  "lon"
@test LoadedData.lon.unit  ==  "deg"
@test LoadedData.lon.values  ==  [68.74, 71.111, 70.757, 70.38, 71.051, 70.55, 68.95, 71.303, 71.567, 70.829]

@test LoadedData.depth.name  ==  "depth"
@test LoadedData.depth.unit  ==  "km"
@test LoadedData.depth.values  ==  [-52.04, -60.51, -60.74, -61.0, -75.5, -81.52, -82.55, -84.36, -90.9, -110.41]

@test LoadedData.values.varnames == ["local_magnitude", "absError(km)"]
@test LoadedData.values.vals == [5.5   9.7;5.6  10.0;5.1  10.7;5.5  20.4;5.1   6.5;5.9  10.1;5.2   5.5;5.1   3.2;5.3   3.6;5.3   3.6]