# Tests data_import.jl

using Test
using GeophysicalModelGenerator

# test the import of a csv file with depth given as positive values and units in km

# test creating structures
LoadedData   =   ReadCSV_LatLon("TestData.csv", "positive");
@test LoadedData.lat.name  ==  "lat"
@test LoadedData.lat.unit  ==  "deg"
@test LoadedData.lat.values  ==  [68.740 71.111 70.757 70.380 71.051 70.550 68.950 71.303 71.567 70.829]

@test LoadedData.lon.name  ==  1000km
@test LoadedData.lon.unit  ==  1000km
@test LoadedData.lon.values  ==  1000km

@test LoadedData.depth.name  ==  1000km
@test LoadedData.depth.unit  ==  1000km
@test LoadedData.depth.values  ==  1000km

@test CharUnits_GEO.Pa      ==  10000000Pa
@test CharUnits_GEO.Mass    ==  1.0e37kg
@test CharUnits_GEO.Time    ==  1.0e12s
@test CharUnits_GEO.Length  ==  1000000m


CharUnits_SI   =   SI_units();
@test CharUnits_SI.length   ==  1000m
@test CharUnits_SI.Pa       ==  10Pa

CharUnits_NO   =   NO_units();
@test CharUnits_NO.length   ==  1
@test CharUnits_NO.Pa       ==  1
