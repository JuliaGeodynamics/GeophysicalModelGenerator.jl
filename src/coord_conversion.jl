# This is coord_conversion.jl
#
# This file contains functions to convert from geographical coordinates to Cartesian coordinates and vice versa.
#
# Author: Marcel Thielmann, 05/2021

function ConvertGeo2Cart(inputdata::GeoData)

    rearth = 3670e3

    lon = inputdata.lon.values
    lat = inputdata.lat.values
    depth = inputdata.depth.values

    R = rearth+inputdata.depth.values
    
    # compute Cartesian coordinates and assign them to a ValueList variable
    X = ValueList("x","m",R.*cosd(inputdata.lon.values).*cosd(inputdata.lat.values))
    Y = ValueList("y","m",R.*sind(inputdata.lon.values).*cosd(inputdata.lat.values))
    Z = ValueList("z","m",R.*sind(inputdata.lat.values))

    # assign all data to the respective struct
    convdata = CartData(X,Y,Z,inputdata.values) 

return convdata