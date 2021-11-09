# This is coord_conversion.jl
#
# This file contains functions to convert from geographical coordinates to Cartesian coordinates and vice versa.
#
# Author: Marcel Thielmann, 05/2021

function ConvertGeo2Cart(inputdata::GeoData)

     if cmp(inputdata.depth.unit,"km")
        rearth = 3671
     elseif cmp(inputdata.depth.unit,"m")
        rearth = 3671e3
     else
        warn("No proper depth unit is given (km or m), using m")
        rearth = 3671e3
     end

    lon = inputdata.lon.values
    lat = inputdata.lat.values
    depth = inputdata.depth.values

    R = rearth+inputdata.depth.values # radius from the center of the earth
    
    # compute Cartesian coordinates and assign them to a ValueList variable
    X = ValueList("x",inputdata.depth.unit,R.*cosd(inputdata.lon.values).*cosd(inputdata.lat.values))
    Y = ValueList("y",inputdata.depth.unit,R.*sind(inputdata.lon.values).*cosd(inputdata.lat.values))
    Z = ValueList("z",inputdata.depth.unit,R.*sind(inputdata.lat.values))

    # assign all data to the respective struct
    convdata = ParaviewData(X,Y,Z,inputdata.values) 

   return convdata

end