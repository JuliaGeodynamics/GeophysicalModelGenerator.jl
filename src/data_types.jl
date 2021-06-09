# This is data_types.jl
# contains type definitions to be used in GeophysicalModelGenerator

import Base: show

export GeoData, CartData, LonLatDepthGrid

# data structure for a list of values - TO BE REMOVED
mutable struct ValueList
    name::String
    unit::String
    values::Vector{Float64}
end

""" 
    GeoData(lon::Any, lat:Any, depth::GeoUnit, fields::NamedTuple)
    
    Data structure that holds one or several fields with longitude, latitude and depth information.

    `depth` can have units of meter, kilometer or be unitless; it will be converted to km.
    `fields` should ideally be a NamedTuple which allows you to specify the names of each of the fields. 
        If case you only pass one array we will convert it to a NamedTuple with default name
        Note that this is added as `(DataFieldName=Data,)` (don't forget the comma at the end)

    `lon`,`lat`,`depth` should all have the same size as each of the `fields`


# Example     
```julia-repl
julia> Lat         =   1.0:10.0;
julia> Lon         =   11.0:20.0;
julia> Depth       =   (-20:-11)*km;
julia> Data        =   zeros(size(Lon));
julia> Data_set    =   GeoData(Lat,Lon,Depth,(DataFieldName=Data,))   
GeoData 
  size  : (10,)
  lon   ϵ [ 1.0 - 10.0]
  lat   ϵ [ 11.0 - 20.0]
  depth ϵ [ -20 km - -11 km]
  fields: (:DataFieldName,)
```
"""
struct GeoData
    lon     ::  GeoUnit
    lat     ::  GeoUnit 
    depth   ::  GeoUnit
    fields  ::  NamedTuple 
    
    # Ensure that the data is of the correct format
    function GeoData(lon,lat,depth,fields)
        
        # check depth & convert it to units of km in case no units are given or it has different length units
        if unit.(depth)[1]==NoUnits 
            depth = depth*km                # in case depth has no dimensions
        end
        depth = uconvert.(km,depth)         # convert to km
        depth = GeoUnit(depth,km)           # convert to GeoUnit structure with units of km

        # fields should be a NamedTuple. In case we simply provide an array, lets transfer it accordingly
        if !(typeof(fields)<: NamedTuple)
            if (typeof(fields)<: Tuple)
                if length(fields)==1
                    fields = (DataSet1=first(fields),)  # The field is a tuple; create a NamedTuple from it
                else
                    error("Please employ a NamedTuple as input, rather than  a Tuple")  # out of luck
                end
            else
                fields = (DataSet1=fields,)
            end
        end

        if !(size(lon)==size(lat)==size(depth)==size(fields[1]))    
            error("The size of Lon/Lat/Depth and the Fields should all be the same!")
        end

        return new(lon,lat,depth,fields)
     end

end

# Print an overview of the Geodata struct:
function Base.show(io::IO, d::GeoData)
    println(io,"GeoData ")
    println(io,"  size  : $(size(d.lon))")
    println(io,"  lon   ϵ [ $(minimum(d.lon.val)) - $(maximum(d.lon.val))]")
    println(io,"  lat   ϵ [ $(minimum(d.lat.val)) - $(maximum(d.lat.val))]")
    println(io,"  depth ϵ [ $(minimum(d.depth.val)) - $(maximum(d.depth.val))]")
    println(io,"  fields: $(keys(d.fields))")
end

"""
    CartData(x::GeoUnit, y::GeoUnit, z::GeoUnit, values::NamedTuple)

Cartesian data in x/y/z coordinates to be used with Paraview
This is usually generated automatically 
"""
mutable struct CartData
    x       ::  GeoUnit
    y       ::  GeoUnit
    z       ::  GeoUnit
    values  ::  NamedTuple
end

# conversion function from GeoData -> CartData
function Base.convert(::Type{CartData}, d::GeoData)  
  
    R   =   Array(ustrip(d.depth.val)) .+ 6371.0;
    lon =   Array(ustrip(d.lon.val));
    lat =   Array(ustrip(d.lat.val));
    
    X = R .* cosd.( lon ) .* cosd.( lat );
    Y = R .* sind.( lon ) .* cosd.( lat );
    Z = R .* sind.( lat );

    return CartData(GeoUnit(X,km),GeoUnit(Y,km),GeoUnit(Z,km),d.fields)
end


"""
    LonLatDepthGrid(Lon::Any, Lat::Any, Depth:Any)

Creates 3D arrays of `Lon`, `Lat`, `Depth` from 1D vectors or numbers

# Example 1: Create 3D grid
```julia-repl
julia> Lon,Lat,Depth =  LonLatDepthGrid(10:20,30:40,(-10:-1)km);
julia> size(Lon)
(11, 11, 10)
```

# Example 2: Create 2D lon/lat grid @ a given depth
```julia-repl
julia> Lon,Lat,Depth =  LonLatDepthGrid(10:20,30:40,-50km);
julia> size(Lon)
(11, 11)
```

# Example 3: Create 2D lon/depth grid @ a given lat
```julia-repl
julia> Lon,Lat,Depth =  LonLatDepthGrid(10:20,30,(-10:-1)km);
julia> size(Lon)
(11, 11)
```
# Example 4: Create 1D vertical line @ a given lon/lat point
```julia-repl
julia> Lon,Lat,Depth =  LonLatDepthGrid(10,30,(-10:-1)km);
julia> size(Lon)
(10, )
```

"""
function LonLatDepthGrid(Lon::Any, Lat::Any, Depth::Any)

    nLon    = length(Lon)
    nLat    = length(Lat)
    nDepth  = length(Depth)

    if nLon==nLat==nDepth==1
        error("Cannot use this routine for a 3D point (no need to create a grid in that case")
    end 
    if maximum([length(size(Lon)), length(size(Lat)), length(size(Depth))])>1
        error("You can only give 1D vectors or numbers as input")
    end

    Lon3D   =   zeros(nLon,nLat,nDepth);
    Lat3D   =   zeros(nLon,nLat,nDepth);
    Depth3D =   zeros(nLon,nLat,nDepth);

    for i=1:nLon
        for j=1:nLat
            for k=1:nDepth
                Lon3D[i,j,k]    =   ustrip(Lon[i]);
                Lat3D[i,j,k]    =   ustrip(Lat[j]);
                Depth3D[i,j,k]  =   ustrip(Depth[k]);
            end
        end
    end

    # drop dimensions that are of length one. Ofcourse, this implies that we cannot have
    Lon3D   = dropdims(Lon3D,   dims = (findall(size(Lon3D  ) .== 1)...,))
    Lat3D   = dropdims(Lat3D,   dims = (findall(size(Lat3D  ) .== 1)...,))
    Depth3D = dropdims(Depth3D, dims = (findall(size(Depth3D) .== 1)...,))
    

    # Add dimensions back
    Lon3D   = Lon3D*unit(Lon[1])
    Lat3D   = Lat3D*unit(Lat[1])
    Depth3D = Depth3D*unit(Depth[1])

    return Lon3D, Lat3D, Depth3D
end

