using NearestNeighbors

export PointData2NearestGrid, CountMap


"""
    Grid_counts = PointData2NearestGrid(Point::CartData, Grid::CartData; radius_factor=1)

Uses nearest neighbour interpolation to count how many points (given by `Point`) are in the vicinity of a 3D `Grid`. 
The search radius is `R=radius_factor*(Δx² + Δy² + Δz²)^(1/3)`

`Point` should have 1D coordinate vectors

`Grid_counts` is `Grid` but with an additional field `Count` that has the number of hits
"""
function PointData2NearestGrid(Point::CartData, Grid::CartData; radius_factor=1)

    @assert length(size(Point.x)) == 1

    # call routine
    Count = PointData2NearestGrid(NumValue(Point.x),NumValue(Point.y), NumValue(Point.z), NumValue(Grid.x),NumValue(Grid.y),NumValue(Grid.z); radius_factor=radius_factor)

    # return CartGrid with added field
    return  AddField(Grid,"Count",Count);
end


"""
    Grid_counts = PointData2NearestGrid(pt_x,pt_y,pt_z, Grid::CartData; radius_factor=1)

Uses nearest neighbour interpolation to count how many points (given by `pt_x`,`pt_y`,`pt_z` coordinate vectors) are in the 
vicinity of 3D `CartGrid` specified by `Grid`. The search radius is `R=radius_factor*(Δx² + Δy² + Δz²)^(1/3)`

`Grid_counts` is `Grid` but with an additional field `Count` that has the number of hits
"""
function PointData2NearestGrid(pt_x,pt_y,pt_z, Grid::CartData; radius_factor=1)

    # call routine
    Count = PointData2NearestGrid(pt_x,pt_y,pt_z, NumValue(Grid.x),NumValue(Grid.y),NumValue(Grid.z); radius_factor=radius_factor)

    # return CartGrid with added field
    return  AddField(Grid,"Count",Count);
end


"""
    Grid_counts = PointData2NearestGrid(Point::GeoData, Grid::GeoData; radius_factor=1)

Uses nearest neighbour interpolation to count how many points (given by `Point`) are in the vicinity of a 3D `Grid`. 
The search radius is `R=radius_factor*(Δx² + Δy² + Δz²)^(1/3)`

`Point` should have 1D coordinate vectors

`Grid_counts` is `Grid` but with an additional field `Count` that has the number of hits
"""
function PointData2NearestGrid(Point::GeoData, Grid::GeoData; radius_factor=1)

    @assert length(size(Point.lon)) == 1

    # call routine
    Count = PointData2NearestGrid(NumValue(Point.lon),NumValue(Point.lat), NumValue(Point.depth), NumValue(Grid.lon),NumValue(Grid.lat),NumValue(Grid.depth); radius_factor=radius_factor)

    # return CartGrid with added field
    return  AddField(Grid,"Count",Count);
end


"""
    Grid_counts = PointData2NearestGrid(pt_x,pt_y,pt_z, Grid::GeoData; radius_factor=1)

Uses nearest neighbour interpolation to count how many points (given by `pt_x`,`pt_y`,`pt_z` coordinate vectors) are in the 
vicinity of 3D `GeoData` specified by `Grid`. The search radius is `R=radius_factor*(Δx² + Δy² + Δz²)^(1/3)`

`Grid_counts` is `Grid` but with an additional field `Count` that has the number of hits
"""
function PointData2NearestGrid(pt_x,pt_y,pt_z, Grid::GeoData; radius_factor=1)

    # call routine
    Count = PointData2NearestGrid(pt_x,pt_y,pt_z, NumValue(Grid.lon),NumValue(Grid.lat),NumValue(Grid.depth); radius_factor=radius_factor)

    # return CartGrid with added field
    return  AddField(Grid,"Count",Count);
end

"""
    count = PointData2NearestGrid(pt_x,pt_y,pt_z, X,Y,Z; radius_factor=1)

This uses nearest neighbour interpolation to count how many points (given by `pt_x`,`pt_y`,`pt_z` coordinate vectors) are in the 
vicinity of 3D grid point specified by `X`,`Y`,`Z` 3D coordinate arrays, with regular spacing `(Δx,Δy,Δz)`.
The search radius is `R=radius_factor*(Δx² + Δy² + Δz²)^(1/3)`

"""
function PointData2NearestGrid(pt_x,pt_y,pt_z, X,Y,Z; radius_factor=1)
    
    data = zeros(3,length(pt_x));
    data[1,:],data[2,:],data[3,:] = pt_x[:], pt_y[:], pt_z[:]
    tree = BallTree(data)                               # Generate tree with EQ data
    
    # Grid spacing
    Δx,Δy,Δz = X[2,1,1]-X[1], Y[1,2,1]-Y[1], Z[1,1,2]-Z[1]

    points     = zeros(3,length(X));
    points[1,:],points[2,:],points[3,:] = X[:], Y[:], Z[:]
 
    radius  =   radius_factor*(Δx^2 + Δy^2 + Δz^2)^(1/3)  # search radius
    idxs    =   inrange(tree, points, radius)           # find points (NearestNeighbors package)
    le      =   length.(idxs)                           # number of points

    count = zeros(Int64,size(X));
    for i in eachindex(X)
        count[i] = le[i]
    end

    return count
end

"""
    DatasetCountMap = CountMap(DataSet::GeoData,field::String,stepslon::Int64,stepslat::Int64)

Takes a 2D GeoData struct and counts entries of `field` per predefined control area. `field` should only consist of 1.0s and 0.0s. The control area is defined by `steplon` and `steplat`.
`steplon` is the number of control areas in longitude direction and `steplat` the number of control areas in latitude direction.
The counts per control area are normalized by the highest count.

```julia
julia> Data_Faults         = GeoData(Lon3D,Lat3D,Faults,(Faults=Faults,))
GeoData 
    size      : (375, 208, 1)
    lon       ϵ [ -9.932408319802885 : 34.93985125012068]
    lat       ϵ [ 35.086096468211394 : 59.919210145128545]
    depth     ϵ [ 0.0 : 1.0]
    fields    : (:Faults,)

julia> steplon  = 125
julia> steplat  = 70
julia> countmap = CountMap(Data_Faults,"Faults",steplon,steplat)

GeoData 
    size      : (124, 69, 1)
    lon       ϵ [ -9.751471789279 : 34.75891471959677]
    lat       ϵ [ 35.26604656731949 : 59.73926004602028]
    depth     ϵ [ 0.0 : 1.0]
    fields    : (:CountMap,)
```julia

"""
function CountMap(DataSet::GeoData,field::String,stepslon::Int64,stepslat::Int64)

    lon      = unique(DataSet.lon.val)
    lat      = unique(DataSet.lat.val)

    # create new lon/lat arrays which hold the boundaries of the control areas
    lonstep  = LinRange(lon[1],lon[end],stepslon)
    latstep  = LinRange(lat[1],lat[end],stepslat)
    dlon     = abs(lonstep[2]-lonstep[1])
    dlat     = abs(latstep[2]-latstep[1])
    loncen   = lonstep[1]+dlon/2:dlon:lonstep[end]-dlon/2
    latcen   = latstep[1]+dlat/2:dlat:latstep[end]-dlat/2
    countmap = zeros(length(loncen),length(latcen))

    expr =   Meta.parse(field)
    if !haskey(DataSet.fields,expr[1])
        error("The GeoData set does not have the field: $(expr[1])")
    end

    # count the ones in every control area
    for i in eachindex(loncen)
        for j in eachindex(latcen)
            indi          = findall((lon .>= lonstep[i]) .& (lon .<= lonstep[i+1]))
            indj          = findall((lat .>= latstep[j]) .& (lat .<= latstep[j+1]))
            dataint       = DataSet.fields[expr[1]][indi,indj,1]
            count         = sum(dataint)
            countmap[i,j] = count
        end
    end

    # normalize count in every control area
    maxcount = maximum(countmap)
    countmap = countmap ./ maxcount

    # create new GeoData
    Lon3D,Lat3D, Data = LonLatDepthGrid(loncen,latcen,0);
    Data[:,:,1]       .= countmap
    DatasetCountMap   = GeoData(Lon3D,Lat3D,Data,(CountMap=Data,))

    return DatasetCountMap
end
