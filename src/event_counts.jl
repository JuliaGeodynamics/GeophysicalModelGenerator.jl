using NearestNeighbors

export PointData2NearestGrid


"""
    Grid_counts = PointData2NearestGrid(Point::CartData, Grid::CartData; radius_factor=1)

Uses nearest neigbour interpolation to count how many points (given by `Point`) are in the vicinity of a 3D `Grid`. 
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

Uses nearest neigbour interpolation to count how many points (given by `pt_x`,`pt_y`,`pt_z` coordinate vectors) are in the 
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

Uses nearest neigbour interpolation to count how many points (given by `Point`) are in the vicinity of a 3D `Grid`. 
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

Uses nearest neigbour interpolation to count how many points (given by `pt_x`,`pt_y`,`pt_z` coordinate vectors) are in the 
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

This uses nearest neigbour interpolation to count how many points (given by `pt_x`,`pt_y`,`pt_z` coordinate vectors) are in the 
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

