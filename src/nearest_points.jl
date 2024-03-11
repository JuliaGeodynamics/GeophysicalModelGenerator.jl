# contains 1D,2D,3D nearest point routines
# mostly useful internally, but are exported as well (may be helpful in some cases)
using NearestNeighbors

export nearest_point_indices

# internal routine that does the computation

"""
    ind = nearest_point_indices(X::Array, X_pt::Vector)

Returns the index of the nearest point in (`X_pt`) to (`X`) and returns the index
"""
function nearest_point_indices(X::Array, X_pt::Vector)

    # use nearest neighbour to interpolate data
    coord   =   [X_pt';]
    kdtree  =   KDTree(coord; leafsize = 10);
    points  =   [vec(X)';];
    idx,_   =   nn(kdtree, points);
    
    # transform to correct shape
    ind     =   zeros(Int64,size(X))
    ind[:]  =   idx

    return ind
end

"""
    ind = nearest_point_indices(X::Array,Y::Array, X_pt::Vector, Y_pt::Vector)

Returns the index of the nearest point in (`X_pt`,`Y_pt`) to (`X`,`Y`) and returns the index
"""
function nearest_point_indices(X::Array,Y::Array, X_pt::Vector, Y_pt::Vector)

    # use nearest neighbour to interpolate data
    coord   =   [X_pt'; Y_pt'];
    kdtree  =   KDTree(coord; leafsize = 10);
    points  =   [vec(X)';vec(Y)'];
    idx,_   =   nn(kdtree, points);
    
    # transform to correct shape
    ind     =   zeros(Int64,size(X))
    ind[:]  =   idx

    return ind
end


"""
    ind = nearest_point_indices(X::Array,Y::Array,Z::Array, X_pt::Vector,Y_pt::Vector,Z_pt::Vector)

Returns the index of the nearest point in (`X_pt`,`Y_pt`,`Z_pt`) to (`X`,`Y`,`Z`) and returns the index
"""
function nearest_point_indices(X::Array,Y::Array,Z::Array, X_pt::Vector,Y_pt::Vector,Z_pt::Vector)

    # use nearest neighbour to interpolate data
    coord   =   [X_pt'; Y_pt'; Z_pt'];
    kdtree  =   KDTree(coord; leafsize = 10);
    points  =   [vec(X)';vec(Y)'; vec(Z)'];
    idx,_   =   nn(kdtree, points);
    
    # transform to correct shape
    ind     =   zeros(Int64,size(X))
    ind[:]  =   idx

    return ind
end
