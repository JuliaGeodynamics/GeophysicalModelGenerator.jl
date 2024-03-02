# This contains a number of routines that deal with surfaces
export remove_NaN_Surface!, drape_on_topo, is_surface
import Base: +,-

"""
    issurf = is_surface(surf::AbstractGeneralGrid)

Returns true if `surf` is a horizontal 3D surface.
"""
function is_surface(surf::AbstractGeneralGrid)
    if size(surf)[3] == 1
        issurf = true
    else
        issurf = false
    end
    return issurf
end


function  +(a::_T, b::_T) where _T<:AbstractGeneralGrid 
    @assert size(a) == size(b)
    return _addSurfaces(a,b)
end

function  -(a::_T, b::_T) where _T<:AbstractGeneralGrid 
    @assert size(a) == size(b)
    return _subtractSurfaces(a,b)
end

# Internal routines
_addSurfaces(a::_T, b::_T) where _T<:GeoData        = GeoData(a.lon.val, a.lat.val, a.depth.val + b.depth.val, merge(a.fields,b.fields))
_addSurfaces(a::_T, b::_T) where _T<:UTMData        = UTMData(a.EW.val, a.NS.val, a.depth.val + b.depth.val, merge(a.fields,b.fields))
_addSurfaces(a::_T, b::_T) where _T<:CartData       = CartData(a.x.val, a.y.val, a.z.val + b.z.val, merge(a.fields,b.fields))
_addSurfaces(a::_T, b::_T) where _T<:ParaviewData   = ParaviewData(a.x.val, a.y.val, a.z.val + b.z.val, merge(a.fields,b.fields))

_subtractSurfaces(a::_T, b::_T) where _T<:GeoData        = GeoData(a.lon.val, a.lat.val, a.depth.val - b.depth.val, merge(a.fields,b.fields))
_subtractSurfaces(a::_T, b::_T) where _T<:UTMData        = UTMData(a.EW.val, a.NS.val, a.depth.val - b.depth.val, merge(a.fields,b.fields))
_subtractSurfaces(a::_T, b::_T) where _T<:CartData       = CartData(a.x.val, a.y.val, a.z.val - b.z.val, merge(a.fields,b.fields))
_subtractSurfaces(a::_T, b::_T) where _T<:ParaviewData   = ParaviewData(a.x.val, a.y.val, a.z.val - b.z.val, merge(a.fields,b.fields))

"""
    remove_NaN_Surface!(Z::Array,X::Array,Y::Array)

Removes NaN's from a grid `Z` by taking the closest points as specified by `X` and `Y`.
"""
function remove_NaN_Surface!(Z,X,Y)
    @assert size(Z) == size(X) == size(Y)

    # use nearest neighbour to interpolate data
    id      = findall(isnan.(Z) .== false)
    id_NaN  = findall(isnan.(Z))

    coord   =   [X[id]'; Y[id]'];
    kdtree  =   KDTree(coord; leafsize = 10);

    points    = [X[id_NaN]'; Y[id_NaN]'];
    idx,dist  = nn(kdtree, points);

    Z[id_NaN] = Z[id[idx]]

    return nothing
end


"""
    Topo = drape_on_topo(Topo::GeoData, Data::GeoData)

This drapes fields of a data set `Data` on the topography `Topo`

"""
function drape_on_topo(Topo::GeoData, Data::GeoData)
    @assert is_surface(Topo)
    @assert is_surface(Data)
   
    Lon,Lat,_   =   LonLatDepthGrid( Topo.lon.val[:,1,1], Topo.lat.val[1,:,1],Topo.depth.val[1,1,:]);

    # use nearest neighbour to interpolate data
    idx         =   nearest_point_indices(Lon,Lat, vec(Data.lon.val), vec(Data.lat.val) ); 

    idx_out     =   findall(  (Lon .<  minimum(Data.lon.val)) .| (Lon .>  maximum(Data.lon.val)) .|
                              (Lat .<  minimum(Data.lat.val)) .| (Lat .>  maximum(Data.lat.val)) )

    fields_new  = Topo.fields;
    field_names = keys(Data.fields);

    for i = 1:length(Data.fields)

        if typeof(Data.fields[i]) <: Tuple

            # vector or anything that contains more than 1 field
            data_tuple = Data.fields[i]      # we have a tuple (likely a vector field), so we have to loop

            data_array = zeros(typeof(data_tuple[1][1]),size(Topo.lon.val,1),size(Topo.lon.val,2),size(Topo.lon.val,3),length(Data.fields[i]));     # create a 3D array that holds the 2D interpolated values
            unit_array = zeros(size(data_array));

            for j=1:length(data_tuple)
                data_field           =   data_tuple[j];
                tmp                  =   data_array[:,:,:,1];
                tmp                  =   data_field[idx]
                data_array[:,:,:,j]  =   tmp
            end

            data_new    = tuple([data_array[:,:,:,c] for c in 1:size(data_array,4)]...)       # transform 4D matrix to tuple

            # remove points outside domain
            for j=1:length(data_tuple)
                tmp           =   data_new[j];
                tmp[idx_out] .= NaN
                data_array[:,:,:,j]  =   tmp
            end
            data_new    = tuple([data_array[:,:,:,c] for c in 1:size(data_array,4)]...)       # transform 4D matrix to tuple

        else

            # scalar field
            data_new        =   zeros(typeof(Data.fields[i][1]), size(Topo.lon.val,1),size(Topo.lon.val,2),size(Topo.lon.val,3));
            data_new        =   Data.fields[i][idx]                                 # interpolate data field

        end

        # replace the one
        new_field   =   NamedTuple{(field_names[i],)}((data_new,))                          # Create a tuple with same name
        fields_new  =   merge(fields_new, new_field);                                       # replace the field in fields_new

    end

    Topo_new        =   GeoData(Topo.lon.val,Topo.lat.val,Topo.depth.val, fields_new)

    return Topo_new
end


"""
    drape_on_topo(Topo::CartData, Data::CartData)

Drapes Cartesian Data on topography
"""
function drape_on_topo(Topo::CartData, Data::CartData)
    @assert is_surface(Topo)
    @assert is_surface(Data)
    
    Topo_lonlat = GeoData(ustrip.(Topo.x.val),ustrip.(Topo.y.val), ustrip.(Topo.z.val), Topo.fields )
    Data_lonlat = GeoData(ustrip.(Data.x.val),ustrip.(Data.y.val), ustrip.(Data.z.val), Data.fields )

    Topo_new_lonlat = drape_on_topo(Topo_lonlat, Data_lonlat)

    Topo_new = CartData(Topo_new_lonlat.lon.val, Topo_new_lonlat.lat.val, Topo_new_lonlat.depth.val, Topo_new_lonlat.fields)

    return Topo_new
end


"""
    surf_new = project_points_to_surface(X::Array, Y::Array, lon_pt::Vector, lat_pt::Vector, depth_pt::Vector)

This changes the depth value of the surface `surf` to the `depth` value of the closest-by-points in (`lon_pt`,`lat_pt`, `depth_pt`) 

"""
function project_points_to_surface(surf::GeoData, lon_pt::Vector, lat_pt::Vector, depth_pt::Vector)

end

