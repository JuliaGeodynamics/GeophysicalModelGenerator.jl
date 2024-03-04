# This contains a number of routines that deal with surfaces
export remove_NaN_Surface!, drape_on_topo, is_surface, fit_surface_to_points
export aboveSurface, belowSurface, interpolateDataOnSurface
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
    surf_new = fit_surface_to_points(surf::GeoData, lon_pt::Vector, lat_pt::Vector, depth_pt::Vector)

This fits the `depth` values of the surface `surf` to the `depth` value of the closest-by-points in (`lon_pt`,`lat_pt`, `depth_pt`) 

"""
function fit_surface_to_points(surf::GeoData, lon_pt::Vector, lat_pt::Vector, depth_pt::Vector)
    @assert is_surface(surf)
    
    idx = nearest_point_indices(NumValue(surf.lon),NumValue(surf.lat),  lon_pt, lat_pt);
    depth = NumValue(surf.depth)
    depth[idx] .= depth_pt[idx];

    surf_new = surf
    surf_new.depth .= depth
    return surf_new
end


"""
    surf_new = fit_surface_to_points(surf::CartData, lon_pt::Vector, lat_pt::Vector, depth_pt::Vector)

This fits the `depth` values of the surface `surf` to the `depth` value of the closest-by-points in (`lon_pt`,`lat_pt`, `depth_pt`) 

"""
function fit_surface_to_points(surf::CartData, X_pt::Vector, Y_pt::Vector, Z_pt::Vector)
    @assert is_surface(surf)
    
    idx = nearest_point_indices(NumValue(surf.x),NumValue(surf.y),  X_pt[:], Y_pt[:]);
    depth = NumValue(surf.z)
    depth = Z_pt[idx]

    surf_new = deepcopy(surf)
    surf_new.z.val .= depth
    return surf_new
end



"""
    aboveSurface(Data::GeoData, DataSurface::GeoData; above=true)

Returns a boolean array of size(Data.Lon), which is true for points that are above the surface DataSurface (or for points below if above=false).

This can be used, for example, to mask points above/below the Moho in a volumetric dataset or in a profile.

# Example
First we create a 3D data set and a 2D surface:
```julia
julia> Lon,Lat,Depth   =   LonLatDepthGrid(10:20,30:40,(-300:25:0)km);
julia> Data            =   Depth*2;
julia> Data_set3D      =   GeoData(Lon,Lat,Depth,(Depthdata=Data,LonData=Lon))
GeoData
  size  : (11, 11, 13)
  lon   ϵ [ 10.0 : 20.0]
  lat   ϵ [ 30.0 : 40.0]
  depth ϵ [ -300.0 km : 0.0 km]
  fields: (:Depthdata, :LonData)
julia> Lon,Lat,Depth   =   LonLatDepthGrid(10:20,30:40,-40km);
julia> Data_Moho       =   GeoData(Lon,Lat,Depth+Lon*km, (MohoDepth=Depth,))
  GeoData
    size  : (11, 11, 1)
    lon   ϵ [ 10.0 : 20.0]
    lat   ϵ [ 30.0 : 40.0]
    depth ϵ [ -30.0 km : -20.0 km]
    fields: (:MohoDepth,)
```
Next, we intersect the surface with the data set:
```julia
julia> Above       =   aboveSurface(Data_set3D, Data_Moho);
```
Now, `Above` is a boolean array that is true for points above the surface and false for points below and at the surface.

"""
function aboveSurface(Data::GeoData, DataSurface::GeoData; above=true)

    if size(DataSurface.lon)[3]!=1
        error("It seems that DataSurface is not a surface")
    end

    # Create interpolation object for surface
    Lon_vec     =  DataSurface.lon.val[:,1,1];
    Lat_vec     =  DataSurface.lat.val[1,:,1];
    interpol    =  linear_interpolation((Lon_vec, Lat_vec), ustrip.(DataSurface.depth.val[:,:,1]));            # create interpolation object

    DepthSurface = interpol.(Data.lon.val,Data.lat.val);
    DepthSurface = DepthSurface*unit(DataSurface.depth.val[1])

    if above
        Above    =   Data.depth.val .> DepthSurface;
    else
        Above    =   Data.depth.val .< DepthSurface;
    end

    return Above
end

"""
    Below = belowSurface(Data::GeoData, DataSurface::GeoData)

Determines if points within the 3D `Data` structure are below the GeoData surface `DataSurface`
"""
function belowSurface(Data::GeoData, DataSurface::GeoData)
    return aboveSurface(Data::GeoData, DataSurface::GeoData; above=false)
end

"""
    Above = aboveSurface(Data_Cart::ParaviewData, DataSurface_Cart::ParaviewData; above=true)

Determines if points within the 3D `Data_Cart` structure are above the Cartesian surface `DataSurface_Cart`
"""
function aboveSurface(Data_Cart::ParaviewData, DataSurface_Cart::ParaviewData; above=true)

    Data            =   GeoData(ustrip.(Data_Cart.x.val),       ustrip.(Data_Cart.y.val),        ustrip.(Data_Cart.z.val), Data_Cart.fields)
    DataSurface     =   GeoData(ustrip.(DataSurface_Cart.x.val),ustrip.(DataSurface_Cart.y.val), ustrip.(DataSurface_Cart.z.val), DataSurface_Cart.fields )

    return Above    =   aboveSurface(Data, DataSurface; above=above)
end

"""
    Above = aboveSurface(Data_Cart::CartData, DataSurface_Cart::CartData; above=true)

Determines if points within the 3D `Data_Cart` structure are above the Cartesian surface `DataSurface_Cart`
"""
function aboveSurface(Data_Cart::CartData, DataSurface_Cart::CartData; above=true)

    Data            =   GeoData(ustrip.(Data_Cart.x.val),       ustrip.(Data_Cart.y.val),        ustrip.(Data_Cart.z.val), Data_Cart.fields)
    DataSurface     =   GeoData(ustrip.(DataSurface_Cart.x.val),ustrip.(DataSurface_Cart.y.val), ustrip.(DataSurface_Cart.z.val), DataSurface_Cart.fields )

    return Above    =   aboveSurface(Data, DataSurface; above=above)
end

"""
    Above = aboveSurface(Grid::CartGrid, DataSurface_Cart::CartData; above=true)

Determines if points described by the `Grid` CartGrid structure are above the Cartesian surface `DataSurface_Cart`
"""
function aboveSurface(Grid::CartGrid, DataSurface_Cart::CartData; above=true)

    X,Y,Z = XYZGrid(Grid.coord1D...)
    Data = CartData(Grid,(Z=Z,))

    return aboveSurface(Data, DataSurface_Cart; above=above)
end


"""
    Below = belowSurface(Grid::CartGrid, DataSurface_Cart::CartData)

    Determines if points described by the `Grid` CartGrid structure are above the Cartesian surface `DataSurface_Cart`
"""
function belowSurface(Grid::CartGrid, DataSurface_Cart::CartData)
    return aboveSurface(Grid, DataSurface_Cart; above=false)
end


"""
    Below = belowSurface(Data_Cart::ParaviewData, DataSurface_Cart::ParaviewData)

Determines if points within the 3D Data_Cart structure are below the Cartesian surface DataSurface_Cart
"""
function belowSurface(Data_Cart::ParaviewData, DataSurface_Cart::ParaviewData)
    return aboveSurface(Data_Cart::ParaviewData, DataSurface_Cart::ParaviewData; above=false)
end

"""
    Below = belowSurface(Data_Cart::CartData, DataSurface_Cart::CartData)

Determines if points within the 3D Data_Cart structure are below the Cartesian surface DataSurface_Cart
"""
function belowSurface(Data_Cart::CartData, DataSurface_Cart::CartData)
    return aboveSurface(Data_Cart::CartData, DataSurface_Cart::CartData; above=false)
end

"""
    Surf_interp = interpolateDataOnSurface(V::ParaviewData, Surf::ParaviewData)

Interpolates a 3D data set `V` on a surface defined by `Surf`.
# Example
```julia
julia> Data
ParaviewData
  size  : (33, 33, 33)
  x     ϵ [ -3.0 : 3.0]
  y     ϵ [ -2.0 : 2.0]
  z     ϵ [ -2.0 : 0.0]
  fields: (:phase, :density, :visc_total, :visc_creep, :velocity, :pressure, :temperature, :dev_stress, :strain_rate, :j2_dev_stress, :j2_strain_rate, :plast_strain, :plast_dissip, :tot_displ, :yield, :moment_res, :cont_res)
julia> surf
ParaviewData
  size  : (96, 96, 1)
  x     ϵ [ -2.9671875 : 3.2671875]
  y     ϵ [ -1.9791666666666667 : 1.9791666666666667]
  z     ϵ [ -1.5353766679763794 : -0.69925457239151]
  fields: (:Depth,)
julia> Surf_interp = interpolateDataOnSurface(Data, surf)
  ParaviewData
    size  : (96, 96, 1)
    x     ϵ [ -2.9671875 : 3.2671875]
    y     ϵ [ -1.9791666666666667 : 1.9791666666666667]
    z     ϵ [ -1.5353766679763794 : -0.69925457239151]
    fields: (:phase, :density, :visc_total, :visc_creep, :velocity, :pressure, :temperature, :dev_stress, :strain_rate, :j2_dev_stress, :j2_strain_rate, :plast_strain, :plast_dissip, :tot_displ, :yield, :moment_res, :cont_res)
```
"""
function interpolateDataOnSurface(V::ParaviewData, Surf::ParaviewData)

    # Create GeoData structure:
    V_geo               =   GeoData(V.x.val, V.y.val, V.z.val, V.fields)
    V_geo.depth.val     =   ustrip(V_geo.depth.val);

    Surf_geo            =   GeoData(Surf.x.val, Surf.y.val, Surf.z.val, Surf.fields)
    Surf_geo.depth.val  =   ustrip(Surf_geo.depth.val);

    Surf_interp_geo     =   interpolateDataOnSurface(V_geo, Surf_geo)
    Surf_interp         =   ParaviewData(Surf_interp_geo.lon.val, Surf_interp_geo.lat.val, ustrip.(Surf_interp_geo.depth.val), Surf_interp_geo.fields)

    return Surf_interp

end

function interpolateDataOnSurface(V::CartData, Surf::CartData)

    # Create GeoData structure:
    V_geo               =   GeoData(V.x.val, V.y.val, V.z.val, V.fields)
    V_geo.depth.val     =   ustrip(V_geo.depth.val);

    Surf_geo            =   GeoData(Surf.x.val, Surf.y.val, Surf.z.val, Surf.fields)
    Surf_geo.depth.val  =   ustrip(Surf_geo.depth.val);

    Surf_interp_geo     =   interpolateDataOnSurface(V_geo, Surf_geo)
    Surf_interp         =   CartData(Surf_interp_geo.lon.val, Surf_interp_geo.lat.val, ustrip.(Surf_interp_geo.depth.val), Surf_interp_geo.fields)

    return Surf_interp

end

"""
    Surf_interp = interpolateDataOnSurface(V::GeoData, Surf::GeoData)

Interpolates a 3D data set `V` on a surface defined by `Surf`
"""
function interpolateDataOnSurface(V::GeoData, Surf::GeoData)

    Surf_interp = InterpolateDataFields(V, Surf.lon.val, Surf.lat.val, Surf.depth.val)

    return Surf_interp
end