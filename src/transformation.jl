# This provides various transformations (GeoData <=> Cartesian)
export GeoData_To_Cartesian, Cartesian_To_GeoData, InterpolateToRectilinear

"""
    data_Cart, refPoint = GeoData_To_Cartesian(data::GeoData; Flatten=true, lonlatDepth=empty, referencePoint="Center")

Projects a GeoData structure (`data`) to a local Cartesian structure around a reference point `refPoint`. 

If `Flatten=true`, the resulting grid will be flattened, which implies that the curvature of the Earth is ignored for the region (and the depth level is the same as in the GeoData struct). 
This is consistent with many 3D Cartesian (geodynamic) codes, but you ofcourse make some errors while doing so which increase for larger model domains.

There are two ways to define the reference point:

    - By specifying the keyword `referencePoint`. This can have several options such as "Center", "CenterBottom", "CenterTop", "FrontLeftTop", "BackRightTop"
    - By specifying the lonlatDepth vector which has [`longitude`, `latitude`, `depth`] values (with depth in km)
    
# Examples 
```julia
julia> Lon,Lat,Depth   =   LonLatDepthGrid(10:20,30:40,(-300:25:0)km);
julia> Data_set3D      =   GeoData(Lon,Lat,Depth,(Depthdata=Depth*2,LonData=Lon))
GeoData 
  size  : (11, 11, 13)
  lon   ϵ [ 10.0 : 20.0]
  lat   ϵ [ 30.0 : 40.0]
  depth ϵ [ -300.0 km : 0.0 km]
  fields: (:Depthdata, :LonData)
```
A case with no flattening is done with
```julia
julia> data_Cart, refPoint     = GeoData_To_Cartesian(Data_set3D, referencePoint="CenterBottom", Flatten=false)
  (CartData 
    size  : (11, 11, 13)
    lon   ϵ [ -481.81931349723936 km : 481.8193134972397 km]
    lat   ϵ [ -553.7761439713204 km : 564.9108783864752 km]
    depth ϵ [ -39.455981623586645 km : 300.00000000000074 km]
    fields: (:Depthdata, :LonData)
  , LLA(lat=35.0°, lon=15.0°, alt=-300000.0))
```
Whereas the default flattens the region:
```julia
julia> data_Cart, refPoint     = GeoData_To_Cartesian(Data_set3D, referencePoint="CenterBottom")
(CartData 
  size  : (11, 11, 13)
  lon   ϵ [ -481.8193134972393 km : 481.81931349723976 km]
  lat   ϵ [ -553.7761439713197 km : 564.9108783864759 km]
  depth ϵ [ -300.0 km : 0.0 km]
  fields: (:Depthdata, :LonData)
, LLA(lat=35.0°, lon=15.0°, alt=0.0))
```
Note the difference in `z`-coordinates between flattening and not.

The resulting Cartesian grid can be passed on to 2D/3D codes. Note, however, that this will NOT be a regular grid in `x/y` direction, but rather slightly distorted, depending on how large the region is. 
So in general you will have to interpolate all fields from `data_Cart` to a regular grid, before using them in such codes.
"""
function GeoData_To_Cartesian(data::GeoData; Flatten=true, lonlatDepth=empty, referencePoint="Center")

    # compute min/max/average lon/lat/z
    lon_data    = [minimum(ustrip.(data.lon.val))   maximum(ustrip.(data.lon.val)) 0]
    lat_data    = [minimum(ustrip.(data.lat.val))   maximum(ustrip.(data.lat.val)) 0]
    z_data      = [minimum(ustrip.(data.depth.val)) maximum(ustrip.(data.depth.val)) 0]
    lon_data[3] = (lon_data[1] + lon_data[2])/2
    lat_data[3] = (lat_data[1] + lat_data[2])/2
    z_data[3]   = (z_data[1]   +   z_data[2])/2

    if lonlatDepth == empty
        if      referencePoint=="Center"
            lon,lat,z = lon_data[3], lat_data[3], z_data[3];
        elseif  referencePoint=="CenterBottom"
            lon,lat,z = lon_data[3], lat_data[3], z_data[1];
        elseif  referencePoint=="CenterTop"
            lon,lat,z = lon_data[3], lat_data[3], z_data[2];
        elseif  referencePoint=="FrontLeftTop"
            lon,lat,z = lon_data[1], lat_data[1], z_data[2];
        elseif  referencePoint=="BackRightTop"
            lon,lat,z = lon_data[2], lat_data[2], z_data[2];
        else
            error("Unknown reference point: $referencePoint")
        end

    else 
        lon,lat,z = lonlatDepth[1], lonlatDepth[2], lonlatDepth[3];
    end
    if Flatten==true
        z = 0
    end
    refPoint    =   LLA(lat,lon,z*1e3)
   
    # transfer to LLA (Geodesy.jl), with depth in meters
    LLA_Data    =   LLA.(data.lat.val, data.lon.val, ustrip(data.depth.val)*1e3)
    trans       =   ENUfromLLA(refPoint, wgs84)
    XYZ         =   trans.(LLA_Data)/1e3               # perform transformation & transfer to km
   
    # convert to cartesian ENU frame
    X,Y,Z       =   zeros(size(XYZ)), zeros(size(XYZ)), zeros(size(XYZ));
    for (i,val) in enumerate(XYZ)
        X[i],Y[i],Z[i] = val.e, val.n, val.u;
    end

    # Convert full structure to cartesian model
    data_Cart = CartData(GeoUnit(X*km,km),GeoUnit(Y*km,km),GeoUnit(Z*km,km),data.fields)

    if Flatten==true
        data_Cart.z = data.depth
    end

    return data_Cart, refPoint
end


"""
    data = Cartesian_To_GeoData(data_Cart::CartData, refPoint::LLA; Flatten=true)

Projects a local Cartesian structure `data_Cart` to a GeoData structure `data_Cart` (with lon/lat/depth), using the reference point `refPoint`.
The refPoint needs to be in LLA  (latitude, longitude, altitude) format, which can be generated with:
```julia
julia> refPoint = LLA(lat,lon,alt)
```
You also need to specify whether the Cartesian domain was `flattened` (by default this is assumed to be the case, as this is how typical geodynamic models work).

"""
function Cartesian_To_GeoData(data_Cart::CartData, refPoint::LLA; Flatten=true)
    
    inverse     =   LLAfromENU(refPoint, wgs84);      
    ENU_Data    =   ENU.(ustrip(data_Cart.x.val)*1e3, ustrip(data_Cart.y.val)*1e3, ustrip(data_Cart.z.val)*1e3); # ENU structure
    LLA_Data    =   inverse.(ENU_Data);

    Lon,Lat,Depth       =   zeros(size(LLA_Data)), zeros(size(LLA_Data)), zeros(size(LLA_Data));
    for (i,val) in enumerate(LLA_Data)
        Lon[i],Lat[i],Depth[i] = val.lon, val.lat, val.alt;
    end

    # in case it was flattened, unflatten it:
    if Flatten==true
        Depth   =   data_Cart.z.val;
    else
        Depth   =   Depth/1e3;      # GeoData wants depth in km
    end

    data    =   GeoData(Lon, Lat, Depth, data_Cart.fields)

    return data
end

"""
    data_Rect = InterpolateToRectilinear(data_Cart::CartData; x_lims=empty, y_lims=empty, z_lims=empty, dims=empty)

This interpolates a Cartesian (but structured) grid `data_Cart`, with (slightly) distorted `x,y` coordinates, to a rectilinear grid with a constant spacing.

`x_lims=[min(x), max(x)]`, `y_lims=[min(y), max(y)]`, `z_lims=[min(z), max(z)]` are optional min/max limits of the new grid and `dims=(nx,ny,nz)` is a tuple that indicates the resolution of the grid. If some of these parameters are nt indicated, they are computed from `data_Cart`. 


NOT FINISHED YET!

"""
function InterpolateToRectilinear(data_Cart::CartData; x_lims=empty, y_lims=empty, z_lims=empty, dims=empty)

    # if we do not set values for the domain, retrieve them such that the new grid is fully covered by `data_Cart` (to avoid interpolation errors)
    if x_lims==empty
        x_lims = [maximum(data_Cart.x.val[1,:,:]), minimum(data_Cart.x.val[end,:,:])]
    end
    if y_lims==empty
        y_lims = [maximum(data_Cart.y.val[:,1,:]), minimum(data_Cart.y.val[:,end,:])]
    end
    if z_lims==empty
        z_lims = [maximum(data_Cart.z.val[:,:,1]), minimum(data_Cart.z.val[:,:,end])]
    end 
    if dims==empty
        dims = size(data_Cart.x)
    end

    x,y,z = LonLatDepthGrid(LinRange(x_lims[1],x_lims[2],dims[1]),
                            LinRange(y_lims[1],y_lims[2],dims[2]),
                            LinRange(z_lims[1],z_lims[2],dims[3]))


    return x,y,z
end
