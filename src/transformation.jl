# This provides various transformations (GeoData <=> Cartesian; UTMData <=> Cartesian)
#

export ProjectCartData


"""
    d_cart = ProjectCartData(d_cart::CartData, d::GeoData, p::ProjectionPoint)

Projects all datafields from the GeoData struct `d` to the CartData struct `d_cart`, around the projection point `p`.
`d_cart` *must* be an orthogonal cartesian grid (deformed doesn't work; use `Convert2CartData(d, proj)`, where `proj` is a projection point in that case).

# Note:    
- If `d_cart` and `d` are horizontal surfaces (3rd dimension has size==1), it also interpolates the depth coordinate.    

"""
function ProjectCartData(d_cart::CartData, d::GeoData, p::ProjectionPoint)
    Data_UTM    = Convert2UTMzone(d_cart, p)
    Data_lonlat = convert(GeoData,Data_UTM)   
    
    # Check whether the data sets have the same sign. If not, we may have to shift one by 360 degrees
    min_lon_cart, max_lon_cart = minimum(Data_lonlat.lon.val), maximum(Data_lonlat.lon.val)
    min_lon, max_lon = minimum(d.lon.val), maximum(d.lon.val)
    
    if (sign(min_lon)!=sign(min_lon_cart)) &&  (sign(max_lon)!=sign(max_lon_cart))
        # the longitude data has a different sign. This can happen if one of them is "West" (and has negative values), whereas the other has 
        if (min_lon_cart<0)
            Data_lonlat = GeoData(Data_lonlat.lon.val .+ 360, Data_lonlat.lat.val, ustrip.(Data_lonlat.depth.val), Data_lonlat.fields)
        end
    end
    
    if size(Data_lonlat.lon.val,3)==1
        z_new, fields_new = InterpolateDataFields2D(d,Data_lonlat.lon.val, Data_lonlat.lat.val)
      
        # Create new struct 
        d_cart = CartData(d_cart.x.val,d_cart.y.val,z_new,fields_new)
    
    else
        d_data = InterpolateDataFields(d, Data_lonlat.lon.val, Data_lonlat.lat.val, Data_lonlat.depth.val)    
        d_cart = CartData(d_cart.x.val,d_cart.y.val,d_cart.z.val,d_data.fields)
    
    end

    return d_cart
end



"""
    d_cart = ProjectCartData(d_cart::CartData, d::UTMData, p::ProjectionPoint)

Projects all datafields from the UTMData struct `d` to the CartData struct `d_cart`, around the projection point `p`.
    `d_cart` *must* be an orthogonal cartesian grid (deformed doesn't work; use `Convert2CartData(d, proj)`, where `proj` is a projection point in that case).
    
    # Note:    
    - If `d_cart` and `d` are horizontal surfaces (3rd dimension has size==1), it also interpolates the depth coordinate.    
        

"""
function ProjectCartData(d_cart::CartData, d::UTMData, p::ProjectionPoint)
    Data_UTM    = Convert2UTMzone(d_cart, p)
    
    if size(Data_UTM.EW.val,3)==1
        z_new, fields_new = InterpolateDataFields2D(d,Data_UTM.EW.val, Data_UTM.NS.val)
        
        # Create new struct 
        d_cart = CartData(d_cart.x.val,d_cart.y.val,z_new,fields_new)
    
    else
        d_data = InterpolateDataFields(d, Data_UTM.EW.val, Data_UTM.NS.val, Data_UTM.depth.val)    
        d_cart = CartData(d_cart.x.val,d_cart.y.val,d_cart.z.val,d_data.fields)
    
    end

    return d_cart
end