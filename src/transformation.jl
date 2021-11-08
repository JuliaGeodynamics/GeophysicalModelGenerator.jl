# This provides various transformations (GeoData <=> Cartesian)
#

export ProjectCartData


"""
    d_cart = ProjectCartData(d_cart::CartData, d::GeoData, p::ProjectionPoint)

Projects all datafields from the GeoData struct `d` to the CartData struct `d_cart`, around the projection point `p`.

# Note:    
- If `d_cart` and `d` are horizontal surfaces (3rd dimension has size==1), it also interpolates the depth coordinate.    

"""
function ProjectCartData(d_cart::CartData, d::GeoData, p::ProjectionPoint)
    Data_UTM    = Convert2UTMzone(d_cart, p)
    Data_lonlat = convert(GeoData,Data_UTM)   
    
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