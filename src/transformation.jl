# This provides various transformations (GeoData <=> Cartesian; UTMData <=> Cartesian)
#
using StaticArrays

export project_CartData, project_FEData_CartData



"""
    d_cart = project_CartData(d_cart::CartData, d::GeoData, p::ProjectionPoint)

Projects all datafields from the GeoData struct `d` to the CartData struct `d_cart`, around the projection point `p`.
`d_cart` *must* be an orthogonal cartesian grid (deformed doesn't work; use `convert2CartData(d, proj)`, where `proj` is a projection point in that case).

# Note:    
- If `d_cart` and `d` are horizontal surfaces (3rd dimension has size==1), it also interpolates the depth coordinate.    

"""
function project_CartData(d_cart::CartData, d::GeoData, p::ProjectionPoint)
    Data_UTM    = convert2UTMzone(d_cart, p)
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
        z_new, fields_new = interpolate_datafields_2D(d,Data_lonlat.lon.val, Data_lonlat.lat.val)
      
        # Create new struct 
        d_cart = CartData(d_cart.x.val,d_cart.y.val,z_new,fields_new)
    
    else
        d_data = interpolate_datafields(d, Data_lonlat.lon.val, Data_lonlat.lat.val, Data_lonlat.depth.val)    
        d_cart = CartData(d_cart.x.val,d_cart.y.val,d_cart.z.val,d_data.fields)
    
    end

    return d_cart
end

"""
    d_cart = project_CartData(d_cart::CartData, d::GeoData, p::ProjectionPoint)

Projects all datafields from the GeoData struct `d` to the CartData struct `d_cart`, around the projection point `p`.
`d_cart` *must* be an orthogonal cartesian grid (deformed doesn't work; use `convert2CartData(d, proj)`, where `proj` is a projection point in that case).

# Note:    
- If `d_cart` and `d` are horizontal surfaces (3rd dimension has size==1), it also interpolates the depth coordinate.    

"""
function project_CartData(d_cart::CartData, d_cart_data0::CartData)
    
    if size(d_cart_data0.x.val,3)==1
        z_new, fields_new = interpolate_datafields_2D(d,d_cart_data0.x.val, d_cart_data0.y.val)
        
        # Create new struct 
        d_cart = CartData(d_cart.x.val,d_cart.y.val,z_new,fields_new)
    
    else
        d_data = interpolate_datafields(d, d_cart_data0.x.val, d_cart_data0.y.val, d_cart_data0.z.val)    
        d_cart = CartData(d_cart.x.val,d_cart.y.val,d_cart.z.val,d_data.fields)
    
    end

    return d_cart
end


"""
    d_cart = project_CartData(d_cart::CartData, d::UTMData, p::ProjectionPoint)

Projects all datafields from the UTMData struct `d` to the CartData struct `d_cart`, around the projection point `p`.
    `d_cart` *must* be an orthogonal cartesian grid (deformed doesn't work; use `convert2CartData(d, proj)`, where `proj` is a projection point in that case).
    
    # Note:    
    - If `d_cart` and `d` are horizontal surfaces (3rd dimension has size==1), it also interpolates the depth coordinate.    
        

"""
function project_CartData(d_cart::CartData, d::UTMData, p::ProjectionPoint)
    Data_UTM    = convert2UTMzone(d_cart, p)
    
    if size(Data_UTM.EW.val,3)==1
        z_new, fields_new = interpolate_datafields_2D(d,Data_UTM.EW.val, Data_UTM.NS.val)
        
        # Create new struct 
        d_cart = CartData(d_cart.x.val,d_cart.y.val,z_new,fields_new)
    
    else
        d_data = interpolate_datafields(d, Data_UTM.EW.val, Data_UTM.NS.val, Data_UTM.depth.val)    
        d_cart = CartData(d_cart.x.val,d_cart.y.val,d_cart.z.val,d_data.fields)
    
    end

    return d_cart
end



"""
    inside = point_in_tetrahedron(p::_T, a::_T, b::_T, c::_T, d::_T, tol=1e-10)
Determines if a point `p` is inside a tetrahedron specified by `a`,`b`,`c`,`d` or not    
"""
function point_in_tetrahedron(p::_T, a::_T, b::_T, c::_T, d::_T, tol=1e-10) where _T<:Vector{Float64}

    # check bounding box
    xmin = min(a[1],b[1],c[1],d[1])
    xmax = max(a[1],b[1],c[1],d[1])
    ymin = min(a[2],b[2],c[2],d[2])
    ymax = max(a[2],b[2],c[2],d[2])
    zmin = min(a[3],b[3],c[3],d[3])
    zmax = max(a[3],b[3],c[3],d[3])
    
    inside = true
    if p[1] < xmin || p[1] > xmax 
        inside = false
    end
    if (p[2] < ymin || p[2] > ymax) && inside  
        inside = false
    end
    if (p[3] < zmin || p[3] > zmax) && inside  
        inside = false
    end
    
    if inside
        v0 = @SVector [d[i] - a[i] for i in 1:3]
        v1 = @SVector [b[i] - a[i] for i in 1:3]
        v2 = @SVector [c[i] - a[i] for i in 1:3]
        v3 = @SVector [p[i] - a[i] for i in 1:3]
    
        denom = dot(v0, cross(v1, v2))
    
        u = dot(v3, cross(v1, v2)) / denom
        v = dot(v0, cross(v3, v2)) / denom
        w = dot(v0, cross(v1, v3)) / denom
    
        inside =  (u >= -tol) && (v >= -tol) && (w >= -tol) && (u + v + w <= 1 + tol)
    end

    return inside
end

"""
    data_cart = project_FEData_CartData(data_cart::CartData, data_fe::FEData)

Projects a FEData object with tetrahedrons (e.g., from Gmsh) to a Cartesian grid
"""
function project_FEData_CartData(data_cart::CartData, data_fe::FEData)

    cellfields_regions = data_fe.cellfields.regions
    regions = zeros(Int64, size(data_cart.x.val))
    
    for i = 1:size(data_fe.connectivity,2) # loop over tetrahedrons
        tetra = data_fe.connectivity[:,i]

        a = data_fe.vertices[:, tetra[1]]
        b = data_fe.vertices[:, tetra[2]]
        c = data_fe.vertices[:, tetra[3]]
        d = data_fe.vertices[:, tetra[4]]

        xmin, xmax = extrema((a[1],b[1],c[1],d[1]))
        ymin, ymax = extrema((a[2],b[2],c[2],d[2]))
        zmin, zmax = extrema((a[3],b[3],c[3],d[3]))

        ind = findall(@. xmax ≥ data_cart.x.val ≥ xmin && 
                         ymax ≥ data_cart.y.val ≥ ymin && 
                         zmax ≥ data_cart.z.val ≥ zmin)

        for I in ind
                x = data_cart.x.val[I]
                y = data_cart.y.val[I]
                z = data_cart.z.val[I]
                p = [x,y,z]
                if point_in_tetrahedron(p,a,b,c,d)
                    regions[I] = cellfields_regions[i]
                end
        end

    end

    return addfield(data_cart, (regions=regions,))
end



