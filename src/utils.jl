# few utils that are useful 

export meshgrid, CrossSection, ExtractSubvolume, SubtractHorizontalMean, Flatten3DData
export ParseColumns_CSV_File, AboveSurface, BelowSurface, VoteMap, InterpolateDataOnSurface
export RotateTranslateScale!

"""
    meshgrid(vx,vy,vz)

Computes an (x,y,z)-grid from the vectors (vx,vy,vz).
For more information, see the MATLAB documentation.
"""
function meshgrid(vx::AbstractVector{T}, vy::AbstractVector{T},
                     vz::AbstractVector{T}) where {T}
    m, n, o = length(vy), length(vx), length(vz)
    vx = reshape(vx, 1, n, 1)
    vy = reshape(vy, m, 1, 1)
    vz = reshape(vz, 1, 1, o)
    om = ones(Int, m)
    on = ones(Int, n)
    oo = ones(Int, o)
    (vx[om, :, oo], vy[:, on, oo], vz[om, on, :])
end

"""
    CrossSection(Volume::GeoData; dims=(100,100), Interpolate=false, Depth_level=nothing; Lat_level=nothing; Lon_level=nothing; Start=nothing, End=nothing )

Creates a cross-section through a volumetric (3D) `GeoData` object. 

- Cross-sections can be horizontal (map view at a given depth), if `Depth_level` is specified
- They can also be vertical, either by specifying `Lon_level` or `Lat_level` (for a fixed lon/lat), or by defining both `Start=(lon,lat)` & `End=(lon,lat)` points.
- `Interpolate` indicates whether we want to simply extract the data from the 3D volume (default) or whether we want to linearly interpolate it on a new grid, which has dimensions as specified in `dims`

# Example:
```julia-repl
julia> Lon,Lat,Depth   =   LonLatDepthGrid(10:20,30:40,(-300:25:0)km);
julia> Data            =   Depth*2;                # some data
julia> Vx,Vy,Vz        =   ustrip(Data*3),ustrip(Data*4),ustrip(Data*5);
julia> Data_set3D      =   GeoData(Lon,Lat,Depth,(Depthdata=Data,LonData=Lon, Velocity=(Vx,Vy,Vz))); 
julia> Data_cross      =   CrossSection(Data_set3D, Depth_level=-100km)  
GeoData 
  size  : (11, 11, 1)
  lon   ϵ [ 10.0 : 20.0]
  lat   ϵ [ 30.0 : 40.0]
  depth ϵ [ -100.0 km : -100.0 km]
  fields: (:Depthdata, :LonData, :Velocity)
```

"""
function CrossSection(V::GeoData; dims=(100,100), Interpolate=false, Depth_level=nothing, Lat_level=nothing, Lon_level=nothing, Start=nothing, End=nothing )

    CheckIsVolume(V);       # Check if it is a volume

    if !isnothing(Depth_level)    # Horizontal slice
        CheckBounds(V.depth, Depth_level)    
        if Interpolate
            Lon,Lat,Depth = LonLatDepthGrid(    LinRange(minimum(V.lon.val), maximum(V.lon.val), dims[1]),
                                                LinRange(minimum(V.lat.val), maximum(V.lat.val), dims[2]),
                                                Depth_level)
        else
            ind_z   =   argmin(abs.(V.depth.val[1,1,:] .- Depth_level))
            iDepth  =   ind_z:ind_z;
            iLon    =   1:size(V.lon.val,1);
            iLat    =   1:size(V.lat.val,2);
        end
    end

    if !isnothing(Lat_level)   # vertical slice @ given latitude
        CheckBounds(V.lat, Lat_level)    
        if Interpolate
            Lon,Lat,Depth = LonLatDepthGrid(    LinRange(minimum(V.lon.val), maximum(V.lon.val), dims[1]),
                                                Lat_level,
                                                LinRange(minimum(V.depth.val), maximum(V.depth.val), dims[2]))
        else
            ind_l   =   argmin(abs.(V.lat.val[1,:,1] .- Lat_level))
            iDepth  =   1:size(V.depth.val,3)
            iLon    =   1:size(V.lon.val,1);
            iLat    =   ind_l:ind_l
        end
    end

    if !isnothing(Lon_level)   # vertical slice @ given longitude
        CheckBounds(V.lon, Lon_level)    
        if Interpolate 
            Lon,Lat,Depth = LonLatDepthGrid(    Lon_level,
                                                LinRange(minimum(  V.lat.val), maximum(  V.lat.val), dims[1]),
                                                LinRange(minimum(V.depth.val), maximum(V.depth.val), dims[2]))
        else
            ind_l   =   argmin(abs.(V.lon.val[:,1,1] .- Lon_level))
            iDepth  =   1:size(V.depth.val,3)
            iLat    =   1:size(V.lat.val,2);
            iLon    =   ind_l:ind_l
        end
    end

    # diagonal profile defined by start and end lon/lat points
    if !isnothing(Start)
        if isnothing(End)
            error("Also define End coordinates if you indicate starting lon/lat value")
        end
        Interpolate = true; # we must interpolate in this case

        Lon_dum,Lat_p,Depth_p = LonLatDepthGrid(    Start[1],
                                                LinRange(Start[2], End[2], dims[1]),
                                                LinRange(minimum(V.depth.val), maximum(V.depth.val), dims[2]))

        Lon_p,Lat_dum,Depth = LonLatDepthGrid(    LinRange(Start[1], End[1], dims[1]),
                                                Start[2],
                                                LinRange(minimum(V.depth.val), maximum(V.depth.val), dims[2]))

        Lon             =   zeros(dims[1],dims[2],1)
        Lat             =   zeros(dims[1],dims[2],1)
        Depth           =   zeros(dims[1],dims[2],1)*Depth_p[1]
        
        # We need 3D matrixes for the paraview writing routine to know we are in 3D
        Lon[:,:,1]      =   Lon_p[:,1,:]
        Lat[:,:,1]      =   Lat_p[1,:,:]
        Depth[:,:,1]    =   Depth_p[1,:,:]
        
    end

    if Interpolate
        # Interpolate data on profile
        DataProfile = InterpolateDataFields(V, Lon, Lat, Depth);    
    else
        # extract data (no interpolation)
        DataProfile = ExtractDataSets(V, iLon, iLat, iDepth);
    end

    return DataProfile

end


"""
    ExtractSubvolume(V::GeoData; Interpolate=false, Lon_level=nothing, Lat_level=nothing, Depth_level=nothing, dims=(50,50,50))

Extract or "cuts-out" a piece of a 2D or 3D GeoData set, defined by `Lon`, `Lat` and `Depth` coordinates.

This is useful if you are only interested in a part of a much bigger larger data set.

- `Lon_level`,`Lat_level` and `Depth_level` should be tuples that indicate `(minimum_value, maximum_value)` along the respective direction. If not specified we use the full range. 
- By default, `Interpolate=false` and we find the closest indices within the data set (so your new data set will not go exactly from minimum to maximum).
- Alternatively, if `Interpolate=true` we interpolate the data onto a new grid that has dimensions `dims`. This can be useful to compare data sets that are originally given in different resolutions.

# 3D Example with no interpolation:
```julia-repl
julia> Lon,Lat,Depth   =   LonLatDepthGrid(10:20,30:40,(-300:25:0)km);
julia> Data            =   Depth*2;                # some data
julia> Vx,Vy,Vz        =   ustrip(Data*3),ustrip(Data*4),ustrip(Data*5);
julia> Data_set3D      =   GeoData(Lon,Lat,Depth,(Depthdata=Data,LonData=Lon, Velocity=(Vx,Vy,Vz)))
GeoData 
  size  : (11, 11, 13)
  lon   ϵ [ 10.0 : 20.0]
  lat   ϵ [ 30.0 : 40.0]
  depth ϵ [ -300.0 km : 0.0 km]
  fields: (:Depthdata, :LonData, :Velocity)
julia> Data_extracted = ExtractSubvolume(Data_set3D,Lon_level=(10,12),Lat_level=(35,40))
GeoData 
  size  : (3, 6, 13)
  lon   ϵ [ 10.0 : 12.0]
  lat   ϵ [ 35.0 : 40.0]
  depth ϵ [ -300.0 km : 0.0 km]
  fields: (:Depthdata, :LonData, :Velocity)
```
By default it extracts the data points closest to the area defined by Lon_level/Lat_level/Depth_level.

# 3D Example with interpolation:
Alternatively, you can also interpolate the data onto a new grid:
```julia
julia> Data_extracted = ExtractSubvolume(Data_set3D,Lon_level=(10,12),Lat_level=(35,40), Interpolate=true, dims=(50,51,52))
GeoData 
  size  : (50, 51, 52)
  lon   ϵ [ 10.0 : 12.0]
  lat   ϵ [ 35.0 : 40.0]
  depth ϵ [ -300.0 km : 0.0 km]
  fields: (:Depthdata, :LonData, :Velocity)
```

"""
function ExtractSubvolume(V::GeoData; Interpolate=false, Lon_level=nothing, Lat_level=nothing, Depth_level=nothing, dims=(50,50,50))

    if isnothing(Lon_level)
        Lon_level = (minimum(V.lon.val), maximum(V.lon.val))
    end
    if isnothing(Lat_level)
        Lat_level = (minimum(V.lat.val), maximum(V.lat.val))
    end
    if isnothing(Depth_level)
        Depth_level = (minimum(V.depth.val), maximum(V.depth.val))
    end
    if Interpolate
        Lon,Lat,Depth   = LonLatDepthGrid(  LinRange(Lon_level[1],      Lon_level[2],   dims[1]),
                                            LinRange(Lat_level[1],      Lat_level[2],   dims[2]),
                                            LinRange(Depth_level[1],    Depth_level[2], dims[3]) );
        Data_extract    =   InterpolateDataFields(V, Lon, Lat, Depth)

    else
        # Don't interpolate
        i_s, i_e    =   argmin(abs.(V.lon.val[:,1,1] .- Lon_level[1])), argmin(abs.(V.lon.val[:,1,1] .- Lon_level[2]))
        iLon        =   i_s:i_e;
        
        i_s, i_e    =   argmin(abs.(V.lat.val[1,:,1] .- Lat_level[1])), argmin(abs.(V.lat.val[1,:,1] .- Lat_level[2]))
        iLat        =   i_s:i_e;
        
        i_s, i_e    =   argmin(abs.(V.depth.val[1,1,:] .- Depth_level[1])), argmin(abs.(V.depth.val[1,1,:] .- Depth_level[2]))
        step        =   1;
        if i_e<i_s
            step=-1
        end
        iDepth      =   i_s:step:i_e;
        Data_extract =  ExtractDataSets(V, iLon, iLat, iDepth);
    end

    return Data_extract
end


function CheckBounds(Data::GeoUnit, Data_Cross)
    
    min_Data, max_Data = minimum(Data.val), maximum(Data.val);
    if Data_Cross < min_Data || Data_Cross>max_Data
        error("Outside bounds [$min_Data : $max_Data]; $Data_Cross")
    end
end


function CheckIsVolume(Volume::GeoData)
    if any(size(Volume.lon).==1)
        error("It appears your GeoData structure is not a volume; can't compute a cross-section")
    end
end

"""
    InterpolateDataFields(V::GeoData, Lon, Lat, Depth)

Interpolates a data field `V` on a grid defined by `Lon,Lat,Depth`
"""
function InterpolateDataFields(V::GeoData, Lon, Lat, Depth)

    Lon_vec     =  V.lon.val[:,1,1];
    Lat_vec     =  V.lat.val[1,:,1];
    Depth_vec   =  V.depth.val[1,1,:];
    if Depth_vec[1]>Depth_vec[end]
        ReverseData = true
    else
        ReverseData = false
    end

    fields_new  = V.fields;
    field_names = keys(fields_new);
    for i = 1:length(V.fields)
        if typeof(V.fields[i]) <: Tuple
            # vector or anything that contains more than 1 field
            data_tuple = fields_new[i]      # we have a tuple (likely a vector field), so we have to loop 
            data_array = zeros(size(Lon,1),size(Lon,2),size(Lon,3),length(data_tuple));     # create a 3D array that holds the 2D interpolated values
            unit_array = zeros(size(data_array));

            for j=1:length(data_tuple)
                if ReverseData
                    ndim        =   length(size(data_tuple[j]))
                    interpol    =   LinearInterpolation((Lon_vec, Lat_vec, reverse(Depth_vec)), reverse(ustrip.(data_tuple[j]), dims=ndim) ,extrapolation_bc = Flat());      # create interpolation object
                else
                    interpol    =   LinearInterpolation((Lon_vec, Lat_vec, Depth_vec), ustrip.(data_tuple[j]),extrapolation_bc = Flat());      # create interpolation object
                end
                data_array[:,:,:,j] =   interpol.(Lon, Lat, Depth);          
            end
            data_new    = tuple([data_array[:,:,:,c] for c in 1:size(data_array,4)]...)     # transform 3D matrix to tuple

        else
            # scalar field
            if ReverseData
                ndim        =   length(size(V.fields[i]))
                interpol    =   LinearInterpolation((Lon_vec, Lat_vec, reverse(Depth_vec)), reverse(V.fields[i], dims=ndim), extrapolation_bc = Flat(),);            # create interpolation object
            else
                interpol    =   LinearInterpolation((Lon_vec, Lat_vec, Depth_vec), V.fields[i], extrapolation_bc = Flat());            # create interpolation object
            end
            data_new    =   interpol.(Lon, Lat, Depth);                                                 # interpolate data field
        end
        
        # replace the one 
        new_field   =   NamedTuple{(field_names[i],)}((data_new,))                          # Create a tuple with same name
        fields_new  =   merge(fields_new, new_field);                                       # replace the field in fields_new
        
    end
    

    # Create a GeoData struct with the newly interpolated fields
    Data_profile = GeoData(Lon, Lat, Depth, fields_new);

    return Data_profile
end


"""
    Surf_interp = InterpolateDataOnSurface(V::CartData, Surf::CartData)

Interpolates a 3D data set `V` on a surface defined by `Surf`. nex
# Example
```julia
julia> Data
CartData 
  size  : (33, 33, 33)
  x     ϵ [ -3.0 : 3.0]
  y     ϵ [ -2.0 : 2.0]
  z     ϵ [ -2.0 : 0.0]
  fields: (:phase, :density, :visc_total, :visc_creep, :velocity, :pressure, :temperature, :dev_stress, :strain_rate, :j2_dev_stress, :j2_strain_rate, :plast_strain, :plast_dissip, :tot_displ, :yield, :moment_res, :cont_res)
julia> surf
CartData 
  size  : (96, 96, 1)
  x     ϵ [ -2.9671875 : 3.2671875]
  y     ϵ [ -1.9791666666666667 : 1.9791666666666667]
  z     ϵ [ -1.5353766679763794 : -0.69925457239151]
  fields: (:Depth,)
julia> Surf_interp = InterpolateDataOnSurface(Data, surf)
  CartData 
    size  : (96, 96, 1)
    x     ϵ [ -2.9671875 : 3.2671875]
    y     ϵ [ -1.9791666666666667 : 1.9791666666666667]
    z     ϵ [ -1.5353766679763794 : -0.69925457239151]
    fields: (:phase, :density, :visc_total, :visc_creep, :velocity, :pressure, :temperature, :dev_stress, :strain_rate, :j2_dev_stress, :j2_strain_rate, :plast_strain, :plast_dissip, :tot_displ, :yield, :moment_res, :cont_res)
```
"""
function InterpolateDataOnSurface(V::CartData, Surf::CartData)
    
    # Create GeoData structure:
    V_geo               =   GeoData(V.x.val, V.y.val, V.z.val, V.fields)
    V_geo.depth.val     =   ustrip(V_geo.depth.val);

    Surf_geo            =   GeoData(Surf.x.val, Surf.y.val, Surf.z.val, Surf.fields)
    Surf_geo.depth.val  =   ustrip(Surf_geo.depth.val);

    Surf_interp_geo     =   InterpolateDataOnSurface(V_geo, Surf_geo)
    Surf_interp         =   CartData(Surf_interp_geo.lon.val, Surf_interp_geo.lat.val, ustrip.(Surf_interp_geo.depth.val), Surf_interp_geo.fields)

    return Surf_interp

end


"""
    Surf_interp = InterpolateDataOnSurface(V::GeoData, Surf::GeoData)

Interpolates a 3D data set `V` on a surface defined by `Surf`
"""
function InterpolateDataOnSurface(V::GeoData, Surf::GeoData)
    
    Surf_interp = InterpolateDataFields(V, Surf.lon.val, Surf.lat.val, Surf.depth.val)

    return Surf_interp
end



# Extracts a sub-data set using indices
function ExtractDataSets(V::GeoData, iLon, iLat, iDepth)

    Lon     =   zeros(typeof(V.lon.val[1]), length(iLon),length(iLat),length(iDepth));
    Lat     =   zeros(typeof(V.lat.val[1]), length(iLon),length(iLat),length(iDepth));
    Depth   =   zeros(typeof(V.depth.val[1]), length(iLon),length(iLat),length(iDepth));
    
    iLo                 =   1:length(iLon);
    iLa                 =   1:length(iLat);
    iDe                 =   1:length(iDepth)
    Lon[iLo,iLa,iDe]    =     V.lon.val[iLon, iLat, iDepth];
    Lat[iLo,iLa,iDe]    =     V.lat.val[iLon, iLat, iDepth];
    Depth[iLo,iLa,iDe]  =   V.depth.val[iLon, iLat, iDepth];

    fields_new  = V.fields;
    field_names = keys(fields_new);
    for i = 1:length(V.fields)
        if typeof(V.fields[i]) <: Tuple
            # vector or anything that contains more than 1 field
            data_tuple = fields_new[i]      # we have a tuple (likely a vector field), so we have to loop 
            data_array = zeros(typeof(data_tuple[1][1]),length(iLon),length(iLat),length(iDepth),length(data_tuple));     # create a 3D array that holds the 2D interpolated values
            unit_array = zeros(size(data_array));

            for j=1:length(data_tuple)
                data_field           =   data_tuple[j];
                data_array[:,:,:,j]  =   data_field[iLon, iLat, iDepth];          
            end
            data_new    = tuple([data_array[:,:,:,c] for c in 1:size(data_array,4)]...)       # transform 4D matrix to tuple

        else
            # scalar field
            data_new                =   zeros(typeof(V.fields[i][1]), length(iLon),length(iLat),length(iDepth));
            data_new[iLo,iLa,iDe]   =   V.fields[i][iLon, iLat, iDepth]                                 # interpolate data field
        end
        
        # replace the one 
        new_field   =   NamedTuple{(field_names[i],)}((data_new,))                          # Create a tuple with same name
        fields_new  =   merge(fields_new, new_field);                                       # replace the field in fields_new
        
    end
    

    # Create a GeoData struct with the newly interpolated fields
    Data_profile = GeoData(Lon, Lat, Depth, fields_new);

end

"""
    V_sub = SubtractHorizontalMean(V::AbstractArray{T, 3}; Percentage=false)

Subtracts the horizontal average of the 3D data array V.

If `Percentage=true`, the result is given as percentage; otherwise absolute values are returned

"""
function SubtractHorizontalMean(V::AbstractArray{T, 3}; Percentage=false) where T

    nx        = size(V,1);
    ny        = size(V,2);
    NumLayers = size(V,3); # get the number of depth levels

    if Percentage
        V_sub     = zeros(size(V));                 # no units
    else
        V_sub     = zeros(typeof(V[1]), size(V));   
    end

    for iLayer = 1:NumLayers
        average             =   mean(filter(!isnan, vec(V[:,:,iLayer])));
        
        if Percentage
            V_sub[:,:,iLayer]   =   ustrip(V[:,:,iLayer]) .- ustrip(average);
            V_sub[:,:,iLayer]   =   V_sub[:,:,iLayer]./ustrip(average)*100.0;     # the result is normalized 
        else
            V_sub[:,:,iLayer]   =   V[:,:,iLayer] .- average;
        end
    end

    return V_sub
end

"""
    V_sub = SubtractHorizontalMean(V::AbstractArray{T, 2}; Percentage=false)

Subtracts the horizontal average of the 2D data array V.

If `Percentage=true`, the result is given as percentage; otherwise absolute values are returned

"""
function SubtractHorizontalMean(V::AbstractArray{T, 2}; Percentage=false) where T

    nx        = size(V,1);
    NumLayers = size(V,2); # get the number of depth levels

    if Percentage
        V_sub     = zeros(size(V));                 # no units
    else
        V_sub     = zeros(typeof(V[1]), size(V));   
    end

    for iLayer = 1:NumLayers
        average             =   mean(filter(!isnan, vec(V[:,iLayer])));
        
        if Percentage
            V_sub[:,iLayer]   =   ustrip(V[:,iLayer]) .- ustrip(average);
            V_sub[:,iLayer]   =   V_sub[:,iLayer]./ustrip(average)*100.0;     # the result is normalized 
        else
            V_sub[:,iLayer]   =   V[:,iLayer] .- average;
        end
    end

    return V_sub
end


# "flatten" a GeoData input to obtain x/y/z values
function Flatten3DData(Data::GeoData)

    ndepth = size(Data.lat.val,3)
    lat = Data.lat.val[:,:,1]
    lon = Data.lon.val[:,:,1]

    # origin
    xo_lla = LLA.(ones(size(lat)).*minimum(vec(lat)), ones(size(lon)).*minimum(vec(lon)), zeros(size(lat))); # convert to LLA format

    # compute x-coordinates by computing the eucledian distance between longitudes, but setting the latitudes to the same value
    x_lla = LLA.(ones(size(lat)).*minimum(vec(lat)), lon, zeros(size(lat))); # convert to LLA format
    x_coord = euclidean_distance.(x_lla, xo_lla);
    x_lla = LLA.(lat, ones(size(lon)).*minimum(vec(lon)), zeros(size(lat))); # convert to LLA format
    y_coord = euclidean_distance.(x_lla, xo_lla);

    # now create matrices from that (convert m to km as this is the internal standard)
    X = repeat(x_coord./1e3,1,1,ndepth);
    Y = repeat(y_coord./1e3,1,1,ndepth);
    Z = ustrip(Data.depth.val);

    DataFlat = CartData(X, Y, Z, Data.fields)
    return DataFlat
end


""" 
    ParseColumns_CSV_File(data_file, num_columns)

This parses numbers from CSV file that is read in with `CSV.File`.
That is useful in case the CSV files has tables that contain both strings (e.g., station names) and numbers (lat/lon/height) and you are only intested in the numbers


# Example
This example assumes that the data starts at line 18, that the colums are separated by spaces, and that it contains at most 4 columns with data:
```julia-repl
julia> using CSV
julia> data_file        =   CSV.File("FileName.txt",datarow=18,header=false,delim=' ')
julia> data = ParseColumns_CSV_File(data_file, 4)
```

"""
function ParseColumns_CSV_File(data_file, num_columns)
    data                =   zeros(size(data_file,1), num_columns);    
    for (row_num,row) in enumerate(data_file)
        num         =   0;
        for i=1:length(row)
            if typeof(row[i])==Float64
                num          +=  1;
                data[row_num,num] = row[i]
            else
            
                try parse(Float64,row[i])
                    num          +=  1;
                    data[row_num,num] = parse(Float64,row[i])
                catch
                end
            end

        end
    end
    return data
end


""" 
    AboveSurface(Data::GeoData, DataSurface::GeoData; above=true)

Returns a boolean array of size(Data.Lon), which is true for points that are above the surface DataSurface (or for points below if above=false).

This can be used, for example, to mask points above/below the Moho in a volumetric dataset or in a profile.

#Example
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
julia> Above       =   AboveSurface(Data_set3D, Data_Moho); 
```
Now, `Above` is a boolean array that is true for points above the surface and false for points below and at the surface.

"""
function AboveSurface(Data::GeoData, DataSurface::GeoData; above=true)
    
    if size(DataSurface.lon)[3]!=1
        error("It seems that DataSurface is not a surface")
    end

    # Create interpolation object for surface
    Lon_vec     =  DataSurface.lon.val[:,1,1];
    Lat_vec     =  DataSurface.lat.val[1,:,1];
    interpol    =  LinearInterpolation((Lon_vec, Lat_vec), ustrip.(DataSurface.depth.val[:,:,1]));            # create interpolation object

    DepthSurface = interpol.(Data.lon.val,Data.lat.val);
    DepthSurface = DepthSurface*unit(DataSurface.depth.val[1])

    if above
        Above       =   Data.depth.val .> DepthSurface;
    else
        Above       =   Data.depth.val .< DepthSurface;
    end

    return Above
end

"""
    Below = BelowSurface(Data::GeoData, DataSurface::GeoData)

Determines if points within the 3D `Data` structure are below the GeoData surface `DataSurface`
"""
function BelowSurface(Data::GeoData, DataSurface::GeoData)
    return AboveSurface(Data::GeoData, DataSurface::GeoData; above=false)
end

"""
    Above = AboveSurface(Data_Cart::CartData, DataSurface_Cart::CartData; above=true)

Determines if points within the 3D `Data_Cart` structure are above the Cartesian surface `DataSurface_Cart`
"""
function AboveSurface(Data_Cart::CartData, DataSurface_Cart::CartData; above=true)

    Data            =   GeoData(ustrip.(Data_Cart.x.val),       ustrip.(Data_Cart.y.val),        ustrip.(Data_Cart.z.val), Data_Cart.fields)
    DataSurface     =   GeoData(ustrip.(DataSurface_Cart.x.val),ustrip.(DataSurface_Cart.y.val), ustrip.(DataSurface_Cart.z.val), DataSurface_Cart.fields )

    return Above    =   AboveSurface(Data, DataSurface; above=above)
end

"""
    Below = BelowSurface(Data_Cart::CartData, DataSurface_Cart::CartData)

Determines if points within the 3D Data_Cart structure are below the Cartesian surface DataSurface_Cart
"""
function BelowSurface(Data_Cart::CartData, DataSurface_Cart::CartData)
    return AboveSurface(Data_Cart::CartData, DataSurface_Cart::CartData; above=false)
end


"""
    VoteMap(DataSets::Vector{GeoData}, criteria::Vector{String}, dims=(50,50,50))

Creates a Vote map which shows consistent features in different 2D/3D tomographic datasets.

The way it works is:
    - Find a common region between the different GeoData sets (overlapping lon/lat/depth regions)
    - Interpolate the fields of all DataSets to common coordinates
    - Filter data points in one model (e.g., areas with a velocity anomaly > 2 percent). Set everything that satisfies this criteria to 1 and everything else to 0.
    - Sum the results of the different datasets

If a feature is consistent between different datasets, it will have larger values. 

# Example
We assume that we have 2 seismic velocity datasets `Data_Zhao_Pwave` and `DataKoulakov_Alps`:
```julia
julia> Data_Zhao_Pwave
GeoData 
  size  : (121, 94, 101)
  lon   ϵ [ 0.0 : 18.0]
  lat   ϵ [ 38.0 : 51.95]
  depth ϵ [ -1001.0 km : -1.0 km]
  fields: (:dVp_Percentage,)
julia> DataKoulakov_Alps
  GeoData 
    size  : (108, 81, 35)
    lon   ϵ [ 4.0 : 20.049999999999997]
    lat   ϵ [ 37.035928143712574 : 49.01197604790419]
    depth ϵ [ -700.0 km : -10.0 km]
    fields: (:dVp_percentage, :dVs_percentage)
```
You can create a VoteMap which combines the two data sets with:
```julia 
julia> Data_VoteMap = VoteMap([Data_Zhao_Pwave,DataKoulakov_Alps],["dVp_Percentage>2.5","dVp_percentage>3.0"])
GeoData 
  size  : (50, 50, 50)
  lon   ϵ [ 4.0 : 18.0]
  lat   ϵ [ 38.0 : 49.01197604790419]
  depth ϵ [ -700.0 km : -10.0 km]
  fields: (:VoteMap,)
```

You can also create a VoteMap of a single dataset:
```julia 
julia> Data_VoteMap = VoteMap(Data_Zhao_Pwave,"dVp_Percentage>2.5", dims=(50,51,52))
GeoData 
  size  : (50, 51, 52)
  lon   ϵ [ 0.0 : 18.0]
  lat   ϵ [ 38.0 : 51.95]
  depth ϵ [ -1001.0 km : -1.0 km]
  fields: (:VoteMap,)
```

"""
function VoteMap(DataSets::Vector{GeoData}, criteria::Vector{String}; dims=(50,50,50))

    numDataSets = length(DataSets)

    if length(criteria) != numDataSets
        error("Need the same number of criteria as the number of data sets")
    end
    
    # Determine the overlapping lon/lat/depth regions of all datasets
    lon_limits  = [minimum(DataSets[1].lon.val);        maximum(DataSets[1].lon.val)];
    lat_limits  = [minimum(DataSets[1].lat.val);        maximum(DataSets[1].lat.val)];
    z_limits    = [minimum(DataSets[1].depth.val);      maximum(DataSets[1].depth.val)];
    for i=1:numDataSets
        lon_limits[1]   =   maximum([lon_limits[1]  minimum(DataSets[i].lon.val)]);
        lon_limits[2]   =   minimum([lon_limits[2]  maximum(DataSets[i].lon.val)]);

        lat_limits[1]   =   maximum([lat_limits[1]  minimum(DataSets[i].lat.val)]);
        lat_limits[2]   =   minimum([lat_limits[2]  maximum(DataSets[i].lat.val)]);
 
        z_limits[1]     =   maximum([z_limits[1]    minimum(DataSets[i].depth.val)]);
        z_limits[2]     =   minimum([z_limits[2]    maximum(DataSets[i].depth.val)]);
    end

    # Loop over all datasets, and interpolate the data set to the new (usually smaller) domain
    VoteMap             =   zeros(Int64,dims)
    for i=1:numDataSets
        VoteMap_Local   =   zeros(Int64,dims)
        
        # Interpolate data set to smaller domain
        DataSet         =   ExtractSubvolume(DataSets[i]; Interpolate=true, Lon_level=lon_limits, Lat_level=lat_limits, Depth_level=z_limits, dims=dims);

        # Extract the criteria to evaluate
        expr            =   Meta.parse(criteria[i]);     # the expression, such as Vs>1.0

        # Extract data field
        if !haskey(DataSet.fields,expr.args[2])
            error("The GeoData set does not have the field: $(expr.args[2])")
        end

        Array3D         =   ustrip.(DataSet.fields[expr.args[2]]);                  # strip units, just in case
        
        # Modify the value, to be Array3D 
        expr_mod        =   Expr(:call, expr.args[1], :($Array3D), expr.args[3]);      # modify the original expression to use Array3D as variable name
        
        # The expression should have a ".", such as Array .> 1.0. If not, it will not apply this in a pointwise manner
        #   Here, we add this dot if it is not there yet
        if cmp(String(expr_mod.args[1])[1],Char('.'))==1
            expr_mod.args[1] = Symbol(".",expr_mod.args[1]);
        end

        ind                 = eval(expr_mod);    # evaluate the modified expression
        VoteMap_Local[ind] .= 1;                 # assign vote-map

        VoteMap = VoteMap + VoteMap_Local;       # Sum 
    end

    DataSet     =   ExtractSubvolume(DataSets[1], Interpolate=true, Lon_level=lon_limits, Lat_level=lat_limits, Depth_level=z_limits, dims=dims);

    # Construct GeoData set that holds the VoteMap (makes it easier to write paraview files)
    VoteData    =   GeoData(DataSet.lon.val,DataSet.lat.val,DataSet.depth.val, (VoteMap=VoteMap,));

    return VoteData
end

# Make this work for single data sets as well
function VoteMap(DataSets::GeoData, criteria::String; dims=(50,50,50))
    VoteMap([DataSets], [criteria]; dims=dims)
end

"""
    RotateTranslateScale!(Data::CartData; Rotate=0, Translate=(0,0,0), Scale=(1.0,1.0,1.0))

Does an in-place rotation, translation and scaling of the Cartesian dataset `Data`. 

# Parameters
Note that we apply the transformations in exactly this order:
-   `Scale`:        scaling applied to the `x,y,z` coordinates of the data set
-   `Rotate`:       rotation around the `x/y` axis (around the center of the box)
-   `Translate`:    translation

# Example
```julia
julia> X,Y,Z   =   XYZGrid(10:20,30:40,-50:-10);
julia> Data_C  =   CartData(X,Y,Z,(Depth=Z,))
CartData 
  size  : (11, 11, 41)
  x     ϵ [ 10.0 : 20.0]
  y     ϵ [ 30.0 : 40.0]
  z     ϵ [ -50.0 : -10.0]
  fields: (:Depth,)
julia> RotateTranslateScale!(Data_C, Rotate=30);
julia> Data_C
CartData 
  size  : (11, 11, 41)
  x     ϵ [ 8.169872981077807 : 21.83012701892219]
  y     ϵ [ 28.16987298107781 : 41.83012701892219]
  z     ϵ [ -50.0 : -10.0]
  fields: (:Depth,)
```
"""
function RotateTranslateScale!(Data::CartData; Rotate=0, Translate=(0,0,0), Scale=(1.0,1.0,1.0))

    X,Y,Z       = Data.x.val,   Data.y.val,     Data.z.val;         # Extract coordinates
    Xr,Yr,Zr    = X,Y,Z;                                            # Rotated coordinates 

    # 1) Scaling
    if length(Scale)==1
        Scale = [Scale Scale Scale];
    end
    Xr .*= Scale[1];
    Yr .*= Scale[2];
    Zr .*= Scale[3];


    # 2) 2D rotation around X/Y axis, around center of box
    Xm,Ym = mean(X), mean(Y);  
    R = [cosd(Rotate[1]) -sind(Rotate[1]); sind(Rotate[1]) cosd(Rotate[1])]; # 2D rotation matrix

    for i in eachindex(X)
        Rot_XY = R*[X[i]-Xm; Y[i]-Ym];
        Xr[i]  = Rot_XY[1] + Xm;
        Yr[i]  = Rot_XY[2] + Ym;
    end 
  
    # 3) Add translation
    Xr .+= Translate[1];
    Yr .+= Translate[2];
    Zr .+= Translate[3];
    
    # Modify original structure
    Data.x.val = Xr;
    Data.y.val = Yr;
    Data.z.val = Zr;
    
end


