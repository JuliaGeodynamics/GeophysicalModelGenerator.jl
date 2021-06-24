# few utils that are useful 

export meshgrid, CrossSection, ExtractSubvolume, SubtractMeanVelocity

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
        
        # we need 3D m[axtrixes for the paraview writing routine to know we are in 3D
        Lon[:,:,1]      =   Lon_p[:,1,:]
        Lat[:,:,1]      =   Lat_p[1,:,:]
        Depth[:,:,1]    =   Depth_p[1,:,:]
        
    end

    if Interpolate
        # Interpolate data on profile
        DataProfile = InterpolateDataFields(V, Lon, Lat, Depth);    
    else
        # extract data (no interplation)
        DataProfile = ExtractDataSets(V, iLon, iLat, iDepth);
    end

    return DataProfile

end


"""
    ExtractSubvolume(V::GeoData; Interpolate=false, Lon_level=nothing, Lat_level=nothing, Depth_level=nothing, dims=(50,50,50))

Extract or "cuts-out" a piece of a 2D or 3D GeoData set, defined by `Lon`, `Lat` and `Depth` coordinates.

This is useful if you are only interested in a part of a much bigger larger data set.

- `Lon_level`,`Lat_level` and `Depth_level` should be tuples that indicate `(minimum_value, maximum_value)` along the respective direction. If not specified 
- By default, `Interpolate=false` and we find the closest indices within the data set
- Alternatively, if `Interpolate=true` we interpolate the data onto a new grid that has dimensions `dims`. This can be useful to compare data sets that are originally given in different resolutions.

# Example:
```julia-repl
julia> Lon,Lat,Depth   =   LonLatDepthGrid(10:20,30:40,(-300:25:0)km);
julia> Data            =   Depth*2;                # some data
julia> Vx,Vy,Vz        =   ustrip(Data*3),ustrip(Data*4),ustrip(Data*5);
julia> Data_set3D      =   GeoData(Lon,Lat,Depth,(Depthdata=Data,LonData=Lon, Velocity=(Vx,Vy,Vz))); 
GeoData 
  size  : (6, 3, 13)
  lon   ϵ [ 10.0 : 15.0]
  lat   ϵ [ 30.0 : 32.0]
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
                                            LinRange(Depth_level[1],    Depth_level[2], dims[3]));

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

function InterpolateDataFields(V::GeoData, Lon, Lat, Depth)

    Lon_vec     =  V.lon.val[:,1,1];
    Lat_vec     =  V.lat.val[1,:,1];
    Depth_vec   =  V.depth.val[1,1,:];

    fields_new  = V.fields;
    field_names = keys(fields_new);
    for i = 1:length(V.fields)
        if typeof(V.fields[i]) <: Tuple
            # vector or anything that contains more than 1 field
            data_tuple = fields_new[i]      # we have a tuple (likely a vector field), so we have to loop 
            data_array = zeros(size(Lon,1),size(Lon,2),size(Lon,3),length(data_tuple));     # create a 3D array that holds the 2D interpolated values
            unit_array = zeros(size(data_array));

            for j=1:length(data_tuple)
                interpol            =   LinearInterpolation((Lon_vec, Lat_vec, Depth_vec), ustrip(data_tuple[j]));      # create interpolation object
                data_array[:,:,:,j] =   interpol.(Lon, Lat, Depth);          
            end
            data_new    = tuple([data_array[:,:,:,c] for c in 1:size(data_array,4)]...)     # transform 3D matrix to tuple

        else
            # scalar field
            interpol    =   LinearInterpolation((Lon_vec, Lat_vec, Depth_vec), V.fields[i]);            # create interpolation object
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

# compute the mean velocity per depth in a 3D dataset and subtract the mean from the given velocities
function SubtractHorizontalMean(V)
    nx        = size(V,1);
    ny        = size(V,2);
    NumLayers = size(V,3); # get the number of depth levels

    V_mean = zeros(size(V))

    for iLayer = 1:NumLayers
        V_sub = tmp[:,:,iLayer] .- mean(filter(!isnan, vec(tmp[:,:,iLayer])));
    end

    return V_sub
end