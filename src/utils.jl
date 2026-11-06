# few utils that are useful

export meshgrid, cross_section, cross_section_volume, cross_section_surface, cross_section_points, extract_subvolume, subtract_horizontalmean
export parse_columns_CSV, votemap, countmap
export interpolate_datafields_2D, interpolate_datafields
export rotate_translate_scale
export lithostatic_pressure!
export flatten_cross_section
export addfield, removefield
export inpoly, inpoly_fast, inpolygon!

using NearestNeighbors

# Lengths in km: a `Quantity` is converted to km, a plain number is taken to be in km already
to_km(x::Quantity) = ustrip(km, x)
to_km(x::Union{Tuple, AbstractArray}) = map(to_km, x)
to_km(x) = x

# Exactly one of the level arguments (or Start) has to be given
function check_single_level(fname; Depth_level, Lat_level, Lon_level, Start)
    return count(!isnothing, (Depth_level, Lat_level, Lon_level, Start)) == 1 || throw(ArgumentError("$fname: give exactly one of Depth_level, Lat_level, Lon_level, Start."))
end

"""
    meshgrid(vx,vy,vz)

Computes an (x,y,z)-grid from the vectors (vx,vy,vz).
For more information, see the MATLAB documentation.
"""
function meshgrid(
        vx::AbstractVector{T}, vy::AbstractVector{T},
        vz::AbstractVector{T}
    ) where {T}
    m, n, o = length(vy), length(vx), length(vz)
    vx = reshape(vx, 1, n, 1)
    vy = reshape(vy, m, 1, 1)
    vz = reshape(vz, 1, 1, o)
    om = ones(Int, m)
    on = ones(Int, n)
    oo = ones(Int, o)
    return (vx[om, :, oo], vy[:, on, oo], vz[om, on, :])
end

"""
    V = addfield(V::AbstractGeneralGrid,field_name::String,data::Any)

Add Fields Data to GeoData or CartData

"""
function addfield(V::AbstractGeneralGrid, field_name::String, data::Any)
    fields_new = V.fields;     new_field = NamedTuple{(Symbol(field_name),)}((data,))
    fields_new = merge(fields_new, new_field)  # replace the field in fields_new

    if isa(V, GeoData)
        V = GeoData(V.lon.val, V.lat.val, V.depth.val, fields_new, V.atts)
    elseif isa(V, CartData)
        V = CartData(V.x.val, V.y.val, V.z.val, fields_new, V.atts)
    else
        error("addfield is only implemented for GeoData and CartData structures")
    end

    return V
end

"""
    V = addfield(V::CartData,new_fields::NamedTuple)

Add `new_fields` fields to a `CartData` dataset
"""
addfield(V::CartData, new_fields::NamedTuple) = CartData(V.x.val, V.y.val, V.z.val, merge(V.fields, new_fields), V.atts)

"""
    V = addfield(V::GeoData,new_fields::NamedTuple)

Add `new_fields` fields to a `GeoData` dataset
"""
addfield(V::GeoData, new_fields::NamedTuple) = GeoData(V.lon.val, V.lat.val, V.depth.val, merge(V.fields, new_fields), V.atts)


"""
    V = addfield(V::Q1Data,new_fields::NamedTuple; cellfield=false)

Add `new_fields` fields to a `Q1Data` dataset; set `cellfield` to `true` if the field is a cell field; otherwise it is a vertex field
"""
function addfield(V::Q1Data, new_fields::NamedTuple; cellfield = false)
    if cellfield
        return Q1Data(V.x.val, V.y.val, V.z.val, V.fields, merge(V.cellfields, new_fields))
    else
        return Q1Data(V.x.val, V.y.val, V.z.val, merge(V.fields, new_fields), V.cellfields)
    end
end


"""
    V = addfield(V::FEData,new_fields::NamedTuple; cellfield=false)

Add `new_fields` fields to a `FEData` dataset; set `cellfield` to `true` if the field is a cell field; otherwise it is a vertex field
"""
function addfield(V::FEData, new_fields::NamedTuple; cellfield = false)
    if cellfield
        return FEData(V.vertices, V.connectivity, V.fields, merge(V.cellfields, new_fields))
    else
        return FEData(V.vertices, V.connectivity, merge(V.fields, new_fields), V.cellfields)
    end
end


# this function is taken from @JeffreySarnoff
function dropnames(namedtuple::NamedTuple, names::Tuple{Vararg{Symbol}})
    keepnames = Base.diff_names(Base._nt_names(namedtuple), names)
    return NamedTuple{keepnames}(namedtuple)
end

"""
    V = removefield(V::AbstractGeneralGrid,field_name::Symbol)

Removes the field with name `field_name` from the GeoData or CartData dataset

"""
function removefield(V::AbstractGeneralGrid, field_name::Symbol)
    fields_new = V.fields
    fields_new = dropnames(fields_new, (field_name,))

    if isa(V, GeoData)
        V = GeoData(V.lon.val, V.lat.val, V.depth.val, fields_new, V.atts)
    elseif isa(V, CartData)
        V = CartData(V.x.val, V.y.val, V.z.val, fields_new, V.atts)
    else
        error("removefield is only implemented for GeoData and CartData structures")
    end

    return V
end

"""
    V = removefield(V::AbstractGeneralGrid,field_name::String)

Removes the field with name `field_name` from the GeoData or CartData dataset

"""
function removefield(V::AbstractGeneralGrid, field_name::String)
    return removefield(V, Symbol(field_name))
end

"""
    V = removefield(V::AbstractGeneralGrid,field_name::NTuple{N,Symbol})

Removes the fields in the tuple `field_name` from the GeoData or CartData dataset

"""
function removefield(V::AbstractGeneralGrid, field_name::NTuple{N, Symbol}) where {N}

    for ifield in 1:N
        V = removefield(V, field_name[ifield])
    end

    return V
end


"""
cross_section_volume(Volume::AbstractGeneralGrid; dims=(100,100), Interpolate=false, Depth_level=nothing; Lat_level=nothing; Lon_level=nothing; Start=nothing, End=nothing, Depth_extent=nothing )

Creates a cross-section through a volumetric (3D) `GeoData` object.

- Cross-sections can be horizontal (map view at a given depth), if `Depth_level` is specified
- They can also be vertical, either by specifying `Lon_level` or `Lat_level` (for a fixed lon/lat), or by defining both `Start=(lon,lat)` & `End=(lon,lat)` points.
- When both `Start=(lon,lat)` & `End=(lon,lat)` are given, one can also provide a the depth extent of the profile by providing Depth_extent=(depth_min,depth_max)
- `Interpolate` indicates whether we want to simply extract the data from the 3D volume (default) or whether we want to linearly interpolate it on a new grid, which has dimensions as specified in `dims`
- `Depth_extent` is an optional parameter that can indicate the depth extent over which you want to interpolate the vertical cross-section. Default is the full vertical extent of the 3D dataset
- `Depth_level` and `Depth_extent` can be given with length units (e.g. `-100km`); plain numbers are in km.

# Example:
```julia-repl
julia> Lon,Lat,Depth   =   lonlatdepth_grid(10:20,30:40,(-300:25:0)km);
julia> Data            =   Depth*2;                # some data
julia> Vx,Vy,Vz        =   ustrip(Data*3),ustrip(Data*4),ustrip(Data*5);
julia> Data_set3D      =   GeoData(Lon,Lat,Depth,(Depthdata=Data,LonData=Lon, Velocity=(Vx,Vy,Vz)));
julia> Data_cross      =   cross_section_volume(Data_set3D, Depth_level=-100km)
GeoData
  size  : (11, 11, 1)
  lon   ϵ [ 10.0 : 20.0]
  lat   ϵ [ 30.0 : 40.0]
  depth ϵ [ -100.0 km : -100.0 km]
  fields: (:Depthdata, :LonData, :Velocity)
```



"""
function cross_section_volume(V::AbstractGeneralGrid; dims = (100, 100), Interpolate = false, Depth_level = nothing, Lat_level = nothing, Lon_level = nothing, Start = nothing, End = nothing, Depth_extent = nothing)

    DataSetType = check_data_set(V)

    if DataSetType != 3
        error("cross_section_volume: the input data set has to be a volume!")
    end

    check_single_level("cross_section_volume"; Depth_level, Lat_level, Lon_level, Start)

    # extract the coordinates
    X, Y, Z = coordinate_grids(V)
    Depth_level, Depth_extent = to_km(Depth_level), to_km(Depth_extent)

    if !isnothing(Depth_level)    # Horizontal slice
        CheckBounds(Z, Depth_level)
        if Interpolate
            Lon, Lat, Depth = lonlatdepth_grid(
                LinRange(minimum(X), maximum(X), dims[1]),
                LinRange(minimum(Y), maximum(Y), dims[2]),
                Depth_level
            )
        else
            ind_z = argmin(abs.(NumValue(Z[1, 1, :]) .- Depth_level))
            iDepth = ind_z:ind_z
            iLon = 1:size(NumValue(X), 1)
            iLat = 1:size(NumValue(Y), 2)
        end
    end

    if !isnothing(Lat_level)   # vertical slice @ given latitude
        CheckBounds(Y, Lat_level)
        if Interpolate
            Lon, Lat, Depth = lonlatdepth_grid(
                LinRange(minimum(X), maximum(X), dims[1]),
                Lat_level,
                LinRange(minimum(Z), maximum(Z), dims[2])
            )
        else
            ind_l = argmin(abs.(Y[1, :, 1] .- Lat_level))
            iDepth = 1:size(Z, 3)
            iLon = 1:size(X, 1)
            iLat = ind_l:ind_l
        end
    end

    if !isnothing(Lon_level)   # vertical slice @ given longitude
        CheckBounds(X, Lon_level)
        if Interpolate
            Lon, Lat, Depth = lonlatdepth_grid(
                Lon_level,
                LinRange(minimum(Y), maximum(Y), dims[1]),
                LinRange(minimum(Z), maximum(Z), dims[2])
            )
        else
            ind_l = argmin(abs.(X[:, 1, 1] .- Lon_level))
            iDepth = 1:size(Z, 3)
            iLat = 1:size(Y, 2)
            iLon = ind_l:ind_l
        end
    end

    # diagonal profile defined by start and end lon/lat points
    if !isnothing(Start)
        if isnothing(End)
            error("Also define End coordinates if you indicate starting lon/lat value")
        end
        Interpolate = true  # we must interpolate in this case

        # if the depth extent is given, modify the Z values to take this into account
        if !isnothing(Depth_extent)
            if length(Depth_extent) != 2
                error("Depth_extent should have length 2")
            end
        end
        if !isnothing(Depth_extent)
            Z = [Depth_extent[1] Depth_extent[2]]
        end

        Lon_dum, Lat_p, Depth_p = lonlatdepth_grid(
            Start[1],
            LinRange(Start[2], End[2], dims[1]),
            LinRange(minimum(Z), maximum(Z), dims[2])
        )

        Lon_p, Lat_dum, Depth = lonlatdepth_grid(
            LinRange(Start[1], End[1], dims[1]),
            Start[2],
            LinRange(minimum(Z), maximum(Z), dims[2])
        )

        Lon = zeros(dims[1], dims[2], 1)
        Lat = zeros(dims[1], dims[2], 1)
        Depth = zeros(dims[1], dims[2], 1) * Depth_p[1]

        # We need 3D matrixes for the paraview writing routine to know we are in 3D
        Lon[:, :, 1] = Lon_p[:, 1, :]
        Lat[:, :, 1] = Lat_p[1, :, :]
        Depth[:, :, 1] = Depth_p[1, :, :]

    end

    if Interpolate
        # Interpolate data on profile
        DataProfile = interpolate_datafields(V, Lon, Lat, NumValue(Depth))
    else
        # extract data (no interpolation)
        DataProfile = ExtractDataSets(V, iLon, iLat, iDepth)
    end

    return DataProfile

end


"""
    cross_section_surface(Surface::AbstractGeneralGrid; dims=(100,), Interpolate=true, Depth_level=nothing, Lat_level=nothing, Lon_level=nothing, Start=nothing, End=nothing)

Creates a cross-section through a surface (2D) `GeoData` object.
The surface is always interpolated onto the profile; `Interpolate=false` throws an error.

- Cross-sections can be horizontal (map view at a given depth), if `Depth_level` is specified
- They can also be vertical, either by specifying `Lon_level` or `Lat_level` (for a fixed lon/lat), or by defining both `Start=(lon,lat)` & `End=(lon,lat)` points. Start and End points will be in km!

- IMPORTANT: The surface to be extracted has to be given as a gridded GeoData object. It may also contain NaNs where it is not defined. Any points lying outside of the defined surface will be considered NaN.

# Example:
```julia-repl
julia> Lon,Lat,Depth   =   lonlatdepth_grid(10:20,30:40,-50km);
julia> Data            =   Depth*2;                # some data
julia> Vx,Vy,Vz        =   ustrip(Data*3),ustrip(Data*4),ustrip(Data*5);
julia> Data_set2D      =   GeoData(Lon,Lat,Depth,(Depth=Depth,));
julia> Data_cross      =   cross_section_surface(Data_set2D, Lat_level =15)
GeoData
  size      : (100,)
  lon       ϵ [ 10.0 : 20.0]
  lat       ϵ [ 15.0 : 15.0]
  depth     ϵ [ NaN : NaN]
  fields    : (:Depth,)
  attributes: ["note"]
```

"""
function cross_section_surface(S::AbstractGeneralGrid; dims = (100,), Interpolate = true, Depth_level = nothing, Lat_level = nothing, Lon_level = nothing, Start = nothing, End = nothing)
    Interpolate || throw(ArgumentError("cross sections through surfaces are always interpolated; Interpolate=false is not supported"))

    DataSetType = check_data_set(S)
    if DataSetType != 2
        error("cross_section_surface: the input data set has to be a surface!")
    end

    check_single_level("cross_section_surface"; Depth_level, Lat_level, Lon_level, Start)

    X, Y, Z = coordinate_grids(S)

    Lon_vec = X[:, 1, 1]
    Lat_vec = Y[1, :, 1]

    if !isnothing(Depth_level)    # horizontal sections require the intersection of two interfaces, which is not supported
        error(" horizontal cross sections not working yet with surface data!")
    end

    if !isnothing(Lat_level)   # vertical slice @ given latitude
        # create a vector that spans the entire dataset @ a given latitude
        Lon = LinRange(minimum(Lon_vec), maximum(Lon_vec), dims[1])
        Lat = ones(size(Lon)) * Lat_level
    end

    if !isnothing(Lon_level)   # vertical slice @ given longitude
        # create a vector that spans the entire dataset @ a given longitude
        Lat = LinRange(minimum(Lat_vec), maximum(Lat_vec), dims[1])
        Lon = ones(size(Lat)) * Lon_level
    end

    # diagonal profile defined by start and end lon/lat points
    if !isnothing(Start)
        if isnothing(End)
            error("Also define End coordinates if you indicate starting lon/lat value")
        end

        Lon = LinRange(Start[1], End[1], dims[1])
        Lat = LinRange(Start[2], End[2], dims[1])

    end

    # now interpolate the depth information of the surface to the profile in question
    interpol = linear_interpolation((Lon_vec, Lat_vec), Z[:, :, 1], extrapolation_bc = NaN)   # create interpolation object, fill with NaNs if outside
    depth_intp = interpol.(Lon, Lat) * km

    # also interpolate any other data that is stored in the GeoData structure on the profile
    fields_new = S.fields
    field_names = keys(fields_new)
    for i in 1:length(S.fields)
        if typeof(S.fields[i]) <: Tuple
            # vector or anything that contains more than 1 field
            data_tuple = fields_new[i]      # we have a tuple (likely a vector field), so we have to loop
            data_array = zeros(size(Lon, 1), size(Lon, 2), length(data_tuple))      # create a 2D array that holds the 2D interpolated values

            for j in 1:length(data_tuple)
                interpol = linear_interpolation((Lon_vec, Lat_vec), dropdims(ustrip.(data_tuple[j]), dims = 3), extrapolation_bc = NaN)       # create interpolation object
                data_array[:, :, j] = interpol.(Lon, Lat)
            end
            data_new = tuple([data_array[:, :, c] for c in 1:size(data_array, 3)]...)     # transform 3D matrix to tuple, do not add unit, as this creates an error in GMG (Issue), to add the unit: *unit(S.fields[i][1][1])

        else
            # scalar field
            interpol = linear_interpolation((Lon_vec, Lat_vec), dropdims(ustrip.(S.fields[i]), dims = 3), extrapolation_bc = NaN)
            data_new = interpol.(Lon, Lat) * unit(S.fields[i][1])                                                  # interpolate data field
        end

        # replace the field
        new_field = NamedTuple{(field_names[i],)}((data_new,))                          # Create a tuple with same name and unit
        fields_new = merge(fields_new, new_field)                                        # replace the field in fields_new

    end

    # create GeoData/CartData structure with the interpolated points
    if isa(S, GeoData)
        Data_profile = GeoData(Lon, Lat, depth_intp, fields_new, S.atts)
    elseif isa(S, CartData)
        Data_profile = CartData(Lon, Lat, depth_intp, fields_new, S.atts)
    else
        error("still to be implemented")
    end
    return Data_profile
end


"""
    cross_section_points(P::GeoData; Depth_level=nothing, Lat_level=nothing, Lon_level=nothing, Start=nothing, End=nothing, section_width=10km)

Creates a projection of separate points (saved as a GeoData object) onto a chosen plane. Only points within a band of width `section_width` centered on the plane are taken into account (for `Start`/`End`, also only points between them).
Exactly one of `Depth_level`, `Lat_level`, `Lon_level`, or `Start` & `End` has to be given. Returns `nothing` if no point is selected.
`Depth_level` and `section_width` can be given with length units (e.g. `10km` or `10000m`); plain numbers are in km.

"""
function cross_section_points(P::GeoData; Depth_level = nothing, Lat_level = nothing, Lon_level = nothing, Start = nothing, End = nothing, section_width = 10km)

    DataSetType = check_data_set(P)
    if DataSetType != 1
        error("cross_section_points: the input data set has to be a pointwise data set!")
    end

    check_single_level("cross_section_points"; Depth_level, Lat_level, Lon_level, Start)
    Depth_level, section_width = to_km(Depth_level), to_km(section_width)
    if !isnothing(Depth_level)
        ind = findall(-0.5 * section_width .< (NumValue(P.depth) .- Depth_level) .< 0.5 * section_width) # find all points around the desired depth level

        # create temporary variables
        lon_tmp = NumValue(P.lon.val[ind])
        lat_tmp = NumValue(P.lat.val[ind])
        depth_tmp = NumValue(P.depth.val[ind])
        depth_proj = ones(size(depth_tmp)) * Depth_level

        # create fields that will be stored additionally on the GeoData structure
        field_tmp = (depth_proj = depth_proj, lat_proj = lat_tmp, lon_proj = lon_tmp) # these are the projected points
    end

    if !isnothing(Lat_level)   # vertical slice @ given latitude
        # distance (km) of each point from the parallel, on a sphere with the mean Earth radius of 6371 km
        ind = findall(abs.(NumValue(P.lat) .- Lat_level) .* (π / 180 * 6371.0) .< 0.5 * section_width)

        # create temporary variables
        lon_tmp = NumValue(P.lon.val[ind])
        lat_tmp = NumValue(P.lat.val[ind])
        depth_tmp = NumValue(P.depth.val[ind])
        lat_proj = ones(size(depth_tmp)) * Lat_level

        # data to be stored on the new GeoData structure
        field_tmp = (depth_proj = depth_tmp, lat_proj = lat_proj, lon_proj = lon_tmp) # these are the projected points

    end

    if !isnothing(Lon_level)   # vertical slice @ given longitude
        # distance (km) of each point from the meridian, on a sphere with the mean Earth radius of 6371 km
        # (longitude difference wrapped to -180..180°, so 0..360° data work with -180..180° levels and vice versa)
        ind = findall(abs.(mod.(NumValue(P.lon) .- Lon_level .+ 180, 360) .- 180) .* cosd.(NumValue(P.lat)) .* (π / 180 * 6371.0) .< 0.5 * section_width)

        # create temporary variables
        lon_tmp = NumValue(P.lon.val[ind])
        lat_tmp = NumValue(P.lat.val[ind])
        depth_tmp = NumValue(P.depth.val[ind])
        lon_proj = ones(size(depth_tmp)) * match_lon_convention(Lon_level, NumValue(P.lon))

        # create fields that will be stored on the GeoData structure
        field_tmp = (depth_proj = depth_tmp, lat_proj = lat_tmp, lon_proj = lon_proj) # these are the projected points

    end

    # vertical profile defined by start and end lon/lat points
    # here we need to compute the distance to a distance_to_plane
    # also, we need to project the points on the profile plane for later plotting
    if !isnothing(Start)
        if isnothing(End)
            error("Also define End coordinates if you indicate starting lon/lat value")
        end
        Start == End && throw(ArgumentError("cross_section_points: Start and End must differ to define a profile, got Start = End = $Start"))

        # choose projection point based on Start and End coordinates of the profile
        p_Point = ProjectionPoint(Lat = 0.5 * (Start[2] + End[2]), Lon = 0.5 * (Start[1] + End[1]))

        # convert P to UTM Data
        P_UTM = convert2UTMzone(P, p_Point) # convert to UTM

        # create a GeoData set containing the points that create the profile plane (we need three points to uniquely define that plane)
        # here, we define the points in a way that the angle between P1-P2 and P1-P3 vectors is 90° --> useful for the cross product
        Profile = GeoData([Start[1] Start[1]  End[1]], [Start[2]  Start[2] End[2]], [0 -200 0] * km, (depth = [0 -200 0] * km,))
        Profile_UTM = convert2UTMzone(Profile, p_Point) # convert to UTM

        # compute the unit normal of the profile plane using the cross product
        # UTM coordinates and depth are all in m
        a1 = Profile_UTM.EW.val[2] - Profile_UTM.EW.val[1]
        a2 = Profile_UTM.NS.val[2] - Profile_UTM.NS.val[1]
        a3 = Profile_UTM.depth.val[2] - Profile_UTM.depth.val[1]

        b1 = Profile_UTM.EW.val[3] - Profile_UTM.EW.val[1]
        b2 = Profile_UTM.NS.val[3] - Profile_UTM.NS.val[1]
        b3 = Profile_UTM.depth.val[3] - Profile_UTM.depth.val[1]

        nx = a2 * b3 - a3 * b2
        ny = a3 * b1 - a1 * b3
        nz = a1 * b2 - a2 * b1

        t = (nx * Profile_UTM.EW.val[1] .- nx * P_UTM.EW.val .+ ny * Profile_UTM.NS.val[1] .- ny * P_UTM.NS.val .+ nz * Profile_UTM.depth.val[1] .- nz * P_UTM.depth.val) / (nx * nx + ny * ny + nz * nz)

        # compute the distance to the plane
        dist = sqrt.((t .* nx) .^ 2 + (t .* ny) .^ 2 + (t .* nz) .^ 2)

        # find the points that are within the required window around the profile and between Start and End
        # (t_along: 0 at Start, 1 at End; distance is in m, section_width in km)
        along = [Profile_UTM.EW.val[3] - Profile_UTM.EW.val[1], Profile_UTM.NS.val[3] - Profile_UTM.NS.val[1]]
        t_along = ((P_UTM.EW.val .- Profile_UTM.EW.val[1]) .* along[1] .+ (P_UTM.NS.val .- Profile_UTM.NS.val[1]) .* along[2]) ./ sum(abs2, along)
        ind = findall((dist .< 500 * section_width) .& (0 .<= t_along .<= 1))

        # project the points on the plane (only the relevant ones)
        px = P_UTM.EW.val[ind] + t[ind] .* nx
        py = P_UTM.NS.val[ind] + t[ind] .* ny
        pz = P_UTM.depth.val[ind] + t[ind] .* nz

        # the projected points are given in UTM coordinates and not in lon/lat/depth
        # therefore we have to recompute the lat/lon/depth values of the projected points
        # then we will return a GeoData structure with all information included
        trans = LLAfromUTM(p_Point.zone, p_Point.isnorth, wgs84) # set up transformation

        plon = zeros(size(ind))
        plat = zeros(size(ind))
        pdepth = zeros(size(ind))

        for i in eachindex(ind)
            utmi = UTM(px[i], py[i], pz[i])
            llai = trans(utmi)

            plon[i] = llai.lon
            plat[i] = llai.lat
            pdepth[i] = llai.alt
        end

        # data to be stored in the GeoData structure
        field_tmp = (depth_proj = pdepth / 1.0e3, lat_proj = plat, lon_proj = plon) # these are the projected points
    end

    # also transfer any other data that is stored in the GeoData structure
    fields_new = P.fields
    field_names = keys(fields_new)
    for i in 1:length(P.fields)
        if typeof(P.fields[i]) <: Tuple
            # vector or anything that contains more than 1 field
            data_tuple = fields_new[i]      # we have a tuple (likely a vector field), so we have to loop
            data_new = Tuple(ustrip.(d[ind]) for d in data_tuple)   # units are dropped, as they cause errors in GMG

        else
            # scalar field
            data_new = fields_new[i][ind]                                                  # interpolate data field
        end

        # replace the field
        new_field = NamedTuple{(field_names[i],)}((data_new,))                          # Create a tuple with same name and unit
        fields_new = merge(fields_new, new_field)                                        # replace the field in fields_new

    end

    # merge old and new fields
    fields_new = merge(fields_new, field_tmp)

    # create a GeoData structure to return
    if length(ind) > 0
        Data_profile = GeoData(P.lon.val[ind], P.lat.val[ind], P.depth.val[ind], fields_new, P.atts)
    else
        Data_profile = nothing
    end


    return Data_profile
end

"""
    cross_section(DataSet::AbstractGeneralGrid; dims=(100,100), Interpolate=nothing, Depth_level=nothing, Lat_level=nothing, Lon_level=nothing, Start=nothing, End=nothing, Depth_extent=nothing, section_width=50km)

Creates a cross-section through a `GeoData` object.

- Cross-sections can be horizontal (map view at a given depth), if `Depth_level` is specified
- They can also be vertical, either by specifying `Lon_level` or `Lat_level` (for a fixed lon/lat), or by defining both `Start=(lon,lat)` & `End=(lon,lat)` points.
- Depending on the type of input data (volume, surface or point data), cross sections will be created in a different manner:
1. Volume data: data will be interpolated or directly extracted from the data set.
2. Surface data: surface data will be interpolated
3. Point data: data will be projected to the chosen profile. Only data within a chosen distance (default is 50 km) will be used

- `Interpolate` indicates whether we want to simply extract the data from a volumetric data set (default) or whether we want to linearly interpolate it on a new grid, which has dimensions as specified in `dims`.
  Surface data are always interpolated (`Interpolate=false` throws an error); point data ignore `Interpolate`.
- `section_width` is the width of the band around the profile within which point data are projected onto it. It can be given with length units; a plain number is in km.
- `Depth_level` and `Depth_extent` can be given with length units; plain numbers are in km.

# Example:
```julia-repl
julia> Lon,Lat,Depth   =   lonlatdepth_grid(10:20,30:40,(-300:25:0)km);
julia> Data            =   Depth*2;                # some data
julia> Vx,Vy,Vz        =   ustrip(Data*3),ustrip(Data*4),ustrip(Data*5);
julia> Data_set3D      =   GeoData(Lon,Lat,Depth,(Depthdata=Data,LonData=Lon, Velocity=(Vx,Vy,Vz)));
julia> Data_cross      =   cross_section(Data_set3D, Depth_level=-100km)
GeoData
  size  : (11, 11, 1)
  lon   ϵ [ 10.0 : 20.0]
  lat   ϵ [ 30.0 : 40.0]
  depth ϵ [ -100.0 km : -100.0 km]
  fields: (:Depthdata, :LonData, :Velocity)
```

"""
function cross_section(DataSet::AbstractGeneralGrid; dims = (100, 100), Interpolate = nothing, Depth_level = nothing, Lat_level = nothing, Lon_level = nothing, Start = nothing, End = nothing, Depth_extent = nothing, section_width = 50km)

    DataSetType = check_data_set(DataSet)  # check which kind of data set we are dealing with

    if DataSetType == 1 # points
        DataProfile = cross_section_points(DataSet; Depth_level, Lat_level, Lon_level, Start, End, section_width)
    elseif DataSetType == 2 # surface
        DataProfile = cross_section_surface(DataSet; dims, Interpolate = something(Interpolate, true), Depth_level, Lat_level, Lon_level, Start, End)
    elseif DataSetType == 3 # volume
        DataProfile = cross_section_volume(DataSet; dims, Interpolate = something(Interpolate, false), Depth_level, Lat_level, Lon_level, Start, End, Depth_extent)

        # add field that has coordinates along the profile
        DataProfile = addfield(DataProfile, "FlatCrossSection", flatten_cross_section(DataProfile))
    end

    return DataProfile
end

"""
    flatten_cross_section(V::CartData)
Takes a diagonal 3D cross_section and flattens it to be converted to a 2D Grid by create_CartGrid
# Example
```julia
Grid                    = create_CartGrid(size=(100,100,100), x=(0.0km, 99.9km), y=(-10.0km, 20.0km), z=(-40km,4km));
X,Y,Z                   = xyz_grid(Grid.coord1D...);
DataSet                 = CartData(X,Y,Z,(Depthdata=Z,));

Data_Cross              = cross_section(DataSet, dims=(100,100), Interpolate=true, Start=(ustrip(Grid.min[1]),ustrip(Grid.max[2])), End=(ustrip(Grid.max[1]), ustrip(Grid.min[2])))

x_new = flatten_cross_section(Data_Cross)

# This flattened cross_section can be added to original Data_Cross by addfield()

Data_Cross = addfield(Data_Cross,"FlatCrossSection", x_new)
CartData
    size    : (100, 100, 1)
    x       ϵ [ 0.0 : 99.9]
    y       ϵ [ -10.0 : 20.0]
    z       ϵ [ -40.0 : 4.0]
    fields  : (:Depthdata, :FlatCrossSection)
  attributes: ["note"]

```
"""
function flatten_cross_section(V::CartData)

    x_new = sqrt.((V.x.val .- V.x.val[1, 1, 1]) .^ 2 .+ (V.y.val .- V.y.val[1, 1, 1]) .^ 2) # NOTE: the result is in km, as V.x and V.y are stored in km

    return x_new

end

"""
    flatten_cross_section(V::GeoData)
    This function takes a 3D cross section through a GeoData structure and computes the distance along the cross section for later 2D processing/plotting
    ```julia-repl
    julia> Lon,Lat,Depth   =   lonlatdepth_grid(10:20,30:40,(-300:25:0)km);
    julia> Data            =   Depth*2;                # some data
    julia> Vx,Vy,Vz        =   ustrip(Data*3),ustrip(Data*4),ustrip(Data*5);
    julia> Data_set3D      =   GeoData(Lon,Lat,Depth,(Depthdata=Data,LonData=Lon, Velocity=(Vx,Vy,Vz)));
    julia> Data_cross      =   cross_section(Data_set3D, Start=(10,30),End=(20,40))
    julia> x_profile        =   flatten_cross_section(Data_cross)
    julia> Data_cross      =   addfield(Data_cross,"x_profile",x_profile)

    ```
"""
function flatten_cross_section(V::GeoData; Start = nothing)

    if isnothing(Start)
        lla_start = LLA(V.lat.val[1][1][1], V.lon.val[1][1][1], 0.0) # start point, at the surface
    else
        lla_start = LLA(Start[2], Start[1], 0.0)
    end
    x_new = zeros(size(V.lon))

    for i in eachindex(x_new)
        x_new[i] = euclidean_distance(LLA(V.lat.val[i], V.lon.val[i], 0.0), lla_start) / 1.0e3 # compute distance as if points were at the surface, CONVERTED TO KM !!!
    end

    return x_new
end

"""
    extract_subvolume(V::GeoData; Interpolate=false, Lon_level=nothing, Lat_level=nothing, Depth_level=nothing, dims=(50,50,50))

Extract or "cuts-out" a piece of a 2D or 3D GeoData set, defined by `Lon`, `Lat` and `Depth` coordinates.

This is useful if you are only interested in a part of a much bigger larger data set.

- `Lon_level`,`Lat_level` and `Depth_level` should be tuples that indicate `(minimum_value, maximum_value)` along the respective direction. If not specified we use the full range.
- `Depth_level` can be given with length units (e.g. `(-200km, -100km)`); plain numbers are in km.
- By default, `Interpolate=false` and we find the closest indices within the data set (so your new data set will not go exactly from minimum to maximum).
- Alternatively, if `Interpolate=true` we interpolate the data onto a new grid that has dimensions `dims`. This can be useful to compare data sets that are originally given in different resolutions.

# 3D Example with no interpolation:
```julia-repl
julia> Lon,Lat,Depth   =   lonlatdepth_grid(10:20,30:40,(-300:25:0)km);
julia> Data            =   Depth*2;                # some data
julia> Vx,Vy,Vz        =   ustrip(Data*3),ustrip(Data*4),ustrip(Data*5);
julia> Data_set3D      =   GeoData(Lon,Lat,Depth,(Depthdata=Data,LonData=Lon, Velocity=(Vx,Vy,Vz)))
GeoData
  size  : (11, 11, 13)
  lon   ϵ [ 10.0 : 20.0]
  lat   ϵ [ 30.0 : 40.0]
  depth ϵ [ -300.0 km : 0.0 km]
  fields: (:Depthdata, :LonData, :Velocity)
julia> Data_extracted = extract_subvolume(Data_set3D,Lon_level=(10,12),Lat_level=(35,40))
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
julia> Data_extracted = extract_subvolume(Data_set3D,Lon_level=(10,12),Lat_level=(35,40), Interpolate=true, dims=(50,51,52))
GeoData
  size  : (50, 51, 52)
  lon   ϵ [ 10.0 : 12.0]
  lat   ϵ [ 35.0 : 40.0]
  depth ϵ [ -300.0 km : 0.0 km]
  fields: (:Depthdata, :LonData, :Velocity)
```

"""
function extract_subvolume(V::GeoData; Interpolate = false, Lon_level = nothing, Lat_level = nothing, Depth_level = nothing, dims = (50, 50, 50))

    if isnothing(Lon_level)
        Lon_level = (minimum(V.lon.val), maximum(V.lon.val))
    end
    if isnothing(Lat_level)
        Lat_level = (minimum(V.lat.val), maximum(V.lat.val))
    end
    if isnothing(Depth_level)
        Depth_level = (minimum(V.depth.val), maximum(V.depth.val))
    end
    Depth_level = to_km(Depth_level)
    (Lon_level[1] <= Lon_level[2] && Lat_level[1] <= Lat_level[2]) || throw(ArgumentError("Lon_level and Lat_level have to be ordered (min, max)."))
    if Interpolate
        Lon, Lat, Depth = lonlatdepth_grid(
            LinRange(Lon_level[1], Lon_level[2], dims[1]),
            LinRange(Lat_level[1], Lat_level[2], dims[2]),
            LinRange(Depth_level[1], Depth_level[2], dims[3])
        )
        Data_extract = interpolate_datafields(V, Lon, Lat, Depth)

    else
        # Don't interpolate
        i_s, i_e = argmin(abs.(V.lon.val[:, 1, 1] .- Lon_level[1])), argmin(abs.(V.lon.val[:, 1, 1] .- Lon_level[2]))
        iLon = i_s:i_e

        i_s, i_e = argmin(abs.(V.lat.val[1, :, 1] .- Lat_level[1])), argmin(abs.(V.lat.val[1, :, 1] .- Lat_level[2]))
        iLat = i_s:i_e

        i_s, i_e = argmin(abs.(V.depth.val[1, 1, :] .- Depth_level[1])), argmin(abs.(V.depth.val[1, 1, :] .- Depth_level[2]))
        step = 1
        if i_e < i_s
            step = -1
        end
        iDepth = i_s:step:i_e
        Data_extract = ExtractDataSets(V, iLon, iLat, iDepth)
    end

    return Data_extract
end


"""
    extract_subvolume(V::CartData; Interpolate=true, X_level=nothing, Y_level=nothing, Z_level=nothing, dims=(50,50,50))

Extract or "cuts-out" a piece of a 2D or 3D CartData set, defined by `X`, `Y` and `Z` coordinates.

This is useful if you are only interested in a part of a much bigger larger data set.

- `X_level`,`Y_level` and `Z_level` should be tuples that indicate `(minimum_value, maximum_value)` along the respective direction. If not specified we use the full range.
  They can be given with length units (e.g. `(-20km, -5km)`); plain numbers are in km.
- By default, `Interpolate=true` and we interpolate the data onto a new grid that has dimensions `dims`. This can be useful to compare data sets that are originally given in different resolutions.
- Alternatively, if `Interpolate=false` we find the closest indices within the data set (so your new data set will not go exactly from minimum to maximum).
- For a vertical cross-section computed with `cross_section(V, Start=.., End=..)`, `X_level` is the distance along the profile and the data are always interpolated.

# 3D Example with no interpolation:
```julia-repl
julia> Lon,Lat,Depth   =   lonlatdepth_grid(10:20,30:40,(-300:25:0)km);
julia> Data            =   Depth*2;                # some data
julia> Vx,Vy,Vz        =   ustrip(Data*3),ustrip(Data*4),ustrip(Data*5);
julia> Data_set3D      =   GeoData(Lon,Lat,Depth,(Depthdata=Data,LonData=Lon, Velocity=(Vx,Vy,Vz)))
GeoData
  size  : (11, 11, 13)
  lon   ϵ [ 10.0 : 20.0]
  lat   ϵ [ 30.0 : 40.0]
  depth ϵ [ -300.0 km : 0.0 km]
  fields: (:Depthdata, :LonData, :Velocity)
julia> Data_extracted = extract_subvolume(Data_set3D,Lon_level=(10,12),Lat_level=(35,40))
GeoData
  size  : (3, 6, 13)
  lon   ϵ [ 10.0 : 12.0]
  lat   ϵ [ 35.0 : 40.0]
  depth ϵ [ -300.0 km : 0.0 km]
  fields: (:Depthdata, :LonData, :Velocity)
```
By default it extracts the data points closest to the area defined by Lon_level/Lat_level/Depth_level.


# 2D Example along a cross-section through 3D data:
```julia-repl
julia> X,Y,Z = xyz_grid(10:20,30:40,-300:25:0);
julia> Data = Z.*2
julia> Data_Int = Int64.(Data)
julia> DataSet_Cart = CartData(X,Y,Z,(Data=Data,Data_Int=Data_Int, Velocity=(X,Y,Z)))

julia> Data_cross = cross_section(DataSet_Cart, Start=(11.0,35), End=(19, 39.0))
CartData
    size    : (100, 100, 1)
    x       ϵ [ 11.0 : 19.0]
    y       ϵ [ 35.0 : 39.0]
    z       ϵ [ -300.0 : 0.0]
    fields  : (:Data, :Data_Int, :Velocity, :FlatCrossSection)
  attributes: ["note"]

julia> Data_extracted = extract_subvolume(Data_cross, X_level=(1,7), Z_level=(-200,-100))
  CartData
      size    : (50, 50, 1)
      x       ϵ [ 11.894427190999917 : 17.260990336999413]
      y       ϵ [ 35.44721359549995 : 38.130495168499706]
      z       ϵ [ -200.0 : -100.0]
      fields  : (:FlatCrossSection, :Data, :Data_Int, :Velocity)
    attributes: ["note"]
julia> typeof(Data_extracted.fields.Data_Int)
    Array{Int64, 3}
```

"""
function extract_subvolume(
        V::CartData;
        Interpolate = true,
        X_level = nothing,
        Y_level = nothing,
        Z_level = nothing,
        dims = (50, 50, 50)
    )

    is_cross_section = haskey(V.fields, :FlatCrossSection)
    X_level, Y_level, Z_level = to_km(X_level), to_km(Y_level), to_km(Z_level)
    is_cross_section && !isnothing(Y_level) && throw(ArgumentError("Y_level cannot be used on a cross section."))
    if isnothing(Z_level)
        Z_level = (minimum(V.z.val), maximum(V.z.val))
    end

    if !is_cross_section
        if isnothing(X_level)
            X_level = (minimum(V.x.val), maximum(V.x.val))
        end
        if isnothing(Y_level)
            Y_level = (minimum(V.y.val), maximum(V.y.val))
        end
    end

    if Interpolate == true && !is_cross_section
        X, Y, Z = lonlatdepth_grid(
            LinRange(X_level[1], X_level[2], dims[1]),
            LinRange(Y_level[1], Y_level[2], dims[2]),
            LinRange(Z_level[1], Z_level[2], dims[3])
        )

        Data_extract = interpolate_datafields(V, X, Y, Z)

    elseif is_cross_section

        # we are dealing with a vertical cross-section through a 3D dataset computed with cross_section(V,Start=.., End=...)
        Xcross = V.fields.FlatCrossSection
        if isnothing(X_level)
            X_level = extrema(Xcross)
        end

        dims_cross = (dims[1], dims[2], 1)

        # we need to interpolate the data onto a new grid given by X_level and Z_level
        X_level_cross = X_level

        interpol_x = linear_interpolation(Xcross[:, 1, 1], V.x.val[:, 1, 1], extrapolation_bc = NaN)       # create interpolation object
        interpol_y = linear_interpolation(Xcross[:, 1, 1], V.y.val[:, 1, 1], extrapolation_bc = NaN)       # create interpolation object

        x = LinRange(X_level_cross[1], X_level_cross[2], dims_cross[1])
        z = LinRange(Z_level[1], Z_level[2], dims_cross[2])

        X, Y, Z = zeros(dims[1], dims[2], 1), zeros(dims[1], dims[2], 1), zeros(dims[1], dims[2], 1)
        X_cross = zero(X)
        for (i, x_val) in enumerate(x), (j, z_val) in enumerate(z)
            X[i, j, 1] = interpol_x(x_val)
            Y[i, j, 1] = interpol_y(x_val)
            Z[i, j, 1] = z_val
            X_cross[i, j, 1] = x_val
        end

        Data_extract = interpolate_data_fields_cross_section(V, X, Y, Z, X_cross)

    else
        # Don't interpolate
        i_s, i_e = argmin(abs.(V.x.val[:, 1, 1] .- X_level[1])), argmin(abs.(V.x.val[:, 1, 1] .- X_level[2]))
        iLon = i_s:i_e

        i_s, i_e = argmin(abs.(V.y.val[1, :, 1] .- Y_level[1])), argmin(abs.(V.y.val[1, :, 1] .- Y_level[2]))
        iLat = i_s:i_e

        i_s, i_e = argmin(abs.(V.z.val[1, 1, :] .- Z_level[1])), argmin(abs.(V.z.val[1, 1, :] .- Z_level[2]))
        step = 1
        if i_e < i_s
            step = -1
        end
        iDepth = i_s:step:i_e
        Data_extract = ExtractDataSets(V, iLon, iLat, iDepth)
    end

    return Data_extract
end


"""
    interpolate_data_fields_cross_section(V::CartData, X,Y,Z,Xcross)

Interpolates data fields along a cross-section defined by `Xcross` and `Z`
"""
function interpolate_data_fields_cross_section(V::CartData, X, Y, Z, X_cross)

    Data_extract = CartData(X, Y, Z, (FlatCrossSection = X_cross,), V.atts)

    fields_new = V.fields
    field_names = keys(fields_new)
    for i in 1:length(V.fields)
        if typeof(V.fields[i]) <: Tuple
            # vector or anything that contains more than 1 field
            data_tuple = fields_new[i]      # we have a tuple (likely a vector field), so we have to loop
            data_array = zeros(size(Data_extract.x)..., length(data_tuple))      # create a 3D array that holds the 2D interpolated values

            for j in 1:length(data_tuple)
                interpol = linear_interpolation((V.fields.FlatCrossSection[:, 1, 1], V.z.val[1, :, 1]), ustrip.(data_tuple[j][:, :, 1]), extrapolation_bc = NaN)       # create interpolation object
                data_array[:, :, :, j] = interpol.(X_cross, Z)
            end
            data_new = tuple([data_array[:, :, :, c] for c in 1:size(data_array, 4)]...)     # transform 3D matrix to tuple

        else
            # scalar field
            if isa(V.fields[i][1], Int64)
                # categorical data: nearest neighbor, held constant outside the domain (an Int cannot be NaN); see #20, #191
                interpol = constant_interpolation((V.fields.FlatCrossSection[:, 1, 1], V.z.val[1, :, 1]), V.fields[i][:, :, 1], extrapolation_bc = Flat())
                data_new = Int64.(interpol.(X_cross, Z))
            else
                interpol = linear_interpolation((V.fields.FlatCrossSection[:, 1, 1], V.z.val[1, :, 1]), V.fields[i][:, :, 1], extrapolation_bc = NaN)             # create interpolation object
                data_new = interpol.(X_cross, Z)                                                  # interpolate data field
            end
        end
        Data_extract = addfield(Data_extract, String(field_names[i]), data_new)
    end

    return Data_extract

end


function CheckBounds(Data::AbstractArray, Data_Cross)

    min_Data, max_Data = NumValue(minimum(Data)), NumValue(maximum(Data))
    return if ustrip(Data_Cross) < min_Data || ustrip(Data_Cross) > max_Data
        error("Outside bounds [$min_Data : $max_Data]; $Data_Cross")
    end
end

# CHECKS FOR VOLUME, SURFACE OR POINTS
function check_data_set(DataSet::GeoData)
    if length(size(DataSet.lon)) == 1 # scattered points
        return 1
    else
        if any(size(DataSet.lon) .== 1) # surface data
            return 2
        else # volume data
            return 3
        end
    end
end

function check_data_set(DataSet::CartData)
    if length(size(DataSet.x)) == 1 # scattered points
        return 1
    else
        if any(size(DataSet.x) .== 1) # surface data
            return 2
        else # volume data
            return 3
        end
    end
end


"""
    Data_interp = interpolate_datafields(V::AbstractGeneralGrid, Lon, Lat, Depth)

Interpolates a data field `V` (`GeoData`, `CartData` or `UTMData`) on a grid defined by `Lon,Lat,Depth`, given in the coordinates of `V` (`EW,NS,Depth` for `UTMData`).
Floating-point fields are interpolated linearly and are `NaN` outside `V`; `Int64` fields (such as phases) take the value of the nearest point.

# Example
```julia
julia> x        =   0:2:10
julia> y        =   -5:5
julia> z        =   -10:2:2
julia> X,Y,Z    =   xyz_grid(x, y, z);
julia> Data     =   Z
julia> Data_set1=   CartData(X,Y,Z, (FakeData=Data,Data2=Data.+1.))
CartData
    size    : (6, 11, 7)
    x       ϵ [ 0.0 km : 10.0 km]
    y       ϵ [ -5.0 km : 5.0 km]
    z       ϵ [ -10.0 km : 2.0 km]
    fields  : (:FakeData, :Data2)
  attributes: ["note"]

julia> X,Y,Z    =   xyz_grid(0:4:10, -1:.1:1, -5:.1:1 );
julia> Data_set2= interpolate_datafields(Data_set1, X,Y,Z)
```

"""
function interpolate_datafields(V::AbstractGeneralGrid, Lon, Lat, Depth)

    X, Y, Z = coordinate_grids(V)

    # evenly spaced axes as ranges, so that a lookup computes the enclosing cell instead of searching for it
    Lon_vec = as_range(NumValue(X[:, 1, 1]))
    Lat_vec = as_range(NumValue(Y[1, :, 1]))
    Depth_vec = as_range(Z[1, 1, :])
    if Depth_vec[1] > Depth_vec[end]
        ReverseData = true
    else
        ReverseData = false
    end

    fields_new = V.fields
    field_names = keys(fields_new)
    for i in 1:length(V.fields)
        if typeof(V.fields[i]) <: Tuple
            # vector or anything that contains more than 1 field
            data_tuple = fields_new[i]      # we have a tuple (likely a vector field), so we have to loop
            data_array = zeros(size(Lon, 1), size(Lon, 2), size(Lon, 3), length(data_tuple))      # create a 3D array that holds the 2D interpolated values

            for j in 1:length(data_tuple)
                if ReverseData
                    ndim = length(size(data_tuple[j]))
                    interpol = linear_interpolation((Lon_vec, Lat_vec, reverse(Depth_vec)), reverse(ustrip.(data_tuple[j]), dims = ndim), extrapolation_bc = NaN)       # create interpolation object
                else
                    interpol = linear_interpolation((Lon_vec, Lat_vec, Depth_vec), ustrip.(data_tuple[j]), extrapolation_bc = NaN)       # create interpolation object
                end
                data_array[:, :, :, j] = interpol.(Lon, Lat, ustrip.(Depth))
            end
            data_new = tuple([data_array[:, :, :, c] for c in 1:size(data_array, 4)]...)     # transform 3D matrix to tuple

        else
            # scalar field
            # categorical data: nearest neighbor, held constant outside the domain (an Int cannot be NaN); see #20, #191
            isint = isa(V.fields[i][1], Int64)
            make_interpolation = isint ? constant_interpolation : linear_interpolation
            bc = isint ? Flat() : NaN
            if ReverseData
                ndim = length(size(V.fields[i]))
                interpol = make_interpolation((Lon_vec, Lat_vec, reverse(Depth_vec)), reverse(V.fields[i], dims = ndim), extrapolation_bc = bc)             # create interpolation object
            else
                interpol = make_interpolation((Lon_vec, Lat_vec, Depth_vec), V.fields[i], extrapolation_bc = bc)             # create interpolation object
            end
            data_new = interpol.(Lon, Lat, ustrip.(Depth))                                                  # interpolate data field
            isint && (data_new = Int64.(data_new))
        end

        # replace the one
        new_field = NamedTuple{(field_names[i],)}((data_new,))                          # Create a tuple with same name
        fields_new = merge(fields_new, new_field)                                        # replace the field in fields_new

    end


    # Create a struct of the same type as V with the newly interpolated fields
    if isa(V, GeoData)
        Data_profile = GeoData(Lon, Lat, Depth, fields_new, V.atts)
    elseif isa(V, CartData)
        Data_profile = CartData(Lon, Lat, Depth, fields_new, V.atts)
    elseif isa(V, UTMData)
        Data_profile = UTMData(Lon, Lat, Depth, V.zone[1], V.northern[1], fields_new, V.atts)
    else
        error("interpolate_datafields is not implemented for $(typeof(V))")
    end

    return Data_profile
end

"""
    interpolate_datafields_2D(V::GeoData, Lon, Lat)

Interpolates a data field `V` on a 2D grid defined by `Lon,Lat`. Typically used for horizontal surfaces.
Returns the interpolated depth and fields; points outside `V` take the value at the nearest edge of `V`.
"""
function interpolate_datafields_2D(V::GeoData, Lon, Lat)
    return InterpolateDataFields2D_vecs(V.lon.val[:, 1, 1], V.lat.val[1, :, 1], V.depth, V.fields, Lon, Lat)
end


"""
    interpolate_datafields_2D(V::UTMData, EW, NS)

Interpolates a data field `V` on a 2D grid defined by `UTM`. Typically used for horizontal surfaces.
Returns the interpolated depth and fields; points outside `V` take the value at the nearest edge of `V`.
"""
function interpolate_datafields_2D(V::UTMData, EW, NS)
    EW_vec = V.EW.val[:, 1, 1]
    NS_vec = V.NS.val[1, :, 1]
    return InterpolateDataFields2D_vecs(EW_vec, NS_vec, V.depth, V.fields, EW, NS)
end

"""
    interpolate_datafields_2D(V::CartData, X, Y)

Interpolates a data field `V` on a 2D CartData grid defined by `X`,`Y`. Typically used for horizontal surfaces.
Returns the interpolated depth and fields; points outside `V` take the value at the nearest edge of `V`.
"""
function interpolate_datafields_2D(V::CartData, X, Y)
    X_vec = V.x.val[:, 1, 1]
    Y_vec = V.y.val[1, :, 1]
    return InterpolateDataFields2D_vecs(X_vec, Y_vec, V.z, V.fields, X, Y)
end

"""
    interpolate_datafields_2D(Original::CartData, New::CartData; Rotate=0.0, Translate=(0,0,0), Scale=(1.0,1.0,1.0))

Interpolates a data field `Original` on a 2D CartData grid `New`.
Typically used for horizontal surfaces. Points of `New` outside `Original` take the value at the nearest edge of `Original`.

Note: `Original` should have orthogonal coordinates. If it has not, e.g., because it was rotated, you'll have to specify the angle `Rotate` that it was rotated by

"""
function interpolate_datafields_2D(Original::CartData, New::CartData; Rotate = 0.0, Translate = (0.0, 0.0, 0.0), Scale = (1.0, 1.0, 1.0))
    if (Rotate != 0.0) || any(Translate .!= (0, 0, 0)) || any(Scale .!= (1.0, 1.0, 1.0))
        Original_r = rotate_translate_scale(Original, Rotate = -1.0 * Rotate, Translate = -1.0 .* Translate, Scale = Scale)
        New_r = rotate_translate_scale(New, Rotate = -1.0 * Rotate, Translate = -1.0 .* Translate, Scale = Scale)
    else
        Original_r = Original
        New_r = New
    end

    X_vec = Original_r.x.val[:, 1, 1]
    Y_vec = Original_r.y.val[1, :, 1]

    Xnew = New_r.x.val
    Ynew = New_r.y.val
    Znew, fields_new = GeophysicalModelGenerator.InterpolateDataFields2D_vecs(X_vec, Y_vec, Original_r.z, Original_r.fields, Xnew, Ynew)

    return CartData(New.x.val, New.y.val, Znew, fields_new, New.atts)
end

"""
    interpolate_datafields_2D(Original::GeoData, New::GeoData; Rotate=0.0, Translate=(0,0,0), Scale=(1.0,1.0,1.0))

Interpolates a data field `Original` on a 2D GeoData grid `New`.
Typically used for horizontal surfaces. Points of `New` outside `Original` take the value at the nearest edge of `Original`.

Note: `Original` should have orthogonal coordinates. If it has not, e.g., because it was rotated, you'll have to specify the angle `Rotate` that it was rotated by

"""
function interpolate_datafields_2D(Original::GeoData, New::GeoData; Rotate = 0.0, Translate = (0.0, 0.0, 0.0), Scale = (1.0, 1.0, 1.0))
    if (Rotate != 0.0) || any(Translate .!= (0, 0, 0)) || any(Scale .!= (1.0, 1.0, 1.0))
        Original_r = rotate_translate_scale(Original, Rotate = -1.0 * Rotate, Translate = -1.0 .* Translate, Scale = Scale)
        New_r = rotate_translate_scale(New, Rotate = -1.0 * Rotate, Translate = -1.0 .* Translate, Scale = Scale)
    else
        Original_r = Original
        New_r = New
    end

    Lon_vec = Original_r.lon.val[:, 1, 1]
    Lat_vec = Original_r.lat.val[1, :, 1]

    Lon_new = New_r.lon.val
    Lat_new = New_r.lat.val
    Znew, fields_new = GeophysicalModelGenerator.InterpolateDataFields2D_vecs(Lon_vec, Lat_vec, Original_r.depth, Original_r.fields, Lon_new, Lat_new)

    return GeoData(New.lon.val, New.lat.val, Znew, fields_new, New.atts)
end

"""
    Surf_interp = interpolate_datafields_2D(V::GeoData, x::AbstractRange, y::AbstractRange;  Lat::Number, Lon::Number)

Interpolates a 3D data set `V` with a projection point `proj=(Lat, Lon)` on a plane defined by `x` and `y`, where `x` and `y` are uniformly spaced.
Returns the 2D array `Surf_interp`.
"""
function interpolate_datafields_2D(V::GeoData, x::AbstractRange, y::AbstractRange; Lat = 49.9929, Lon = 8.2473)
    # Default: Lat=49.9929, Lon=8.2473 => Mainz (center of universe)
    proj = ProjectionPoint(; Lat = Lat, Lon = Lon)
    return interpolate_datafields_2D(V::GeoData, proj, x, y)
end

function interpolate_datafields_2D(LonLat::GeoData, proj::ProjectionPoint, x::AbstractRange, y::AbstractRange)
    cart_grid = CartData(xyz_grid(x, y, 0))
    tproj = project_CartData(cart_grid, LonLat, proj)
    return tproj.z.val[:, :, 1]
end

"""
    InterpolateDataFields2D_vecs(x_vec, y_vec, depth, fields_new, X, Y)

Interpolates a data field `V` on a 2D grid defined by `UTM`. Typically used for horizontal surfaces
"""
function InterpolateDataFields2D_vecs(EW_vec, NS_vec, depth, fields_new, EW, NS)
    EW_vec, NS_vec = as_range(EW_vec), as_range(NS_vec)     # evenly spaced axes: computed rather than searched lookups

    field_names = keys(fields_new)
    for i in 1:length(fields_new)
        data_tuple = fields_new[i]

        if (typeof(data_tuple) <: Tuple) & (!contains(String(field_names[i]), "colors"))
            # vector or anything that contains more than 1 field; each component takes the shape of EW, like scalar fields
            data_new = map(data_tuple) do d
                data_j = ndims(d) == 3 ? d[:, :, 1] : d
                interpol = linear_interpolation((EW_vec, NS_vec), ustrip.(data_j), extrapolation_bc = Flat())       # create interpolation object
                interpol.(EW, NS)
            end

        elseif contains(String(field_names[i]), "colors")
            # This is a 3D matrix. We need to interpolate each color channel separately, while using nearest neighbour
            # as you don't want to average colors

            # use nearest neighbour to interpolate data
            X, Y, _ = xyz_grid(EW_vec, NS_vec, depth.val[1])

            coord = [vec(X)'; vec(Y)']
            kdtree = KDTree(coord; leafsize = 10)
            points = [vec(EW)';vec(NS)']
            idx, dist = nn(kdtree, points)
            I_idx = CartesianIndices(axes(X))[idx]
            data_new = map(c -> reshape(c[I_idx], size(EW)), data_tuple)     # each channel takes the shape of EW, like the other fields

        else
            # scalar field
            if length(size(data_tuple)) == 3
                interpol = linear_interpolation((EW_vec, NS_vec), data_tuple[:, :, 1], extrapolation_bc = Flat())             # create interpolation object
            else
                interpol = linear_interpolation((EW_vec, NS_vec), data_tuple, extrapolation_bc = Flat())             # create interpolation object
            end

            data_new = interpol.(EW, NS)                                                  # interpolate data field
        end

        # replace the one
        new_field = NamedTuple{(field_names[i],)}((data_new,))                          # Create a tuple with same name
        fields_new = merge(fields_new, new_field)                                        # replace the field in fields_new

    end

    # Interpolate z-coordinate as well
    if length(size(depth)) == 3
        interpol = linear_interpolation((EW_vec, NS_vec), depth.val[:, :, 1], extrapolation_bc = Flat())             # create interpolation object
    else
        interpol = linear_interpolation((EW_vec, NS_vec), depth.val, extrapolation_bc = Flat())             # create interpolation object
    end
    depth_new = interpol.(EW, NS)

    return depth_new, fields_new
end

# Extracts a sub-data set using indices
function ExtractDataSets(V::AbstractGeneralGrid, iLon, iLat, iDepth)

    X, Y, Z = coordinate_grids(V)


    Lon = zeros(typeof(X[1]), length(iLon), length(iLat), length(iDepth))
    Lat = zeros(typeof(Y[1]), length(iLon), length(iLat), length(iDepth))
    Depth = zeros(typeof(Z[1]), length(iLon), length(iLat), length(iDepth))

    iLo = 1:length(iLon)
    iLa = 1:length(iLat)
    iDe = 1:length(iDepth)
    Lon[iLo, iLa, iDe] = X[iLon, iLat, iDepth]
    Lat[iLo, iLa, iDe] = Y[iLon, iLat, iDepth]
    Depth[iLo, iLa, iDe] = Z[iLon, iLat, iDepth]

    fields_new = V.fields
    field_names = keys(fields_new)
    for i in 1:length(V.fields)
        if typeof(V.fields[i]) <: Tuple
            # vector or anything that contains more than 1 field
            data_tuple = fields_new[i]      # we have a tuple (likely a vector field), so we have to loop
            data_array = zeros(typeof(data_tuple[1][1]), length(iLon), length(iLat), length(iDepth), length(data_tuple))      # create a 3D array that holds the 2D interpolated values

            for j in 1:length(data_tuple)
                data_field = data_tuple[j]
                data_array[:, :, :, j] = data_field[iLon, iLat, iDepth]
            end
            data_new = tuple([data_array[:, :, :, c] for c in 1:size(data_array, 4)]...)       # transform 4D matrix to tuple

        else
            # scalar field
            data_new = zeros(typeof(V.fields[i][1]), length(iLon), length(iLat), length(iDepth))
            data_new[iLo, iLa, iDe] = V.fields[i][iLon, iLat, iDepth]                                 # interpolate data field
        end

        # replace the one
        new_field = NamedTuple{(field_names[i],)}((data_new,))                          # Create a tuple with same name
        fields_new = merge(fields_new, new_field)                                        # replace the field in fields_new

    end


    # Create a GeoData struct with the newly interpolated fields
    return if isa(V, GeoData)
        Data_profile = GeoData(Lon, Lat, Depth, fields_new, V.atts)
    elseif isa(V, CartData)
        Data_profile = CartData(Lon, Lat, Depth, fields_new, V.atts)
    else
        error("Not yet implemented")
    end

end

"""
    V_sub = subtract_horizontalmean(V::AbstractArray{T, 3}; Percentage=false)

Subtracts the horizontal average of the 3D data array V.

If `Percentage=true`, the result is given as percentage; otherwise absolute values are returned

"""
function subtract_horizontalmean(V::AbstractArray{T, 3}; Percentage = false) where {T}
    V_sub = similar(V, Percentage ? Float64 : float(T))     # a percentage has no units

    for iLayer in axes(V, 3)
        layer = view(V, :, :, iLayer)
        average = mean(filter(!isnan, vec(layer)))

        if Percentage
            V_sub[:, :, iLayer] .= (ustrip.(layer) .- ustrip(average)) ./ ustrip(average) .* 100.0      # the result is normalized
        else
            V_sub[:, :, iLayer] .= layer .- average
        end
    end

    return V_sub
end

"""
    V_sub = subtract_horizontalmean(V::AbstractArray{T, 2}; Percentage=false)

Subtracts the horizontal average of the 2D data array V.

If `Percentage=true`, the result is given as percentage; otherwise absolute values are returned

"""
function subtract_horizontalmean(V::AbstractArray{T, 2}; Percentage = false) where {T}
    V_sub = similar(V, Percentage ? Float64 : float(T))     # a percentage has no units

    for iLayer in axes(V, 2)
        layer = view(V, :, iLayer)
        average = mean(filter(!isnan, layer))

        if Percentage
            V_sub[:, iLayer] .= (ustrip.(layer) .- ustrip(average)) ./ ustrip(average) .* 100.0      # the result is normalized
        else
            V_sub[:, iLayer] .= layer .- average
        end
    end

    return V_sub
end

"""
    parse_columns_CSV(data_file, num_columns)

This parses numbers from CSV file that is read in with `CSV.File`.
That is useful in case the CSV files has tables that contain both strings (e.g., station names) and numbers (lat/lon/height) and you are only interested in the numbers


# Example
This example assumes that the data starts at line 18, that the columns are separated by spaces, and that it contains at most 4 columns with data:
```julia-repl
julia> using CSV
julia> data_file        =   CSV.File("FileName.txt",datarow=18,header=false,delim=' ')
julia> data = parse_columns_CSV(data_file, 4)
```

"""
function parse_columns_CSV(data_file, num_columns)
    data = zeros(size(data_file, 1), num_columns)
    for (row_num, row) in enumerate(data_file)
        num = 0
        for i in 1:length(row)
            ismissing(row[i]) && continue # empty cells from repeated delimiters
            v = row[i] isa Real ? Float64(row[i]) : tryparse(Float64, row[i])
            if !isnothing(v)
                num += 1
                data[row_num, num] = v
            end

        end
    end
    return data
end

"""
    votemap(DataSets::Vector{GeoData}, criteria::Vector{String}, dims=(50,50,50))

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
You can create a votemap which combines the two data sets with:
```julia
julia> Data_VoteMap = votemap([Data_Zhao_Pwave,DataKoulakov_Alps],["dVp_Percentage>2.5","dVp_percentage>3.0"])
GeoData
  size  : (50, 50, 50)
  lon   ϵ [ 4.0 : 18.0]
  lat   ϵ [ 38.0 : 49.01197604790419]
  depth ϵ [ -700.0 km : -10.0 km]
  fields: (:votemap,)
```

You can also create a votemap of a single dataset:
```julia
julia> Data_VoteMap = votemap(Data_Zhao_Pwave,"dVp_Percentage>2.5", dims=(50,51,52))
GeoData
  size  : (50, 51, 52)
  lon   ϵ [ 0.0 : 18.0]
  lat   ϵ [ 38.0 : 51.95]
  depth ϵ [ -1001.0 km : -1.0 km]
  fields: (:votemap,)
```

"""
function votemap(DataSets::Vector{GeoData}, criteria::Vector{String}; dims = (50, 50, 50))

    numDataSets = length(DataSets)

    if length(criteria) != numDataSets
        error("Need the same number of criteria as the number of data sets")
    end

    # Determine the overlapping lon/lat/depth regions of all datasets
    lon_limits = [minimum(DataSets[1].lon.val);        maximum(DataSets[1].lon.val)]
    lat_limits = [minimum(DataSets[1].lat.val);        maximum(DataSets[1].lat.val)]
    z_limits = [minimum(DataSets[1].depth.val);      maximum(DataSets[1].depth.val)]
    for i in 1:numDataSets
        lon_limits[1] = maximum([lon_limits[1]  minimum(DataSets[i].lon.val)])
        lon_limits[2] = minimum([lon_limits[2]  maximum(DataSets[i].lon.val)])

        lat_limits[1] = maximum([lat_limits[1]  minimum(DataSets[i].lat.val)])
        lat_limits[2] = minimum([lat_limits[2]  maximum(DataSets[i].lat.val)])

        z_limits[1] = maximum([z_limits[1]    minimum(DataSets[i].depth.val)])
        z_limits[2] = minimum([z_limits[2]    maximum(DataSets[i].depth.val)])
    end

    for (name, limits) in (("lon", lon_limits), ("lat", lat_limits), ("depth", z_limits))
        limits[1] <= limits[2] || error("votemap: the data sets do not overlap in $name.")
    end

    # Loop over all datasets, and interpolate the data set to the new (usually smaller) domain
    votemap = zeros(Int64, dims)
    lon = lat = depth = nothing   # coordinates of the common grid, taken from the first interpolated data set
    for i in 1:numDataSets
        VoteMap_Local = zeros(Int64, dims)

        # Interpolate data set to smaller domain
        DataSet = extract_subvolume(DataSets[i]; Interpolate = true, Lon_level = lon_limits, Lat_level = lat_limits, Depth_level = z_limits, dims = dims)
        if i == 1
            lon, lat, depth = DataSet.lon.val, DataSet.lat.val, DataSet.depth.val
        end

        # Extract the criteria to evaluate
        expr = Meta.parse(criteria[i])      # the expression, such as Vs>1.0

        # Extract data field
        if !haskey(DataSet.fields, expr.args[2])
            error("The GeoData set does not have the field: $(expr.args[2])")
        end

        Array3D = ustrip.(DataSet.fields[expr.args[2]])                   # strip units, just in case

        # Modify the value, to be Array3D
        expr_mod = Expr(:call, expr.args[1], :($Array3D), expr.args[3])       # modify the original expression to use Array3D as variable name

        # The expression should have a ".", such as Array .> 1.0. If not, it will not apply this in a pointwise manner
        #   Here, we add this dot if it is not there yet
        if !startswith(String(expr_mod.args[1]), ".")
            expr_mod.args[1] = Symbol(".", expr_mod.args[1])
        end

        ind = eval(expr_mod)     # evaluate the modified expression
        VoteMap_Local[ind] .= 1                  # assign vote-map

        votemap = votemap + VoteMap_Local        # Sum
    end

    # Construct GeoData set that holds the votemap (makes it easier to write paraview files)
    VoteData = GeoData(lon, lat, depth, (votemap = votemap,))

    return VoteData
end

# Make this work for single data sets as well
function votemap(DataSets::GeoData, criteria::String; dims = (50, 50, 50))
    return votemap([DataSets], [criteria]; dims = dims)
end

"""
    Data_R = rotate_translate_scale(Data::Union{ParaviewData, CartData}; Rotate=0, Translate=(0,0,0), Scale=(1.0,1.0,1.0), Xc=(0.0,0.0))

Returns a copy of the Cartesian dataset `Data` that is scaled, rotated and translated; `Data` is not modified.

# Parameters
Note that we apply the transformations in exactly this order:
-   `Scale`:        scaling applied to the `x,y,z` coordinates of the data set
-   `Rotate`:       rotation around the `x/y` axis, around the point `Xc`
-   `Translate`:    translation
-   `Xc`:           center of rotation (default: the origin)

# Example
```julia
julia> X,Y,Z   =   xyz_grid(10:20,30:40,-50:-10);
julia> Data_C  =   CartData(X,Y,Z,(Depth=Z,));
julia> Data_R = rotate_translate_scale(Data_C, Rotate=30)
CartData
    size    : (11, 11, 41)
    x       ϵ [ -11.339745962155614 : 2.320508075688771]
    y       ϵ [ 30.980762113533157 : 44.64101615137754]
    z       ϵ [ -50.0 : -10.0]
    fields  : (:Depth,)
```
"""
function rotate_translate_scale(Data::Union{ParaviewData, CartData}; Rotate::Number = 0.0, Translate = (0, 0, 0), Scale = (1.0, 1.0, 1.0), Xc = (0.0, 0.0))

    X, Y, Z = copy(Data.x.val), copy(Data.y.val), copy(Data.z.val)      # Extract coordinates
    Xr, Yr, Zr = X, Y, Z                                                     # Rotated coordinates

    # 1) Scaling
    if length(Scale) == 1
        Scale = [Scale, Scale, Scale]
    end
    Xr .*= Scale[1]
    Yr .*= Scale[2]
    Zr .*= Scale[3]


    # 2) 2D rotation around X/Y axis, around Xc
    R = [cosd(Rotate) -sind(Rotate); sind(Rotate) cosd(Rotate)]  # 2D rotation matrix

    for i in eachindex(X)
        Rot_XY = R * [X[i] - Xc[1]; Y[i] - Xc[2]]
        Xr[i] = Rot_XY[1] + Xc[1]
        Yr[i] = Rot_XY[2] + Xc[2]
    end

    # 3) Add translation
    Xr .+= Translate[1]
    Yr .+= Translate[2]
    Zr .+= Translate[3]

    # Modify original structure
    if isa(Data, ParaviewData)
        Data_r = ParaviewData(Xr, Yr, Zr, Data.fields)
    else
        Data_r = CartData(Xr, Yr, Zr, Data.fields, Data.atts)
    end

    return Data_r
end


"""
    lithostatic_pressure!(Plithos::Array, Density::Array, dz::Number; g=9.81)

Computes lithostatic pressure from a 3D density array, assuming constant soacing `dz` in vertical direction. Optionally, the gravitational acceleration `g` can be specified.

"""
function lithostatic_pressure!(Plithos::Array{T, N}, Density::Array{T, N}, dz::Number; g = 9.81) where {T, N}

    Plithos[:] = Density * dz * g

    selectdim(Plithos, N, size(Plithos)[N]) .= 0      # set upper row to zero

    Plithos[:] = reverse!(cumsum(reverse!(Plithos), dims = N))

    return nothing
end

"""
    inpolygon!(INSIDE::Matrix, PolyX::Vector, PolyY::Vector, X::Matrix, Y::Matrix; fast=false)

Checks if points given by matrices `X` and `Y` are in or on (both cases return true) a polygon given by `PolyX` and `PolyY`. Boolean `fast` will trigger faster version that may miss points that are exactly on the edge of the polygon. Speedup is a factor of 3.

"""
function inpolygon!(INSIDE::Matrix{Bool}, PolyX::Vector{T}, PolyY::Vector{T}, X::Matrix{T}, Y::Matrix{T}; fast = false) where {T <: Real}
    return if fast
        for j in 1:size(X, 2)
            for i in 1:size(X, 1)
                INSIDE[i, j] = inpoly_fast(PolyX, PolyY, X[i, j], Y[i, j])
            end
        end
    else
        for j in 1:size(X, 2)
            for i in 1:size(X, 1)
                INSIDE[i, j] = inpoly_either(PolyX, PolyY, X[i, j], Y[i, j])
            end
        end
    end
end

"""
    inpolygon!(inside::Vector, PolyX::Vector, PolyY::Vector, x::Vector, y::Vector; fast=false)

Same as above but `inside`, `X` and `Y` and are vectors.

"""
function inpolygon!(inside::Vector{Bool}, PolyX::AbstractVector{T}, PolyY::AbstractVector{T}, x::Vector{T}, y::Vector{T}; fast = false) where {T <: Real}
    return if fast
        for i in eachindex(x)
            inside[i] = inpoly_fast(PolyX, PolyY, x[i], y[i])
        end
    else
        for i in eachindex(x)
            inside[i] = inpoly_either(PolyX, PolyY, x[i], y[i])
        end
    end
end

# In or on the polygon for either order of the axes: the test of `inpolygon!` without `fast`
inpoly_either(PolyX, PolyY, x, y) = inpoly(PolyX, PolyY, x, y) || inpoly(PolyY, PolyX, y, x)

"""
    inpoly(PolyX::AbstractVector, PolyY::AbstractVector, x, y)

Checks if a point given by x and y is in or on (both cases return true) a polygon given by PolyX and PolyY. This function should be used through inpolygon!().

"""
function inpoly(PolyX::AbstractVector, PolyY::AbstractVector, x, y)
    Base.require_one_based_indexing(PolyX, PolyY)
    inside1, inside2, inside3, inside4 = false, false, false, false
    n = length(PolyX)
    for i in eachindex(PolyX, PolyY)
        j = i - 1
        j += n * (j < 1)
        xi, yi, xj, yj = PolyX[i], PolyY[i], PolyX[j], PolyY[j]

        con1 = ((yi > y) != (yj > y))
        con2 = ((yi >= y) != (yj >= y))
        (con1 || con2) || continue
        xc = (xj - xi) * (y - yi) / (yj - yi + eps()) + xi       # x where the edge crosses the height y
        inside1 ⊻= con1 & (x > xc)
        inside2 ⊻= con1 & (x >= xc)
        inside3 ⊻= con2 & (x > xc)
        inside4 ⊻= con2 & (x >= xc)
    end
    return inside1 || inside2 || inside3 || inside4
end

"""
    inpoly_fast(PolyX::Vector, PolyY::Vector, x::Number, y::Number, iSteps::Vector, jSteps::)

Faster version of inpoly() but will miss some points that are on the edge of the polygon.

"""
function inpoly_fast(PolyX::Vector{T}, PolyY::Vector{T}, x::T, y::T) where {T <: Real}
    inside = false
    n = length(PolyX)
    for i in eachindex(PolyX)
        j = i - 1
        j += n * (j < 1)
        xi = PolyX[i]
        yi = PolyY[i]
        xj = PolyX[j]
        yj = PolyY[j]

        if ((yi > y) != (yj > y)) && (x > (xj - xi) * (y - yi) / (yj - yi + eps()) + xi)
            inside = !inside
        end
    end
    return inside
end
