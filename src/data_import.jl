# This is data_import.jl
#
# This file contains functions to import different data types.
#
# Author: Marcel Thielmann, 05/2021

using LightXML

export Screenshot_To_GeoData, Screenshot_To_CartData, Screenshot_To_UTMData, GetLonLatDepthMag_QuakeML

# import CSV data using standard library functions
# here we assume that the data is indeed comma separated and that comments are preceded with a "#"

function ReadCSV_LatLon(filename::AbstractString,DepthCon::AbstractString)
    # import data from file with coordinates given in lat/lon/depth format and additional data given in additional columns
    # the idea here is to assign the data to a structure of the type GeoData which will then be used for further processing
    data,hdr = readdlm(filename,',', Float64,'\n'; header=true, skipblanks=true, comments=true, comment_char='#')

    # initialize array of structures to store the data
    # while doing so, separate the unit from the variable name
    ndata   = size(data,1) # number of entries
    nfields = size(data,2) # number of fields

    # declare some variables as local, otherwise they are not known outside of the following for loop
    local LonData
    local LatData
    local DepthData
    local vals_range
    vals_range = zeros(Int64,nfields-3)
    ivals = 1;
    # get the fields for lon/lat/depth
    for ifield = 1:nfields
        if occursin("lon",hdr[ifield])
            lon_ind = ifield;
            varname = GetVariableName(hdr[ifield])# get variable name
            varunit = GetVariableUnit(hdr[ifield])# get variable unit
            LonData = data[1:end,ifield]
        elseif occursin("lat",hdr[ifield])
            lat_ind = ifield;
            varname = GetVariableName(hdr[ifield])# get variable name
            varunit = GetVariableUnit(hdr[ifield])# get variable unit
            LatData = data[1:end,ifield]
        elseif occursin("depth",hdr[ifield])
            # ISSUE: WE DEFINE DEPTH AS NEGATIVE, BUT HOW DO WE SET THAT?
            # WE COULD ADD A FLAG THAT INDICATES THE DEPTH CONVENTION AND
            # TREAT IT ACCORDINGLY
            depth_ind = ifield;
            varname = GetVariableName(hdr[ifield])# get variable name
            varunit = GetVariableUnit(hdr[ifield])# get variable unit

            # take care of positive or negative convection for depth (here we use negative)
            if cmp(DepthCon,"positive")==0 # if depth is given as positive values, convert to negative
                DepthData = -1*data[1:end,ifield]
            elseif cmp(DepthCon,"negative")==0
                DepthData = data[1:end,ifield]
            else # default behaviour assumes that dpeth is negative
                DepthData = data[1:end,ifield]
            end

            # if depth is given in m, convert to km
            if cmp(varunit,"m")==0
                DepthData = DepthData./1e3;
            end

        else
            vals_range[ivals] = ifield
            ivals = ivals+1
        end
    end


    # create named tuple for additional data
    tmp_hdr  = hdr[vals_range];
    tmp_data = data[1:end,vals_range];

    nhdr = size(tmp_hdr,1)
    tmp_vec = Vector{Vector{Float64}}(undef, nhdr) # this is used for later tuple creation, I haven't found a better way around

    for ihdr = 1:nhdr

        # take care of the header strings
        varname = GetVariableName(tmp_hdr[ihdr])# get variable name
        varunit = GetVariableUnit(tmp_hdr[ihdr])# get variable unit
        if cmp(varunit,"%")==0
            tmp_hdr[ihdr] = string(varname,"_percentage")
        else
            tmp_hdr[ihdr] = string(varname,"_",varunit)
        end

        # take care of the matrix columns
        tmp_vec[ihdr] = tmp_data[1:end,ihdr];
    end

    hdr_tpl  = Tuple(Symbol(x) for x in tmp_hdr) # convert header to tuple
    data_tpl = Tuple.(tmp_vec for i in size(tmp_vec,1)) # convert data to tuple
    tmp = NamedTuple{hdr_tpl}(data_tpl)

    println(typeof(tmp))

    # initialize data structure
    importdata = GeoData(LonData,LatData,DepthData,tmp)

    # assign data to output
    return importdata

end

function GetVariableName(inputstring::SubString{String})
    # convert input to normal String
    inputstring = String(inputstring)
    # assume that if the string contains a unit, it is given in brackets (),[],{}
    indfirst = nothing
    iloop     = 1
    str2find = ["(","[","{"]
    # find first occurrence of one of the brackets
    while isnothing(indfirst)
        indfirst = findfirst(str2find[iloop],inputstring)
        iloop = iloop + 1
        if iloop>length(str2find)
            break
        end
    end

    # either return the whole inputstring or only the part before the unit
    if isnothing(indfirst)
        return inputstring
    else
        indfirst = indfirst[1]-1
        return inputstring[1:indfirst]
    end
end

function GetVariableUnit(inputstring::SubString{String})
    # convert input to normal String
    inputstring = String(inputstring)
    # assume that if the string contains a unit, it is given in brackets (),[],{}
    indfirst = nothing
    iloop     = 1;
    firststr2find = ["(","[","{"]
    laststr2find = [")","]","}"]
    while isnothing(indfirst)
        indfirst = findfirst(firststr2find[iloop],inputstring)
        iloop = iloop + 1
        if iloop>length(firststr2find)
            break
        end
    end

    if isnothing(indfirst)
        return "none"
    else
        indlast = findfirst(laststr2find[iloop-1],inputstring)
        indfirst = indfirst[1]+1
        indlast = indlast[1]-1
        return inputstring[indfirst:indlast]
    end

end



"""
    Screenshot_To_GeoData(filename::String, Corner_LowerLeft, Corner_UpperRight; Corner_LowerRight=nothing, Corner_UpperLeft=nothing, Cartesian=false, UTM=false, UTMzone, isnorth=true, fieldname::Symbol=:colors)

Take a screenshot of Georeferenced image either a `lat/lon`, `x,y` (if `Cartesian=true`) or in UTM coordinates (if `UTM=true`) at a given depth or along profile and converts it to a `GeoData`, `CartData` or `UTMData` struct, which can be saved to Paraview

The lower left and upper right coordinates of the image need to be specified in tuples of `(lon,lat,depth)` or `(UTM_ew, UTM_ns, depth)`, where `depth` is negative inside the Earth (and in km).

The lower right and upper left corners can be specified optionally (to take non-orthogonal images into account). If they are not specified, the image is considered orthogonal and the corners are computed from the other two.

*Note*: if your data is in `UTM` coordinates you also need to provide the `UTMzone` and whether we are on the northern hemisphere or not (`isnorth`).
"""
function Screenshot_To_GeoData(filename::String, Corner_LowerLeft, Corner_UpperRight; Corner_LowerRight=nothing, Corner_UpperLeft=nothing, Cartesian=false, UTM=false, UTMzone=nothing, isnorth::Bool=true, fieldname::Symbol=:colors)

    img         =   load(filename)      # load image

    # Define lon/lat/depth of lower left corner

    # try to determine if this is a horizontal profile or not
    if abs(Corner_UpperRight[3]-Corner_LowerLeft[3])>0.0
        DepthProfile = true
    else
        DepthProfile = false
    end

    # We should be able to either define 4 corners or only 2 and reconstruct the other two from the
    if isnothing(Corner_LowerRight) || isnothing(Corner_UpperLeft)
        if DepthProfile
            Corner_LowerRight  = (Corner_UpperRight[1], Corner_UpperRight[2], Corner_LowerLeft[3])
            Corner_UpperLeft   = (Corner_LowerLeft[1],  Corner_LowerLeft[2], Corner_UpperRight[3])
        else
            Corner_LowerRight  = (Corner_UpperRight[1], Corner_LowerLeft[2],  Corner_LowerLeft[3])
            Corner_UpperLeft   = (Corner_LowerLeft[1],  Corner_UpperRight[2], Corner_UpperRight[3])
        end
    end

    # Print overview of the 4 corners here:
    if Cartesian
        println("Extracting CartData from: $(filename)")
        println("           └ Corners:         x         y         z")
        println("              └ lower left  = ($(rpad( Corner_LowerLeft[1],7)), $(rpad( Corner_LowerLeft[2],7)),  $(rpad( Corner_LowerLeft[3],7)))")
        println("              └ lower right = ($(rpad(Corner_LowerRight[1],7)), $(rpad(Corner_LowerRight[2],7)),  $(rpad(Corner_LowerRight[3],7)))")
        println("              └ upper left  = ($(rpad( Corner_UpperLeft[1],7)), $(rpad( Corner_UpperLeft[2],7)),  $(rpad( Corner_UpperLeft[3],7)))")
        println("              └ upper right = ($(rpad(Corner_UpperRight[1],7)), $(rpad(Corner_UpperRight[2],7)),  $(rpad(Corner_UpperRight[3],7)))")
    end
    if UTM
        if isnothing(UTMzone)
            error("You need to specify UTMzone and isnorth if reading in UTM data.")
        end
        println("Extracting UTMData from: $(filename)")
        if isnorth
        println("       UTM Zone $(UTMzone) Northern Hemisphere")
        else
        println("       UTM Zone $(UTMzone) Southern Hemisphere")
        end
        println("           └ Corners:         E-W (x)  | N-S (y) | depth (z)")
        println("              └ lower left  = ($(rpad( Corner_LowerLeft[1],7)), $(rpad( Corner_LowerLeft[2],7)),  $(rpad( Corner_LowerLeft[3],7)))")
        println("              └ lower right = ($(rpad(Corner_LowerRight[1],7)), $(rpad(Corner_LowerRight[2],7)),  $(rpad(Corner_LowerRight[3],7)))")
        println("              └ upper left  = ($(rpad( Corner_UpperLeft[1],7)), $(rpad( Corner_UpperLeft[2],7)),  $(rpad( Corner_UpperLeft[3],7)))")
        println("              └ upper right = ($(rpad(Corner_UpperRight[1],7)), $(rpad(Corner_UpperRight[2],7)),  $(rpad(Corner_UpperRight[3],7)))")
    end
    if (!Cartesian) && (!UTM)
        println("Extracting GeoData from: $(filename)")
        println("           └ Corners:         lon       lat       depth")
        println("              └ lower left  = ($(rpad( Corner_LowerLeft[1],7)), $(rpad( Corner_LowerLeft[2],7)),  $(rpad( Corner_LowerLeft[3],7)))")
        println("              └ lower right = ($(rpad(Corner_LowerRight[1],7)), $(rpad(Corner_LowerRight[2],7)),  $(rpad(Corner_LowerRight[3],7)))")
        println("              └ upper left  = ($(rpad( Corner_UpperLeft[1],7)), $(rpad( Corner_UpperLeft[2],7)),  $(rpad( Corner_UpperLeft[3],7)))")
        println("              └ upper right = ($(rpad(Corner_UpperRight[1],7)), $(rpad(Corner_UpperRight[2],7)),  $(rpad(Corner_UpperRight[3],7)))")
    end

    # Reconstruct the 4 corners into a matrix
    i = 1; Corners_lon     = [Corner_UpperLeft[i] Corner_UpperRight[i]; Corner_LowerLeft[i] Corner_LowerRight[i]; ]
    i = 2; Corners_lat     = [Corner_UpperLeft[i] Corner_UpperRight[i]; Corner_LowerLeft[i] Corner_LowerRight[i]; ]
    i = 3; Corners_depth   = [Corner_UpperLeft[i] Corner_UpperRight[i]; Corner_LowerLeft[i] Corner_LowerRight[i]; ]

   # i = 1; Corners_lon     = [Corner_LowerLeft[i] Corner_LowerRight[i]; Corner_UpperLeft[i] Corner_UpperRight[i];]
   # i = 2; Corners_lat     = [Corner_LowerLeft[i] Corner_LowerRight[i]; Corner_UpperLeft[i] Corner_UpperRight[i];]
   # i = 3; Corners_depth   = [Corner_LowerLeft[i] Corner_LowerRight[i]; Corner_UpperLeft[i] Corner_UpperRight[i];]


    # Extract the colors from the grid
    img_RGB     =   convert.(RGB, img)     # convert to  RGB data

    # extract the red-green-blue values from the image
    r           =   zeros(size(img_RGB))
    g           =   zeros(size(img_RGB))
    b           =   zeros(size(img_RGB))
    for i in eachindex(g)
        r[i] = Float64(img_RGB[i].r)
        g[i] = Float64(img_RGB[i].g)
        b[i] = Float64(img_RGB[i].b)
    end

    # Construct depth, lon and lat 2D grids from the corner points through linear interpolation
    grid_size               =   size(r)
    xs                      =   [1,grid_size[1]];
    zs                      =   [1,grid_size[2]];
    interp_linear_lon       =   linear_interpolation((xs, zs), Corners_lon)      # create interpolation object
    interp_linear_lat       =   linear_interpolation((xs, zs), Corners_lat)       # create interpolation object
    interp_linear_depth     =   linear_interpolation((xs, zs), Corners_depth)     # create interpolation object

    # Interpolate
    X_int,Y_int,Depth       =   XYZGrid(1:grid_size[1],1:grid_size[2],0)
    X                       =   interp_linear_lon.(X_int,   Y_int);
    Y                       =   interp_linear_lat.(X_int,   Y_int);
    Depth                   =   interp_linear_depth.(X_int, Y_int);

    # Transfer to 3D arrays (check if needed or not; if yes, redo error message in struct routine)
    red                     =   zeros(size(Depth)); red[:,:,1]   = r;
    green                   =   zeros(size(Depth)); green[:,:,1] = g;
    blue                    =   zeros(size(Depth)); blue[:,:,1]  = b;

    # Create GeoData structure - NOTE: RGB data must be 2D matrixes, not 3D!
    color_data = NamedTuple{(fieldname,)}(((red,green,blue),));

    if Cartesian
        data_Image              =   CartData(X, Y, Depth, color_data)
    end
    if UTM
        data_Image              =   UTMData(X, Y, Depth, UTMzone, isnorth, color_data)
    end
    if (!Cartesian) && (!UTM)
        data_Image              =   GeoData(X, Y, Depth, color_data)
    end
    return data_Image
end


"""
    Data = Screenshot_To_CartData(filename::String, Corner_LowerLeft, Corner_UpperRight; Corner_LowerRight=nothing, Corner_UpperLeft=nothing)

Does the same as `Screenshot_To_GeoData`, but returns a `CartData` structure
"""
function Screenshot_To_CartData(filename::String, Corner_LowerLeft, Corner_UpperRight; Corner_LowerRight=nothing, Corner_UpperLeft=nothing, fieldname::Symbol=:colors)


    # first create a GeoData struct
    Data_Cart = Screenshot_To_GeoData(filename, Corner_LowerLeft, Corner_UpperRight; Corner_LowerRight=Corner_LowerRight, Corner_UpperLeft=Corner_UpperLeft, Cartesian=true, fieldname=fieldname)

    return Data_Cart

end

"""
    Data = Screenshot_To_UTMData(filename::String, Corner_LowerLeft, Corner_UpperRight; Corner_LowerRight=nothing, Corner_UpperLeft=nothing, UTMzone::Int64=nothing, isnorth::Bool=true, fieldname=:colors)

Does the same as `Screenshot_To_GeoData`, but returns for UTM data
Note that you have to specify the `UTMzone` and `isnorth`
"""
function Screenshot_To_UTMData(filename::String, Corner_LowerLeft, Corner_UpperRight; Corner_LowerRight=nothing, Corner_UpperLeft=nothing, UTMzone::Int64=nothing, isnorth::Bool=true, fieldname::Symbol=:colors)

      # first create a GeoData struct
      Data_UTM = Screenshot_To_GeoData(filename, Corner_LowerLeft, Corner_UpperRight; Corner_LowerRight=Corner_LowerRight, Corner_UpperLeft=Corner_UpperLeft, Cartesian=false, UTM=true, UTMzone=UTMzone, isnorth=isnorth, fieldname=fieldname)
      return Data_UTM
end

"""
    Data = GetLonLatDepthMag_QuakeML(filename::String)

Extracts longitude, latitude, depth and magnitude from a QuakeML file that has been e.g. downloaded from ISC. The data is then returned in GeoData format.
"""
function GetLonLatDepthMag_QuakeML(filename::String)
    # The QuakeML format consists of a tree with quite a lot of branches, so we have to traverse it to quite some extent to get the desired values
    # using LightXML: extension???
    xdoc = parse_file(filename); # parse the whole file
    xroot =root(xdoc);
    catalogues = get_elements_by_tagname(xroot,"eventParameters");
    catalogue  = catalogues[1];
    events     = get_elements_by_tagname(catalogue,"event"); # now those are all events
    num_events = size(events,1);

    # allocate, lat,lon,depth,magnitude
    lon    = zeros(num_events,1);
    lat    = zeros(num_events,1);
    depth  = zeros(num_events,1);
    mag    = zeros(num_events,1);

    # now loop over the events and assign the respective values
    for ievent = 1:num_events
        tmp_event    = events[ievent];
        origin = get_elements_by_tagname(events[ievent], "origin");
        magnitude = get_elements_by_tagname(events[ievent], "magnitude");

        # this is a bit dirty, if you find a better/cleaner way, be my guest...
        lon[ievent]   = parse(Float64,string(collect(child_nodes(collect(child_elements(get_elements_by_tagname(origin[1], "longitude")[1]))[1]))[1]))
        lat[ievent]   = parse(Float64,string(collect(child_nodes(collect(child_elements(get_elements_by_tagname(origin[1], "latitude")[1]))[1]))[1]))
        depth[ievent] = parse(Float64,string(collect(child_nodes(collect(child_elements(get_elements_by_tagname(origin[1], "depth")[1]))[1]))[1]))
        mag[ievent]   = parse(Float64,string(collect(child_nodes(get_elements_by_tagname(get_elements_by_tagname(magnitude[1],"mag")[1],"value")[1]))[1]));
    end
    
    Data_ISC = GeoData(lon,lat,-1*depth/1e3,(Magnitude=mag,Depth=-1*depth/1e3*km));
    return Data_ISC
end
