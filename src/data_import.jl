# This is data_import.jl
#
# This file contains functions to import different data types. 
#
# Author: Marcel Thielmann, 05/2021


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
    # find first occurence of one of the brackets
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