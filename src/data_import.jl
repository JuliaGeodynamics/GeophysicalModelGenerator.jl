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

    # get the fields for lon/lat/depth
    for ifield = ifield = 1:nfields
        if occursin("lon",hdr[ifield])
            lon_ind = ifield;
            varname = GetVariableName(hdr[ifield])# get variable name
            varunit = GetVariableUnit(hdr[ifield])# get variable unit
            LonData = ValueList(varname,varunit,data[1:end,ifield])
        elseif occursin("lat",hdr[ifield])
            lat_ind = ifield;
            varname = GetVariableName(hdr[ifield])# get variable name
            varunit = GetVariableUnit(hdr[ifield])# get variable unit
            LatData = ValueList(varname,varunit,data[1:end,ifield])
        elseif occursin("depth",hdr[ifield])
            # ISSUE: WE DEFINE DEPTH AS NEGATIVE, BUT HOW DO WE SET THAT?
            # WE COULD ADD A FLAG THAT INDICATES THE DEPTH CONVENTION AND 
            # TREAT IT ACCORDINGLY
            depth_ind = ifield;
            varname = GetVariableName(hdr[ifield])# get variable name
            varunit = GetVariableUnit(hdr[ifield])# get variable unit

            if cmp(DepthCon,"positive") # if depth is given as positive values, convert to negative
                DepthData = ValueList(varname,varunit,-1*data[1:end,ifield])
            elseif cmp(DepthCon,"negative")
                DepthData = ValueList(varname,varunit,data[1:end,ifield])
            else # default behaviour assumes that dpeth is negative
                DepthData = ValueList(varname,varunit,data[1:end,ifield])
            end
            DepthData = ValueList(varname,varunit,data[1:end,ifield])
        end
    end

    # assign rest of the data to the values tuple
    tmp = (varnames = hdr[3:end],vals = data[1:end,3:end])

    # initialize data structure
    importdata = GeoData(LonData,LatData,DepthData,tmp)

    # assign data to output
    return importdata 

end

function GetVariableName(inputstring::String)
    # assume that if the string contains a unit, it is given in brackets (),[],{}
    indfirst = nothing
    iloop     = 1;
    str2find = ["(","[","{"]
    # find first occurence of one of the brackets
    while isnothing(ind_first)
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
        return inputstring(1:indfirst-1)
    end
end

function GetVariableUnit(inputstring::String)
    # assume that if the string contains a unit, it is given in brackets (),[],{}
    indfirst = nothing
    iloop     = 1;
    firststr2find = ["(","[","{"]
    laststr2find = [")","]","}"]
    while isnothing(ind_first)
        indfirst = findfirst(firststr2find[iloop],inputstring)
        iloop = iloop + 1
        if iloop>length(str2find)
            break
        end
    end

    if isnothing(indfirst)
        return "none"
    else
        indlast = findfirst(laststr2find[iloop-1],inputstring)
        return inputstring(indfirst+1:indlast-1)
    end

end