# This is data_import.jl
#
# This file contains functions to import different data types. 
#
# Author: Marcel Thielmann, 05/2021


# import CSV data using standard library functions
# here we assume that the data is indeed comma separated and that comments are preceded with a "#"





function ReadCSV_LatLon(filename::String)
    # import data from file with coordinates given in lat/lon/depth format and additional data given in additional columns
    # the idea here is to assign the data to a structure of the type ScatteredPoints
    data,hdr = readdlm(filename,',', Float64,'\n'; header=true, skipblanks=true, comments=true, comment_char='#')
    
    # initialize array of structures to store the data
    # while doing so, separate the unit from the variable name 
    ndata   = size(data,1) # number of entries
    nfields = size(data,2) # number of fields

    importdata = Vector{ScatteredPoints}(undef, ndata)
     for ifield = 1:nfields
        varname = GetVariableName(hdr[ifield])# get variable name
        varunit = GetVariableUnit(hdr[ifield])# get variable unit
        # make sure that lon, lat and depth are used, not longitude,latitude,depth
        v[i] = ScatteredPoints(varname, varunit,data[1:end,1])
    end

    # ISSUE: we may have additional data that we want to extract from the csv, but how to do that?
    # We could define a structure with three fields: variable name, unit and values
    # using this structure, we could treat all possible variables

    # assign data to output tuple
    return v

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