# This is data_import.jl
#
# This file contains functions to import different data types. 
#
# Author: Marcel Thielmann, 05/2021


# import CSV data using standard library functions
# here we assume that the data is indeed comma separated and that comments are preceded with a "#"





function ReadCSV_LatLon(filename)
    # import data from file with coordinates given in lat/lon/depth format and additional data given in additional columns
    data,hdr = readdlm(filename,',', Float64,'\n'; header=true, skipblanks=true, comments=true, comment_char='#')
    
    # initialize data vectors -> to be redone with ScatteredPoints structure
    ndata   = size(data,1) # number of entries
    nfields = size(data,2) # number of fields
    lon     = zeros(Float64, n, 1)
    lat     = zeros(Float64, n, 1)
    depth   = zeros(Float64, n, 1)

    # ISSUE: we may have additional data that we want to extract from the csv, but how to do that?
    # We could define a structure with three fields: variable name, unit and values
    # using this structure, we could treat all possible variables

    # assign data to output tuple
    return [lon lat depth]

end