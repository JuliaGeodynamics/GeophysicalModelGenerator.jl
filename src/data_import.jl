# This is data_import.jl
#
# This file contains functions to import different data types. 
#
# Author: Marcel Thielmann, 05/2021


# import CSV data using standard library functions
# here we assume that the data is indeed comma separated and that comments are preceded with a "#"
# this is basically just using a function from the CSV package

function ReadCSV_LatLon(filename)
    # import data from file with coordinates given in lat/lon/depth format and additional data given in additional columns
    data,hdr = readdlm(filename,',', Float64,'\n'; header=true, skipblanks=true, comments=true, comment_char='#')
    
    # initialize data
    ndata   = size(data,1) # number of entries
    nfields = size(data,2) # number of fields
    lon     = zeros(Float64, n, 1)
    lat     = zeros(Float64, n, 1)
    depth   = zeros(Float64, n, 1)

    # assign data to vectors


end