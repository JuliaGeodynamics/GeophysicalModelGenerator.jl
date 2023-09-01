# 
# this is ProfileProcessing.jl
# It contains functions and type definitions to gather selected data for given profiles

export ProfileData, ExtractProfileData, CreateProfileData

# load packages
using DelimitedFiles

"""
    struct ProfileData
        lon_start::Float64
        lat_start::Float64
        lon_end::Float64
        lon_end::Float64
        VolData::GeoData
        SurfData::Vector{GeophysicalModelGenerator.GeoData}
        PointData::Vector{GeophysicalModelGenerator.GeoData}
    end

    Structure to store cross section data
"""
mutable struct ProfileData
    start_point::Tuple{Float64,Float64}
    end_point::Tuple{Float64,Float64}
    VolData::GeophysicalModelGenerator.GeoData
    SurfData::Vector{GeophysicalModelGenerator.GeoData}
    PointData::Vector{GeophysicalModelGenerator.GeoData}

    function ProfileData(;kwargs...) # this constructor allows to define only certain fields and leave the others blank
        K = new()
        for (key, value) in kwargs
            # make sure that start and end point are given as tuples of Float64
            if key==Symbol("start_point")
                setfield!(K, key, convert(Tuple{Float64,Float64},value))
            elseif key==Symbol("end_point")
                setfield!(K, key, convert(Tuple{Float64,Float64},value))
            else
                setfield!(K, key, value)
            end
        end
        return K
    end
end

### function to process volume data
function CreateProfileVolume!(Profile,DataSetName,DataSetFile,DimsVolCross,DepthVol)
    num_datasets = length(DataSetName)
    fields_vol = NamedTuple()
    local lon_vol
    local lat_vol
    local depth_vol

    println("Number of volume datasets ", num_datasets)
    for idata = 1:num_datasets
        # load data set --> each data set should have been saved in a single GeoData structure, so we'll only have to get the respective key to load the correct type
        println("processing ",DataSetFile[idata],"...")
        
        tmp_load = load(DataSetFile[idata])  # this gives us a dict with a key that is the name if the data set and the values as the GeoData structure
        tmp_load = collect(values(tmp_load))      # this gives us a vector with a single GeoData entry
        data_tmp = tmp_load[1]               # and now we extract that entry...
        cross_tmp = CrossSection(data_tmp,dims=DimsVolCross,Start=Profile.start_point,End=Profile.end_point,Depth_extent=DepthVol)        # create the cross section

        # store profile coordinates and field data on first go
        if idata==1
            # get lon,lat and depth
            # as these are in GeoUnits and our GeoData structure does not take them as input, we need to only take the value
            lon_vol = cross_tmp.lon.val
            lat_vol = cross_tmp.lat.val
            depth_vol = cross_tmp.depth.val # this will be in km
            
            # extract fields
            tmp_fields  = cross_tmp.fields;
            tmp_key     = keys(cross_tmp.fields) # get the key of all the fields

            # create a named tuple to store the fields with changed field name
            fieldname   = DataSetName[idata]*"_"*String(tmp_key[1])
            fielddata   = cross_tmp.fields[1]
            fields_vol = NamedTuple{(Symbol(fieldname),)}((fielddata,))

            for ifield = 2:length(tmp_fields) 
                fieldname   = DataSetName[idata]*"_"*String(tmp_key[ifield])
                fielddata   = cross_tmp.fields[ifield]
                new_field   = NamedTuple{(Symbol(fieldname),)}((fielddata,))
                fields_vol = merge(fields_vol,new_field) # add to the existing NamedTuple
            end
        else # only store fields
            tmp_fields  = cross_tmp.fields;
            tmp_key     = keys(cross_tmp.fields) # get the key of all the fields
            for ifield = 1:length(tmp_fields) 
                fieldname   = DataSetName[idata]*"_"*String(tmp_key[ifield])
                fielddata   = cross_tmp.fields[ifield]
                new_field   = NamedTuple{(Symbol(fieldname),)}((fielddata,))
                fields_vol = merge(fields_vol,new_field) # add to the existing NamedTuple
            end
        end

    end

    tmp = GeophysicalModelGenerator.GeoData(lon_vol,lat_vol,depth_vol,fields_vol)
    # flatten cross section and add this data to the structure
    x_profile = FlattenCrossSection(tmp,Start=Profile.start_point)
    tmp = AddField(tmp,"x_profile",x_profile)

    Profile.VolData = tmp # assign to Profile data structure
    return
end

### function to process surface data - contrary to the volume data, we here have to save lon/lat/depth pairs for every surface data set, so we create a vector of GeoData data sets
function CreateProfileSurface!(Profile,DataSetName,DataSetFile,DimsSurfCross)
    num_datasets = length(DataSetName)
    println("Number of surface datasets ", num_datasets)

    tmp = Vector{GeophysicalModelGenerator.GeoData}(undef,num_datasets)

    for idata = 1:num_datasets
        println("processing ",DataSetFile[idata],"...")
        # load data set --> each data set should have been saved in a single GeoData structure, so we'll only have to get the respective key to load the correct type
        tmp_load = load(DataSetFile[idata])  # this gives us a dict with a key that is the name if the data set and the values as the GeoData structure
        tmp_load = collect(values(tmp_load))      # this gives us a vector with a single GeoData entry
        data_tmp = tmp_load[1]               # and now we extract that entry...
        tmp[idata] = CrossSection(data_tmp, dims=DimsSurfCross,Start=Profile.start_point,End=Profile.end_point)        # create the cross section
        # flatten cross section and add this data to the structure
        x_profile = FlattenCrossSection(tmp[idata],Start=Profile.start_point)
        tmp[idata]      = AddField(tmp[idata],"x_profile",x_profile)
        # add the data set name as an attribute (not required if there is proper metadata, but odds are that there is not)
        tmp[idata].atts["dataset"] = DataSetName[idata]
    end

    Profile.SurfData = tmp # assign to profile data structure
    return 
end

### function to process point data - contrary to the volume data, we here have to save lon/lat/depth pairs for every point data set
function CreateProfilePoint!(Profile,DataSetName,DataSetFile,WidthPointProfile)
    num_datasets = length(DataSetName)
    tmp = Vector{GeophysicalModelGenerator.GeoData}(undef,num_datasets)
    println("Number of point datasets ", num_datasets)

    for idata = 1:num_datasets
        println("processing ",DataSetFile[idata],"...")
        # load data set --> each data set should have been saved in a single GeoData structure, so we'll only have to get the respective key to load the correct type
        tmp_load = load(DataSetFile[idata])  # this gives us a dict with a key that is the name if the data set and the values as the GeoData structure
        tmp_load = collect(values(tmp_load))      # this gives us a vector with a single GeoData entry
        data_tmp = tmp_load[1]               # and now we extract that entry...
        tmp[idata] = CrossSection(data_tmp,Start=Profile.start_point,End=Profile.end_point,section_width = WidthPointProfile)        # create the cross section
        # flatten cross section and add this data to the structure
        x_profile = FlattenCrossSection(tmp[idata],Start=Profile.start_point)
        tmp[idata]       = AddField(tmp[idata],"x_profile",x_profile)
        # add the data set name as an attribute (not required if there is proper metadata, but odds are that there is not)
        tmp[idata].atts["dataset"] = DataSetName[idata]
    end

    Profile.PointData = tmp # assign to profile data structure

    return 
end

### wrapper function to extract data for a single profile
function ExtractProfileData(ProfileCoordFile,ProfileNumber,DataSetName,DataSetFile,DataSetType,DimsVolCross,DepthVol,DimsSurfCross,WidthPointProfile)

    # start and end points are saved in a text file
    profile_data = readdlm(ProfileCoordFile,skipstart=1,',')

    NUM = profile_data[ProfileNumber,1]
    LON_START = profile_data[ProfileNumber,2]
    LAT_START = profile_data[ProfileNumber,3]
    LON_END   = profile_data[ProfileNumber,4]
    LAT_END   = profile_data[ProfileNumber,5]

    println("lon start ",LON_START)
    println("lat start ",LAT_START)
    println("lon end ",LON_END)
    println("lat end ",LAT_END)

    # create the cross section data set with the given lat and lon data (rest will be added later)
    Profile = ProfileData(start_point=(LON_START,LAT_START),end_point=(LON_END,LAT_END))

    # Determine the number of volume, surface and point data sets
    ind_vol    = findall( x -> x .== "Volume", DataSetType)
    ind_surf   = findall( x -> x .== "Surface", DataSetType)
    ind_point  = findall( x -> x .== "Point", DataSetType)

    # extract volume data
    CreateProfileVolume!(Profile,DataSetName[ind_vol],DataSetFile[ind_vol],DimsVolCross,DepthVol)

    # extract surface data
    CreateProfileSurface!(Profile,DataSetName[ind_surf],DataSetFile[ind_surf],DimsSurfCross)

    # extract point data
    CreateProfilePoint!(Profile,DataSetName[ind_point],DataSetFile[ind_point],WidthPointProfile)

    return Profile
end

#= 
### wrapper function to read the profile numbers+coordinates from a text file, the dataset names+locations+types from another text file
### once this is done, the different datasets are projected onto the profiles
function CreateProfileData(file_profiles,file_datasets,Depth_extent=(-300,0),DimsVolCross=(500,300),DimsSurfCross = (100,),WidthPointProfile = 20km,SaveMatlab = false)
    # get the number of profiles
    profile_data = readdlm(file_profiles,skipstart=1,',')
    NUM          = convert.(Int,profile_data[:,1]);

    ProfileNumber = NUM; # profile number, can also be a sequence of numbers

    # get dataset info
    datasets = readdlm(file_datasets,',',skipstart =1); # read information on datasets to be used from text file

    DataSetName = rstrip.(datasets[:,1]);
    DataSetFile = rstrip.(datasets[:,2]);
    DataSetType = rstrip.(datasets[:,3]);

    for iprofile = 1:length(ProfileNumber)

        # 2. process the profiles
        ExtractedData = ExtractProfileData(file_profiles,ProfileNumber[iprofile],DataSetName,DataSetFile,DataSetType,DimsVolCross,Depth_extent,DimsSurfCross,WidthPointProfile)
    
        # 3. save data 
        fn = "Profile"*string(ProfileNumber[iprofile])
        jldsave(fn*".jld2";ExtractedData)
    
    end

end
 =#

