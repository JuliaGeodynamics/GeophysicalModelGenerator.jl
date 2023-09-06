# 
# this is ProfileProcessing.jl
# It contains functions and type definitions to gather selected data for given profiles

export ProfileData, ExtractProfileData, CreateProfileData, GMG_Dataset, Load_Dataset_file, combine_VolData
import Base: show

"""
Structure that holds profile data (interpolated/projected on the profile) 

    struct ProfileData
        vertical        ::  Bool # vertical:true, horizontal:false
        start_lonlat    ::  Union{Nothing,Tuple{Float64,Float64}}
        end_lonlat      ::  Union{Nothing,Tuple{Float64,Float64}}
        depth           ::  Union{Nothing,Float64}
        VolData         ::  GeophysicalModelGenerator.GeoData
        SurfData        ::  Vector{GeophysicalModelGenerator.GeoData}
        PointData       ::  Vector{GeophysicalModelGenerator.GeoData}
        ScreenshotData  ::  Union{Nothing,Vector{GeophysicalModelGenerator.GeoData}}
    end

    Structure to store cross section data
"""
mutable struct ProfileData
    vertical        ::  Bool # vertical:true, horizontal:false
    start_lonlat    ::  Union{Nothing, Tuple{Float64,Float64}}
    end_lonlat      ::  Union{Nothing, Tuple{Float64,Float64}}
    depth           ::  Union{Nothing, Float64}
    VolData         ::  Union{Nothing, GeophysicalModelGenerator.GeoData}
    SurfData        ::  Union{Nothing, Vector{GeophysicalModelGenerator.GeoData}}
    PointData       ::  Union{Nothing, Vector{GeophysicalModelGenerator.GeoData}}
    ScreenshotData  ::  Union{Nothing, Vector{GeophysicalModelGenerator.GeoData}}

    function ProfileData(;kwargs...) # this constructor allows to define only certain fields and leave the others blank
        K = new()
        for (key, value) in kwargs
            # make sure that start and end point are given as tuples of Float64
            if key==Symbol("start_lonlat")
                setfield!(K, key, convert(Tuple{Float64,Float64},Float64.(value)))
                setfield!(K, :vertical, true)
            elseif key==Symbol("end_lonlat")
                setfield!(K, key, convert(Tuple{Float64,Float64},Float64.(value)))
                setfield!(K, :vertical, true)
            elseif key==Symbol("depth")
                setfield!(K, key, convert(Float64,value))
                setfield!(K, :vertical, false)
            else
                setfield!(K, key, value)
            end
        end
        return K
    end
end


##### THE STRUCUTRE BELOW MAY BE REPLACED WITH E.G. THE NAMED TUPLE THAT IS USED IN 
# DONE - Boris
#=
"""
    struct DataSet
        Name::Vector{String}
        Type::Vector{String}
        Location::Vector{String}
        GeoData::Vector{GeophysicalModelGenerator.GeoData}
    end

    Structure to store all datasets
"""
mutable struct GeoDataSet
    Name::Vector{String}
    Type::Vector{String}
    Location::Vector{String}
    GeoData::Vector{GeophysicalModelGenerator.GeoData}

    function DataSet(num_datasets;kwargs...) # this constructor allows to define only certain fields and leave the others blank
        K = new()
        Name                = Vector{GeophysicalModelGenerator.GeoData}(undef,num_datasets)
        Type                = Vector{GeophysicalModelGenerator.GeoData}(undef,num_datasets)
        Location            = Vector{GeophysicalModelGenerator.GeoData}(undef,num_datasets)
        tmp                 = Vector{GeophysicalModelGenerator.GeoData}(undef,num_datasets)
        for (key, value) in kwargs
                setfield!(K, key, value)
        end
        return K
    end
end
=#


"""

Structure that stores info about a GMG Dataset, which is useful to collect a wide variety of datasets.

- Name    :: String          # Name of the dataset
- Type    :: String          # Volumetric, Surface, Point, Screenshot    
- DirName :: String          # Directory name or url of dataset 
- active  :: Bool            # should this data be loaded or not?

"""
mutable struct GMG_Dataset
    Name    :: String          # Name of the dataset
    Type    :: String          # Volumetric, Surface, Point, Screenshot    
    DirName :: String          # Directory name or url of dataset 
    active  :: Bool            # active in the GUI or not?

    function GMG_Dataset(Name::String,Type::String,DirName::String,active::Bool=false) 
        Type = strip(Type)
        Name = strip(Name)
        DirName = strip(DirName)
        
        if !any(occursin.(Type,["Volume","Surface","Point","Screenshot","Topography"]))
            error("Type should be either: Volume,Surface,Point,Topography or Screenshot. Is: $Type.")
        end
    
        if DirName[end-4:end] == ".jld2"
            DirName = DirName[1:end-5]
        end
        new(Name,Type,DirName,active)
    end

end


# Print info 
function show(io::IO, g::GMG_Dataset)
    if g.active
        str_act = "(active)  :"
    else
        str_act = "(inactive):"
    end    
    print(io, "GMG $(g.Type) Dataset $str_act $(g.Name) @ $(g.DirName)")
    
    return nothing
end


"""
    data::NamedTuple = load_GMG(data::GMG_Dataset)

Loads a dataset specified in GMG_Dataset `data` and returns it as a named tuple
"""
function load_GMG(data_input::GMG_Dataset)
    data = load_GMG(data_input.DirName)
    name = Symbol(data_input.Name)
    return NamedTuple{(name,)}((data,))
end

"""

    Datasets = Load_Dataset_file(file_datasets::String)

This loads a CSV textfile that lists datasets, which is expected to have the following format:

- `Name`,`Location`,`Type`, `[Active]`
-  AlpArray,./Seismicity/ALPARRAY/AlpArraySeis.jld2,Point, true
-  Plomerova2022,https://seafile.rlp.net/f/abccb8d3302b4ef5af17/?dl=1,Volume
Note that the first line of the file is skipped.

Here, the meaning of the variables is:
- `Name`: The name of the dataset to be loaded
- `Location`: the location of the file (directory and filename) on your local machine, or an url where we can download the file from the web. The url is expected to start with "http". 
- `Type`: type of the dataset (Volume, Surface, Point, Screenshot)
- `Active`: Do we want this file to be loaded or not? Optional parameter that defaults to `true`

"""
function Load_Dataset_file(file_datasets::String)
    datasets    = readdlm(file_datasets,',',skipstart =1); # read information on datasets to be used from text file
    n           = size(datasets,1)

    # Deal with last collumn (in case it is not specified or not specified everywhere)
    if size(datasets,2)==4
        active      = datasets[:,4]
        active      = replace(active,""=>true)
        active      = Bool.(active)
    elseif size(datasets,2)==3
        active      = ones(Bool,n)
    end

    Datasets = Vector{GMG_Dataset}()
    for i=1:n
        push!(Datasets, GMG_Dataset( String(datasets[i,1]), String(datasets[i,3]), String(datasets[i,2]), active[i]))
    end

    return Datasets
end

""" 
    DataVol, DataSurf, DataPoint, DataScreenshot, DataTopo = load_GMG(Datasets::Vector{GMG_Dataset})

This loads all the active datasets in `Datasets`, and returns named tuples with Volumetric, Surface, Point, Screenshot or Topography data

"""
function load_GMG(Datasets::Vector{GMG_Dataset})

    DataPoint       =   NamedTuple();
    DataSurf        =   NamedTuple();
    DataScreenshot  =   NamedTuple();
    DataVol         =   NamedTuple();
    DataTopo        =   NamedTuple();
    for data in Datasets
        if data.active
            # load into NamedTuple (I'm sure this can be done more compact somehow..)
            loaded_data = load_GMG(data)   
            if data.Type=="Volume"
                DataVol         =   merge(DataVol,loaded_data)
            elseif data.Type=="Surface"
                DataSurf        =   merge(DataSurf,loaded_data)
            elseif data.Type=="Point"
                DataPoint       =   merge(DataPoint,loaded_data)
            elseif data.Type=="Screenshot"
                DataScreenshot  =   merge(DataScreenshot,loaded_data)
            elseif data.Type=="Topography"
                DataTopo        =   merge(DataTopo,loaded_data)
            end
        end
    end

    return DataVol, DataSurf, DataPoint, DataScreenshot, DataTopo
end

#=

# BORIS: I think the functionality is now implemented in load_GMG(data::Vector{GMG_Dataset}) above

### function to load all datasets and put them in a GeoData vector
function MergeDataSets()
    # get dataset info
    datasets    = readdlm(file_datasets,',',skipstart =1); # read information on datasets to be used from text file
    DataSetName = rstrip.(datasets[:,1]);
    DataSetFile = rstrip.(datasets[:,2]);
    DataSetType = rstrip.(datasets[:,3]);

    num_datasets = length(DataSetName)

    # preallocate data set
    DataSet = GeoDataSet ## DEFINE A CONSTRUCTOR FIRST ABOVE!!!
    # loop over all datasets and load them in a single Vector
    for idata = 1:num_datasets
        println("processing ",DataSetFile[idata],"...")
        tmp_load = load(DataSetFile[idata])  # this gives us a dict with a key that is the name if the data set and the values as the GeoData structure
        tmp_load = collect(values(tmp_load))      # this gives us a vector with a single GeoData entry
        data_tmp = tmp_load[1]               # and now we extract that entry...

        ### PUT THE DATA SET IN THE GLOBAL DATA SET STRUCTURE
    end

    return DataSet
end
=#

"""

    VolData_combined = combine_VolData(VolData::NamedTuple; lat=nothing, lon=nothing, depth=nothing, dims=(100,100,100), dataset_preferred = 1)

This takes different volumetric datasets (specified in `VolData`) & merges them into a single one. 
You need to either provide the "reference" dataset within the NamedTuple (`dataset_preferred`), or the lat/lon/depth and dimensions of the new dataset.

"""
function combine_VolData(VolData::NamedTuple; lat=nothing, lon=nothing, depth=nothing, dims=(100,100,100), dataset_preferred = 1)

    # Get dimensions of new Data_set
    i = dataset_preferred
    if isnothing(lon);   lon   = extrema(VolData[i].lon.val);   end
    if isnothing(lat);   lat   = extrema(VolData[i].lat.val);   end
    if isnothing(depth); depth = extrema(VolData[i].depth.val); end
    if isnothing(dims);  dims  = size(VolData[i].depth.val);    end

    # Create reference dataset
    lon1D   = range(lon...,     dims[1])
    lat1D   = range(lat...,     dims[2])
    z1D     = range(depth...,   dims[3])
    Lon,Lat,Z  =   XYZGrid(lon1D, lat1D, z1D);
    DataSetRef =   GeoData(Lon, Lat, Z, (Temporary=Z,))

    # Loop through all datasets
    DataSet_Names = String.(keys(VolData))
    for (i,DataSet) in enumerate(VolData)
        DataSet_interp  = InterpolateDataFields(DataSet, Lon,Lat,Z)
        names_fields    = String.(keys(DataSet_interp.fields))
        for (j,name) in enumerate(names_fields)
            name_new_field = DataSet_Names[i]*"_"*name # name of new field includes name of dataset
            # Note: we use ustrip here, and thereby remove the values, as the cross-section routine made problems 
            DataSetRef = AddField(DataSetRef,name_new_field, ustrip.(DataSet_interp.fields[j]))
        end
    end

    # remove fake field
    DataSetRef = RemoveField(DataSetRef, "Temporary")

    return DataSetRef
end


"""
    CreateProfileVolume!(Profile::ProfileData, VolData::AbstractGeneralGrid; DimsVolCross::NTuple=(100,100), DepthVol=nothing)

Creates a cross-section through a volumetric 3D dataset `VolData` with the data supplied in `Profile` 
"""
function CreateProfileVolume!(Profile::ProfileData, VolData::AbstractGeneralGrid; DimsVolCross::NTuple=(100,100), DepthVol=nothing)

    if Profile.vertical 
        # take a vertical cross section
        cross_tmp = CrossSection(VolData,dims=DimsVolCross, Start=Profile.start_lonlat,End=Profile.end_lonlat,Depth_extent=DepthVol)        # create the cross section

        # flatten cross section and add this data to the structure
        x_profile = FlattenCrossSection(cross_tmp,Start=Profile.start_lonlat)
        cross_tmp = AddField(cross_tmp,"x_profile",x_profile)

    else
        # take a horizontal cross section
        cross_tmp = CrossSection(VolData, Depth_level=Profile.depth, Interpolate=true, dims=DimsVolCross)
    end
 
    Profile.VolData = cross_tmp # assign to Profile data structure
    return nothing
end




#=
### function to process volume data --> datasets are given as GeoDataSet
function CreateProfileVolume!(Profile,DataSet,DataSetFile,DimsVolCross,DepthVol)
    num_datasets = length(DataSet)
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
        
        if Profile.profile_type # if this field is true, we create a vertical profile
            # take a vertical cross section
            cross_tmp = CrossSection(data_tmp,dims=DimsVolCross,Start=Profile.start_point,End=Profile.end_point,Depth_extent=DepthVol)        # create the cross section
        else
            error("horizontal profiles not yet implemented")
        end
        
        # IT MAY BE THAT THE LINES BELOW EQUALLY WORK FOR HORIZONTAL PROFILES, SO THEY ARE NOT IN THE CONDITIONAL LOOP ABOVE
        # THIS HAS TO BE TESTED THOUGH

        # store profile coordinates and field data on first go
        if idata==1
            # get lon,lat and depth
            # as these are in GeoUnits and our GeoData structure does not take them as input, we need to only take the value
            lon_vol = cross_tmp.lon.val
            lat_vol = cross_tmp.lat.val
            depth_vol = cross_tmp.depth.val # this will be in km
            
            # MERGE FIELDS
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
=#


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

        if profile_type # if this field is true, we create a vertical profile
            # take a vertical cross section
            tmp[idata] = CrossSection(data_tmp, dims=DimsSurfCross,Start=Profile.start_point,End=Profile.end_point)        # create the cross section
            # flatten cross section and add this data to the structure
            x_profile = FlattenCrossSection(tmp[idata],Start=Profile.start_point)
            tmp[idata]      = AddField(tmp[idata],"x_profile",x_profile)
            # add the data set name as an attribute (not required if there is proper metadata, but odds are that there is not)
            tmp[idata].atts["dataset"] = DataSetName[idata]
        else
            error("horizontal profiles not yet implemented")
        end

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

### wrapper function to read the profile numbers+coordinates from a text file, the dataset names+locations+types from another text file
### once this is done, the different datasets are projected onto the profiles

### currently, the function is quite slow, as the different data sets are reloaded for each profile. 
### a faster way would be to load one data set and create the profiles from it and then move on to the next one. However, this would require to hold all the profile data in memory, which may be a bit much...

function CreateProfileData(file_profiles,file_datasets;Depth_extent=(-300,0),DimsVolCross=(500,300),DimsSurfCross = (100,),WidthPointProfile = 20km)
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
    
        # 3. save data as JLD2
        fn = "Profile"*string(ProfileNumber[iprofile])
        jldsave(fn*".jld2";ExtractedData)
    
    end

end


