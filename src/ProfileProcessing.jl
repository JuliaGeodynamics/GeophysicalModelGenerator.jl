#
# this is ProfileProcessing.jl
# It contains functions and type definitions to gather selected data for given profiles

export ProfileData, extractProfileData, createProfileData, GMG_Dataset, load_Dataset_file, combine_VolData
export extractProfileData!, readPickedProfiles
import Base: show

"""
Structure that holds profile data (interpolated/projected on the profile)

    struct ProfileData
        vertical        ::  Bool # vertical:true, horizontal:false
        start_lonlat    ::  Union{Nothing,Tuple{Float64,Float64}}
        end_lonlat      ::  Union{Nothing,Tuple{Float64,Float64}}
        depth           ::  Union{Nothing,Float64}
        VolData         ::  GeophysicalModelGenerator.GeoData
        SurfData        ::  Union{Nothing, NamedTuple}
        PointData       ::  Union{Nothing, NamedTuple}
    end

    Structure to store cross section data
"""
mutable struct ProfileData
    vertical        ::  Bool # vertical:true, horizontal:false
    start_lonlat    ::  Union{Nothing, Tuple{Float64,Float64}}
    end_lonlat      ::  Union{Nothing, Tuple{Float64,Float64}}
    depth           ::  Union{Nothing, Float64}
    VolData         ::  Union{Nothing, GeophysicalModelGenerator.GeoData}
    SurfData        ::  Union{Nothing, NamedTuple}
    PointData       ::  Union{Nothing, NamedTuple}

    function ProfileData(;kwargs...) # this constructor allows to define only certain fields and leave the others blank
        K = new(true,nothing,nothing,nothing,nothing,nothing,nothing)
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


function show(io::IO, g::ProfileData)
    if g.vertical
        println(io, "Vertical ProfileData")
        println(io, "  lon/lat    : $(g.start_lonlat)-$(g.end_lonlat) ")
    else
        println(io, "Horizontal ProfileData ")
        println(io, "  depth      : $(g.depth) ")
    end
    if !isnothing(g.VolData)
        println(io, "    VolData  : $(keys(g.VolData.fields)) ")
    end
    if !isnothing(g.SurfData)
        println(io, "    SurfData : $(keys(g.SurfData)) ")
    end
    if !isnothing(g.PointData)
        println(io, "        PointData: $(keys(g.PointData)) ")
    end

    return nothing
end

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

    Datasets = load_Dataset_file(file_datasets::String)

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
function load_Dataset_file(file_datasets::String)
    datasets    = readdlm(file_datasets,',',skipstart =1); # read information on datasets to be used from text file
    n           = size(datasets,1)

    # Deal with last column (in case it is not specified or not specified everywhere)
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
    Data = load_GMG(Datasets::Vector{GMG_Dataset})

This loads all the active datasets in `Datasets`, and returns a NamedTuple with Volume, Surface, Point, Screenshot and Topography data

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

    Data = (Volume=DataVol, Surface=DataSurf, Point=DataPoint, Screenshot=DataScreenshot, Topography=DataTopo)

    return Data
end

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
    Lon,Lat,Z  =   xyzGrid(lon1D, lat1D, z1D);
    DataSetRef =   GeoData(Lon, Lat, Z, (Temporary=Z,))

    # Loop through all datasets
    DataSet_Names = String.(keys(VolData))
    for (i,DataSet) in enumerate(VolData)
        DataSet_interp  = interpolateDataFields(DataSet, Lon,Lat,Z)
        names_fields    = String.(keys(DataSet_interp.fields))
        for (j,name) in enumerate(names_fields)
            name_new_field = DataSet_Names[i]*"_"*name # name of new field includes name of dataset
            # Note: we use ustrip here, and thereby remove the values, as the cross-section routine made problems
            DataSetRef = addField(DataSetRef,name_new_field, ustrip.(DataSet_interp.fields[j]))
        end
    end

    # remove fake field
    DataSetRef = removeField(DataSetRef, "Temporary")

    return DataSetRef
end


"""
    CreateProfileVolume!(Profile::ProfileData, VolData::AbstractGeneralGrid; DimsVolCross::NTuple=(100,100), Depth_extent=nothing)

Creates a cross-section through a volumetric 3D dataset `VolData` with the data supplied in `Profile`. `Depth_extent` can be the minimum & maximum depth for vertical profiles
"""
function CreateProfileVolume!(Profile::ProfileData, VolData::AbstractGeneralGrid; DimsVolCross::NTuple=(100,100), Depth_extent=nothing)

    if Profile.vertical
        # take a vertical cross section
        cross_tmp = crossSection(VolData,dims=DimsVolCross, Start=Profile.start_lonlat,End=Profile.end_lonlat,Depth_extent=Depth_extent)        # create the cross section

        # flatten cross section and add this data to the structure
        x_profile = flattenCrossSection(cross_tmp,Start=Profile.start_lonlat)
        cross_tmp = addField(cross_tmp,"x_profile",x_profile)

    else
        # take a horizontal cross section
        cross_tmp = crossSection(VolData, Depth_level=Profile.depth, Interpolate=true, dims=DimsVolCross)
    end

    Profile.VolData = cross_tmp # assign to Profile data structure
    return nothing
end


### internal function to process surface data - contrary to the volume data, we here have to save lon/lat/depth pairs for every surface data set, so we create a NamedTuple of GeoData data sets
function CreateProfileSurface!(Profile::ProfileData, DataSet::NamedTuple; DimsSurfCross=(100,))
    num_datasets = length(DataSet)

    tmp = NamedTuple()             # initialize empty one
    DataSetName = keys(DataSet)    # Names of the datasets
    for idata = 1:num_datasets

        # load data set --> each data set is a single GeoData structure, so we'll only have to get the respective key to load the correct type
        data_tmp = DataSet[idata]

        if Profile.vertical
            # take a vertical cross section
            data = crossSectionSurface(data_tmp, dims=DimsSurfCross, Start=Profile.start_lonlat, End=Profile.end_lonlat)        # create the cross section

            # flatten cross section and add this data to the structure
            x_profile   = flattenCrossSection(data,Start=Profile.start_lonlat)
            data        = addField(data,"x_profile",x_profile)

            # add the data set as a NamedTuple
            data_NT     = NamedTuple{(DataSetName[idata],)}((data,))
            tmp         = merge(tmp,data_NT)

        else
            # we do not have this implemented
            #error("horizontal profiles not yet implemented")
        end
    end

    Profile.SurfData = tmp # assign to profile data structure
    return
end



### function to process point data - contrary to the volume data, we here have to save lon/lat/depth pairs for every point data set
function CreateProfilePoint!(Profile::ProfileData, DataSet::NamedTuple; section_width=50km)
    num_datasets = length(DataSet)

    tmp = NamedTuple()             # initialize empty one
    DataSetName = keys(DataSet)    # Names of the datasets
    for idata = 1:num_datasets

        # load data set --> each data set is a single GeoData structure, so we'll only have to get the respective key to load the correct type
        data_tmp = DataSet[idata]

        if Profile.vertical
            # take a vertical cross section
            data    = crossSectionPoints(data_tmp, Start=Profile.start_lonlat, End=Profile.end_lonlat, section_width = section_width)        # create the cross section

            # flatten cross section and add this data to the structure
            x_profile   = flattenCrossSection(data,Start=Profile.start_lonlat)
            data        = addField(data,"x_profile",x_profile)

            # add the data set as a NamedTuple
            data_NT     = NamedTuple{(DataSetName[idata],)}((data,))
            tmp         = merge(tmp,data_NT)

        else
            # take a horizontal cross section
            data    = crossSection(data_tmp, Depth_level=Profile.depth, section_width = section_width)        # create the cross section

            # add the data set as a NamedTuple
            data_NT     = NamedTuple{(DataSetName[idata],)}((data,))
            tmp         = merge(tmp,data_NT)
        end
    end

    Profile.PointData = tmp # assign to profile data structure
    return
end


"""
    extractProfileData!(Profile::ProfileData,VolData::GeoData, SurfData::NamedTuple, PointData::NamedTuple; DimsVolCross=(100,100),Depth_extent=nothing,DimsSurfCross=(100,),section_width=50)

Extracts data along a vertical or horizontal profile
"""
function extractProfileData!(Profile::ProfileData,VolData::Union{Nothing,GeoData}=nothing, SurfData::NamedTuple=NamedTuple(), PointData::NamedTuple=NamedTuple(); DimsVolCross=(100,100),Depth_extent=nothing,DimsSurfCross=(100,),section_width=50km)

    if !isnothing(VolData)
        CreateProfileVolume!(Profile, VolData; DimsVolCross=DimsVolCross, Depth_extent=Depth_extent)
    end
    CreateProfileSurface!(Profile, SurfData, DimsSurfCross=DimsSurfCross)
    CreateProfilePoint!(Profile, PointData, section_width=section_width)

    return nothing
end

"""
This reads the picked profiles from disk and returns a vector of ProfileData
"""
function readPickedProfiles(ProfileCoordFile::String)

    profiles = Vector{ProfileData}()
    profile_data = readdlm(ProfileCoordFile,skipstart=1,',')

    for i=1:size(profile_data,1)
        start_lonlat = (profile_data[i,2:3]...,)
        end_lonlat   = (profile_data[i,4:5]...,)
        profile      = ProfileData(start_lonlat=start_lonlat, end_lonlat=end_lonlat)
        push!(profiles,profile)
    end

    return profiles
end

# this is mostly for backwards compatibility
"""
    extractProfileData(ProfileCoordFile::String,ProfileNumber::Int64,DataSetFile::String; DimsVolCross=(100,100),DepthVol=nothing,DimsSurfCross=(100,),WidthPointProfile=50km)

This is a convenience function (mostly for backwards compatibility with the MATLAB GUI) that loads the data from file & projects it onto a profile
"""
function extractProfileData(ProfileCoordFile::String,ProfileNumber::Int64,DataSetFile::String; DimsVolCross=(100,100),DepthVol=nothing,DimsSurfCross=(100,),WidthPointProfile=50km)

    # read profile
    profile_list = readPickedProfiles(ProfileCoordFile)
    profile = profile_list[ProfileNumber]

    println("lon start ",   profile.start_lonlat[1])
    println("lat start ",   profile.start_lonlat[2])
    println("lon end ",     profile.end_lonlat[1])
    println("lat end ",     profile.end_lonlat[2])

    # read all datasets:
    Datasets_all = load_Dataset_file(DataSetFile)

    # load all Data
    VolData, SurfData, PointData, ScreenshotData, TopoData = load_GMG(Datasets_all)

    # merge VolData:
    VolData_combined = combine_VolData(VolData)

    # project data onto profile:
    extractProfileData!(profile, VolData_combined, SurfData, PointData,
                        DimsVolCross=DimsVolCross, DimsSurfCross=DimsSurfCross,
                        Depth_extent=DepthVol, section_width=WidthPointProfile)

    return profile
end


#=

# Boris: I don't know exactly in which format you have your current files;
### wrapper function to extract data for a single profile
function extractProfileData(ProfileCoordFile,ProfileNumber,DataSetName,DataSetFile,DataSetType,DimsVolCross,DepthVol,DimsSurfCross,WidthPointProfile)

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
    Profile = ProfileData(start_lonlat=(LON_START,LAT_START),end_lonlat=(LON_END,LAT_END))

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

function createProfileData(file_profiles,file_datasets;Depth_extent=(-300,0),DimsVolCross=(500,300),DimsSurfCross = (100,),WidthPointProfile = 20km)
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
        ExtractedData = extractProfileData(file_profiles,ProfileNumber[iprofile],DataSetName,DataSetFile,DataSetType,DimsVolCross,Depth_extent,DimsSurfCross,WidthPointProfile)

        # 3. save data as JLD2
        fn = "Profile"*string(ProfileNumber[iprofile])
        jldsave(fn*".jld2";ExtractedData)

    end

end
=#
