# This shows how to plot a 3D seismic tomography model 
# The paper that describes it is:
#
# Zhao, L., Paul, A., Malusà, M.G., Xu, X., Zheng, T., Solarino, S., Guillot, S., Schwartz, S., Dumont, T., Salimbeni, S., Aubert, C., Pondrelli, S., Wang, Q., Zhu, R., 2016. Continuity of the Alpine slab unraveled by high-resolution P wave tomography. Journal of Geophysical Research: Solid Earth 121, 8720–8737. doi:10.1002/2016JB013310
#
# You will need to download the file `Zhao_etal_JGR_2016_Pwave_Alps_3D_k60.txt` from:
#   https://seafile.rlp.net/d/a50881f45aa34cdeb3c0/
#
# And make sure that you are in the same directory as the data file

using DelimitedFiles, GeophysicalModelGenerator

# Load data:
data=readdlm("Zhao_etal_JGR_2016_Pwave_Alps_3D_k60.txt",' ',Float64,'\n', skipstart=0,header=false)
lon        = data[:,1];
lat        = data[:,2];
depth      = data[:,3];
dVp_perc   = data[:,4];

# Create 3D regular grid:
resolution =  (length(unique(lon)), length(unique(lat)), length(unique(depth)))
Lon          = reshape(lon,      resolution);
Lat          = reshape(lat,      resolution);
Depth        = reshape(depth,    resolution);
dVp_perc_3D  = reshape(dVp_perc, resolution);

# save paraview file
Data_set    =   GeoData(Lon,Lat,Depth,(dVp_Percentage=dVp_perc_3D,))
Write_Paraview(Data_set, "Zhao_etal_2016_dVp_percentage")

# extract cross-sections 
Data_cross  =   CrossSection(Data_set, Depth_level=-100km)  
Write_Paraview(Data_cross, "Zhao_CrossSection_100km")

Data_cross  =   CrossSection(Data_set, Lon_level=10)
Write_Paraview(Data_cross, "Zhao_CrossSection_Lon10")

Data_cross  =   CrossSection(Data_set, Lon_level=10, Interpolate=true)
Write_Paraview(Data_cross, "Zhao_CrossSection_Lon10_interpolated");

Data_cross  =   CrossSection(Data_set, Start=(1.0,39), End=(18,50))
Write_Paraview(Data_cross, "Zhao_CrossSection_diagonal")

# Extract a 3D subset of the data
Data_subset     =   ExtractSubvolume(Data_set,Lon_level=(5,12), Lat_level=(40,45))
Write_Paraview(Data_subset, "Zhao_Subset")

Data_subset_interp     =   ExtractSubvolume(Data_set,Lon_level=(5,12), Lat_level=(40,45), Interpolate=true)
Write_Paraview(Data_subset, "Zhao_Subset_interp")
