# This shows how to plot a 3D seismic tomography model 

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