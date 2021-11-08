"""
This is Tutorial_Flegrei.jl. It contains all the necessary commands to import the data files from
geophysical models at the volcano and represent them in cartesian and UTM coordinates.
The following steps are performed:
1. Export earthquake data with a time label in both formats. 
2. Export a velocity models in both formats. 
3. Scale and interpolate a shear-wave model and export in both formats.

You will need to download the zipped folder containing all files from:
[https://seafile.rlp.net/f/ff2c8424274c4d56b1f7/](https://ngdc.noaa.gov/mgg/global/global.html)
Remember to work in the downloaded directory.
"""
using DelimitedFiles, GeophysicalModelGenerator, Glob, GeoStats

# 1. define where the two files are located on your computer
data_80s            = readdlm("SeismicLocations/Seismicity_UTM_1983_1984.txt", '\t', skipstart=0, header=false);
data_00s            = readdlm("SeismicLocations/Seismicity_UTM_2005_2016.txt", ' ', skipstart=0, header=false);

# 1.1 create a single file and detine the timing and locations
data                = vcat(data_80s,data_00s)        
time                = data[:,1];
WE                  = data[:,2];
SN                  = data[:,3];
depth               = data[:,4];

# 1.2 save in both formats
EQ_Data_Cart        = CartData(WE,SN,depth,(Depth=depth,Time=time*y,));
Write_Paraview(EQ_Data_Cart, "CF_Earthquakes_Cartesian", PointsData=true)
EQ_Data_UTM         = UTMData(WE, SN, depth, 33, true, (Depth=depth,Time=time,));
Data_set_UTM        =   convert(GeoData,EQ_Data_UTM)
Write_Paraview(Data_set_UTM, "CF_Earthquakes_UTM", PointsData=true)

# 2 load the velocity model and assign the locations of nodes and velocities + Vp/Vs
data            =   readdlm("TravelTimeTomography/modvPS.dat", '\t', Float64, skipstart=0, header=false);
WE              =   data[:,1];
SN              =   data[:,2];
depth           =   data[:,3]*m;
depthUTM        =   data[:,3];
Vp              =   data[:,4];
Vs              =   data[:,5];
VpVs            =   data[:,6];

# 2 Reshape everything and save in both formats.
resolution      =   (length(unique(depthUTM)),  length(unique(SN)), length(unique(WE)))
dim_perm        =   [3 2 1]
we              =   permutedims(reshape(WE, resolution), dim_perm);
sn              =   permutedims(reshape(SN, resolution), dim_perm);
depth           =   permutedims(reshape(depthUTM, resolution), dim_perm);
Vp3d            =   permutedims(reshape(Vp, resolution), dim_perm);
Vs3d            =   permutedims(reshape(Vs, resolution), dim_perm);
Vp_Vs3d         =   permutedims(reshape(VpVs, resolution), dim_perm);
Data_set_Cartesian  =   CartData(we, sn, depth, (vp = Vp3d * (km / s), vs = Vs3d * (km / s), vpvs = Vp_Vs3d,))
Write_Paraview(Data_set_Cartesian, "CF_Velocity_Cartesian")
Data_set        =   UTMData(we, sn, depth, 33, true, (vp = Vp3d * (km / s), vs = Vs3d * (km / s), vpvs = Vp_Vs3d,))
Data_set_UTM    =   convert(GeoData,Data_set)
Write_Paraview(Data_set_UTM, "CF_Velocity_UTM")

# 3 Creae a list of text files containing sections and loop through them
list_files        = glob("AmbientNoiseTomography/*.txt");
li                = size(list_files, 1);
for i = 1:li
  name            = list_files[i];
  name_vts        = name[1:3];
  data            = readdlm(name, '\t', Float64);
  WE              = data[:,1];
  SN              = data[:,2];
  depthUTM        = data[:,3];
  Vs              = data[:,4];

# 3.1 Constrain them on a smaller grid:
  indWE           = findall(x -> 419000<=x<=435000, WE);
  WE              = WE[indWE];
  SN              = SN[indWE];
  depthUTM        = depthUTM[indWE];
  Vs              = Vs[indWE];
  indSN           = findall(x -> 4514000<=x<=4528000, SN);
  WE              = WE[indSN];
  SN              = SN[indSN];
  depthUTM        = depthUTM[indSN];
  Vs              = Vs[indSN];

  # 3.2 Make the grid regular and interpolate, then write the result:
l                 = length(WE);
 n_WE             = minimum(WE):100:maximum(WE);
 n_SN             = minimum(SN):100:maximum(SN);
 we, sn, Depth    = XYZGrid(n_WE, n_SN, depthUTM[1]);
 Vs_3D            = zeros(size(Depth));
 Cgrid            = CartesianGrid((size(we, 1), size(we, 2)), (minimum(we), minimum(sn)), (we[2,2,1] - we[1,1,1], sn[2,2,1] - sn[1,1,1]))
 coord            = PointSet([WE[:]'; SN[:]']);
 Geo              = georef((Vs = Vs[ind],), coord);
 P                = EstimationProblem(Geo, Cgrid, :Vs);
 S                = IDW(:Vs => (;neighbors=2));
 sol              = solve(P, S);
 sol_Vs           = values(sol).Vs;
 Vs_2D            = reshape(sol_Vs, size(domain(sol)));
 Vs_3D[:,:,1]     = Vs_2D;
 Data_set_Cartesian =   CartData(we, sn, Depth, (Vs_m_s = Vs_3D,))
 Write_Paraview(Data_set_Cartesian, "CF_Noise_"*name_vts*"_Cartesian")
 Data_set         =   UTMData(we, sn, Depth, 33, true, (Vs = Vs_3D*(km / s),));
 Data_set_UTM     =   convert(GeoData,Data_set);
 Write_Paraview(Data_set_UTM, "CF_Noise_"*name_vts*"_UTM")
 end
