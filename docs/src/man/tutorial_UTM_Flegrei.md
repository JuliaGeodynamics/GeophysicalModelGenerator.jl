# Km-scale volcano tutorial using UTM coordinates

## Goal

Interpreting geological and geophysical data jointly in a volcano is crucial. This tutorial
visualizes available 3D data at a local volcano (Campi Flegrei caldera, Italy) using UTM coordinates.
It includes geological and geophysical data in UTM format from the following papers:

- Two shape files containing coastline and faults:
    - Vilardo, G., Ventura, G., Bellucci Sessa, E. and Terranova, C., 2013. Morphometry of the Campi Flegrei caldera (southern Italy). Journal of maps, 9(4), pp.635-640. doi:10.1080/17445647.2013.842508

- Earthquake data for two volcanic unrests, in 1983-84 and 2005-2016:
    - De Siena, L., Chiodini, G., Vilardo, G., Del Pezzo, E., Castellano, M., Colombelli, S., Tisato, N. and Ventura, G., 2017. Source and dynamics of a volcanic caldera unrest: Campi Flegrei, 1983–84. Scientific reports, 7(1), pp.1-13. doi:10.1038/s41598-017-08192-7

    - De Siena, L., Sammarco, C., Cornwell, D.G., La Rocca, M., Bianco, F., Zaccarelli, L. and Nakahara, H., 2018. Ambient seismic noise image of the structurally controlled heat and fluid feeder pathway at Campi Flegrei caldera. Geophysical Research Letters, 45(13), pp.6428-6436. doi:10.1029/2018GL078817

- Travel time tomography model:
    - Battaglia, Jean, Aldo Zollo, Jean Virieux, and Dario Dello Iacono, 2008. Merging active and passive data sets in traveltime tomography: the case study of Campi Flegrei caldera (Southern Italy). Geophysical Prospecting 56, no. 4: 555-573.  doi:10.1111/j.1365-2478.2007.00687.x

- Ambient noise tomography model:
    - Battaglia, Jean, Aldo Zollo, Jean Virieux, and Dario Dello Iacono, 2008. Merging active and passive data sets in traveltime tomography: the case study of Campi Flegrei caldera (Southern Italy). Geophysical Prospecting 56, no. 4: 555-573.  doi:10.1111/j.1365-2478.2007.00687.x


## Steps

#### 1. Download all data for region

You will need to download the zipped folder containing all files from: https://seafile.rlp.net/f/ff2c8424274c4d56b1f7/
Make sure that you are in the unzipped directory.

#### 2. Geomorphology

Load both the the shape (.shp) files contained in "./Geomorphology/*.shp" inside Paraview:

![Tutorial_Flegrei_Geomorphology](../assets/img/Flegrei_Geomorphology.png)

To reproduce it, represent the coastline as data points with black solid color and assign your favourite color map to the morphology.
Note that each block color number corresponds to a different morphology.

#### 3. Earthquakes

Now let's plot earthquake data provided as text files. Start loading the data contained in "./SeismicLocations/*.txt".
The first column gives us a temporal marker we can use to plot earthquakes in different periods

```julia
julia> using DelimitedFiles, GeophysicalModelGenerator
julia> data_80s            = readdlm("Seismicity_UTM_1983_1984.txt", '\t', skipstart=0, header=false);
julia> data_00s            = readdlm("Seismicity_UTM_2005_2016.txt", ' ', skipstart=0, header=false);
julia> data                = vcat(data_80s,data_00s)        
julia> time                = data[:,1];
julia> WE                  = data[:,2];
julia> SN                  = data[:,3];
julia> depth               = data[:,4];
```
Use the Geodesy.jl package to convert also in Latitude/Longitude - this will allow plotting the earthquakes together with the ISC data

```julia
julia> lat=zeros(length(WE));
julia> lon=zeros(length(WE));
julia> using Geodesy
julia> for i=1:length(WE)
julia>     utmz_i=UTMZ(WE[i],SN[i],depth[i],33,true)
julia>     lla_i=LLA(utmz_i,wgs84)
julia>     lat[i]=lla_i.lat
julia>     lon[i]=lla_i.lon
julia> end
```
Save in paraview with both UTM and LATLON formats - you need the CartData utility for this

```julia
julia> using GeophysicalModelGenerator
julia> EQ_Data_UTM             = CartData(WE,SN,depth,(Depth=depth,Time=time,));
julia> Write_Paraview(EQ_Data_UTM, "EQ_Flegrei_UTM", PointsData=true)
julia> EQ_Data_LATLON          = GeoData(lon,lat,depth/1e3,(Depth=depth,Time=time*y,));
julia> Write_Paraview(EQ_Data_LATLON, "EQ_Flegrei_LATLON", PointsData=true)
```
The final seismicity map looks like this - the colour scale distinguishes earthquakes of different decades. Notice the progressive migration of recent seismicity (black dots) towards East:

![Tutorial_Flegrei_seismicity](../assets/img/Flegrei_Seismicity.png)

#### 4. Velocity model

Using the Alps tutorial it is easy to create a paraview file from the Vp, Vs and Vp/Vs model in "./TravelTmeTomography/modvPS.dat".

```julia
julia> using DelimitedFiles, GeophysicalModelGenerator
julia> data            =   readdlm("modvPS.dat", '\t', Float64, skipstart=0, header=false);
julia> WE              =   data[:,1];
julia> SN              =   data[:,2];
julia> depth           =   data[:,3]*m;
julia> depthUTM        =   data[:,3];
julia> Vp              =   data[:,4];
julia> Vs              =   data[:,5];
julia> VpVs            =   data[:,6];
julia> resolution      =   (length(unique(depthUTM)),  length(unique(SN)), length(unique(WE)))
julia> dim_perm        =   [3 2 1]
julia> we              =   permutedims(reshape(WE, resolution), dim_perm);
julia> sn              =   permutedims(reshape(SN, resolution), dim_perm);
julia> depth           =   permutedims(reshape(depthUTM, resolution), dim_perm);
julia> Vp3d            =   permutedims(reshape(Vp, resolution), dim_perm);
julia> Vs3d            =   permutedims(reshape(Vs, resolution), dim_perm);
julia> Vp_Vs3d         =   permutedims(reshape(VpVs, resolution), dim_perm);
julia> Data_set_UTM    =   CartData(we, sn, depth, (vp = Vp3d * (km / s), vs = Vs3d * (km / s), vpvs = Vp_Vs3d,))
julia> Write_Paraview(Data_set_UTM, "CF_Velocity_UTM")
```
Including the Vp/Vs model in the previous Paraview file workspace:

![Tutorial_Flegrei_VpVs](../assets/img/Flegrei_VpVs.png)

#### 5. Horizontal slices of shear velocity on irregular grid

Using ambient noise you can map shear wave velocity at different depths. The models at each depth are contained in the files "./NoiseTomography/*.txt. We  read them consecutively in a "for" loop:

```julia
julia> using DelimitedFiles, GeophysicalModelGenerator, Glob
julia> list_files          = glob("*.txt");
julia> li                  = size(list_files, 1);
julia> for i = 1:li
julia>   name            = list_files[i];
julia>   name_vts        = name[1:3];
julia>   data            = readdlm(name, '\t', Float64);
julia>   WE              = data[:,1];
julia>   SN              = data[:,2];
julia>   depthUTM        = data[:,3];
julia>   Vs              = data[:,4];
```

However these models extend much further than the area imaged so far, so it is better to constrain them:

```julia
julia>   indWE           = findall(x -> 419000<=x<=435000, WE);
julia>   WE              = WE[indWE];
julia>   SN              = SN[indWE];
julia>   depthUTM        = depthUTM[indWE];
julia>   Vs              = Vs[indWE];
julia>   indSN           = findall(x -> 4514000<=x<=4528000, SN);
julia>   WE              = WE[indSN];
julia>   SN              = SN[indSN];
julia>   depthUTM        = depthUTM[indSN];
julia>   Vs              = Vs[indSN];
```
Also, nodes are irregular, hence we create a 3D lat,long regular grid from UTM - measurements are all in the same time zone:

```julia
julia>  l               = length(WE);
julia>  n_WE           = minimum(WE):100:maximum(WE);
julia>  n_SN           = minimum(SN):100:maximum(SN);
julia>  we, sn, Depth  = XYZGrid(n_WE, n_SN, depthUTM[1]);
julia>  Vs_3D          = zeros(size(Depth));

julia>  using GeoStats
julia>  Cgrid = CartesianGrid((size(we, 1), size(we, 2)), (minimum(we), minimum(sn)), (we[2,2,1] - we[1,1,1], sn[2,2,1] - sn[1,1,1]))
julia>  coord = PointSet([WE[:]'; SN[:]']);
julia>  Geo   = georef((Vs = Vs[ind],), coord);
julia>  P     = EstimationProblem(Geo, Cgrid, :Vs);
julia>  S     = IDW(:Vs => (;neighbors=2));
julia>  sol   = solve(P, S);
julia>  sol_Vs = values(sol).Vs;
julia>  Vs_2D = reshape(sol_Vs, size(domain(sol)));
julia>  Vs_3D[:,:,1] = Vs_2D;
julia>  Data_set_UTM        =   CartData(we, sn, Depth, (Vs_m_s = Vs_3D,))
julia>  Write_Paraview(Data_set_UTM, name_vts*"_UTM")
julia>  end
```
This is one of the horizontal sections created by the code in the previous model:

![Tutorial_Flegrei_VpVs](../assets/img/Flegrei_Noise.png)
