"""
This is Tutorial_CSV.jl. It contains all the necessary commands to import a CSV file from a seismic tomography, 
to convert it to the GeoData format and to export it to a Paraview format. The following steps are performed:
1. Read data from this file.
2. Put the data in a GeoData format (this is the format that is used internally in the GMG).
3. Export that data to a format readable by Paraview.

You will need to download the data from:
[https://www.seismologie.ifg.uni-kiel.de/en/research/research-data/mere2020model](https://www.seismologie.ifg.uni-kiel.de/en/research/research-data/mere2020model).
"""

using DelimitedFiles, GeophysicalModelGenerator, GeoStats

# 1. define the filename and path
filename = "./El-Sharkawy-etal-G3.2020_MeRE2020_Mediterranean.csv" # this only works if you are in the directory where the data is located

# 2. load desired data
data        =   readdlm(filename,'|',Float64,'\n', skipstart=23,header=false)
lat         =   data[:,1];
lon         =   data[:,2];
depth       =  -data[:,3];
Vs          =   data[:,4];

# Create 3D regular grid:
Depth_vec       =   unique(depth)
Lon,Lat,Depth   =   LonLatDepthGrid(-10:0.5:40,32:0.25:50,Depth_vec);

# Employ GeoStats to interpolate irregular data points to a regular grid
dLon = Lon[2,1,1]-Lon[1,1,1]
dLat = Lat[1,2,1]-Lat[1,1,1]
Cgrid = CartesianGrid((size(Lon,1),size(Lon,2)),(minimum(Lon),minimum(Lat)),(dLon,dLat))
Vs_3D = zeros(size(Depth));
for iz=1:size(Depth,3)
    println("Depth = $(Depth[1,1,iz])")
    ind   = findall(x -> x==Depth[1,1,iz], depth)
    coord = PointSet([lon[ind]'; lat[ind]'])
    Geo   = georef((Vs=Vs[ind],), coord)
    P     = EstimationProblem(Geo, Cgrid, :Vs)
    S     = IDW(:Vs => (distance=Euclidean(),neighbors=2)); 
    sol   = solve(P, S)
    sol_Vs= values(sol).Vs
    Vs_2D = reshape(sol_Vs, size(domain(sol)))
    Vs_3D[:,:,iz] = Vs_2D;
end

# Save data to paraview:
Data_set    =   GeophysicalModelGenerator.GeoData(Lon,Lat,Depth,(Vs_km_s=Vs_3D,))   # the GeoStats package defines its own GeoData structure, so you have to choose the correct one here
Write_Paraview(Data_set, "MeRe_ElSharkawy")
