# This shows how to  

# You will need to download the data from:
# https://www.seismologie.ifg.uni-kiel.de/en/research/research-data/mere2020model.
#
# And make sure that you are in the same directory as the data file

using DelimitedFiles, GeophysicalModelGenerator, GeoStats

# Load data:
data        =   readdlm("El-Sharkawy-etal-G3.2020_MeRE2020_Mediterranean.csv",'|',Float64,'\n', skipstart=23,header=false)
lat         =   data[:,1];
lon         =   data[:,2];
depth       =  -data[:,3];
Vs          =   data[:,4];

# Create 3D regular grid:
Depth_vec       =   unique(depth)
Lon,Lat,Depth   =   LonLatDepthGrid(-10:0.5:40,32:0.25:50,Depth_vec);

# Employ GeoStats to interpolate irregular data points to a regular grid
Cgrid = CartesianGrid((size(Lon,1),size(Lon,2)),(minimum(Lon),minimum(Lat)),(Lon[2,2,2]-Lon[1,1,1],Lat[2,2,2]-Lat[1,1,1]))
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
Data_set    =   GeoData(Lon,Lat,Depth,(Vs_km_s=Vs_3D,))   
Write_Paraview(Data_set, "MeRe_ElSharkawy")