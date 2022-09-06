```@meta
EditURL = "<unknown>/LaPalma_example.jl"
```

# Create a 3D LaMEM setup for La Palma

## Aim
In this tutorial, your will learn how to use real data to create a geodynamic model setup with LaMEM. We will use the data of La Palma, which is a volcanic island that started erupting in mid september 2021.
LaMEM is a cartesian geodynamic model, which implies that you will have to transfer the data from `GeoData` to `CartData`.


## 1. Load data
We will use two types of data to create the model
 1) Topography
 2) Earthquake locations

We start with loading the required packages

````julia
using GeophysicalModelGenerator, JLD2
using GMT
````

In case you managed to install GMT on your machine, you can automatically download the topography with

````julia
Topo = ImportTopo(lon = [-18.7, -17.1], lat=[28.0, 29.2], file="@earth_relief_03s.grd")
````

````
GeoData
  size  : (1921, 1441, 1)
  lon   ϵ [ 341.3 : 342.9]
  lat   ϵ [ 28.0 : 29.2]
  depth ϵ [ -4.38297802734375 km : 2.414 km]
  fields: (:Topography,)

````

In case you did not manage that, we have prepared a JLD2 file [here](https://seafile.rlp.net/f/03286a46a9f040a69b0f/?dl=1).
Download it, and doublecheck that you are in the same directory as the data file with:

````julia
pwd()
````

````
"/Users/kausb/.julia/dev/GeophysicalModelGenerator/docs/src/scripts"
````

Load the data:

````julia
Topo=load("Topo_LaPalma.jld2","Topo")
````

````
GeoData
  size  : (1921, 1441, 1)
  lon   ϵ [ 341.3 : 342.9]
  lat   ϵ [ 28.0 : 29.2]
  depth ϵ [ -4.38297802734375 km : 2.414 km]
  fields: (:Topography,)

````

Next, lets load the seismicity. We have prepared a file with earthquake locations up to early November (from [https://www.ign.es/web/ign/portal/vlc-catalogo](https://www.ign.es/web/ign/portal/vlc-catalogo)).
Download that [here](https://seafile.rlp.net/f/03286a46a9f040a69b0f/?dl=1)

````julia
data_all_EQ = load("EQ_Data.jld2","data_all_EQ")
````

````
GeoData
  size  : (6045,)
  lon   ϵ [ -18.0341 : -17.6671]
  lat   ϵ [ 28.3102 : 28.8144]
  depth ϵ [ NaN km : NaN km]
  fields: (:Magnitude, :TimeSinceJan1_2021)

````

Write the data to paraview with

````julia
Write_Paraview(data_all_EQ,"data_all_EQ",PointsData=true)
Write_Paraview(Topo,"Topo")
````

````
Saved file: data_all_EQ.vtu
Saved file: Topo.vts

````

As earthquakes are point-wise data, you have to specify this.
![LaPalma_EQTopo_GeoData](../assets/img/TopoEQs_LaPalma_GeoData.png)
Note that this data is not in "easy" coordinates (see coordinate axis in the plot, where `z` is *not* pointing upwards).

## 2. Convert data
In order to create model setups, it is helpful to first transfer the data to Cartesian.
This requires us to first determine a *projection point*, that is fixed. Often, it is helpful to use the center of the topography for this. In the present example, we will center the model around La Palma itself:

````julia
proj = ProjectionPoint(Lon=-17.84, Lat=28.56)
````

````
ProjectionPoint(28.56, -17.84, 222158.69478000276, 3.1625327286383654e6, 28, true)
````

Once this is done you can convert the topographic data to the cartesian reference frame

````julia
EQ_cart   = Convert2CartData(data_all_EQ, proj);
Topo_cart = Convert2CartData(Topo, proj)
````

````
CartData
    size   : (1921, 1441, 1)
    x      ϵ [ -86.09445705828863 km : 73.67229892155609 km]
    y      ϵ [ -63.5531883197492 km : 73.28446155584604 km]
    z      ϵ [ -4.38297802734375 km : 2.414 km]
    fields : (:Topography,)

````

It is important to realize that the cartesian coordinates of the topographic grid is no longer strictly orthogonal after this conversion. You don't notice that in the current example, as the model domain is rather small.
In other cases, however, this is quite substantial (e.g., India-Asia collision zone).
LaMEM needs an orthogonal grid of topography, which we can create with:

````julia
Topo_LaMEM = CartData(XYZGrid(-70:.2:70,-60:.2:70,0));
nothing #hide
````

In a next step, the routine `ProjectCartData` projects a `GeoData` structure to a `CartData` struct

````julia
Topo_LaMEM = ProjectCartData(Topo_LaMEM, Topo, proj)
````

````
CartData
    size   : (701, 651, 1)
    x      ϵ [ -70.0 km : 70.0 km]
    y      ϵ [ -60.0 km : 70.0 km]
    z      ϵ [ -4.29541702331831 km : 2.3840784607152257 km]
    fields : (:Topography,)

````

Let's have a look at the data:

````julia
Write_Paraview(EQ_cart,"EQ_cart",PointsData=true)
Write_Paraview(Topo_LaMEM,"Topo_LaMEM")
````

````
Saved file: EQ_cart.vtu
Saved file: Topo_LaMEM.vts

````

## 3. Create LaMEM setup

In a next step, we need to read the LaMEM input file of the simulation. In particular, this will read the lateral dimensions of the grid and
the number of control volumes (elements), you want to apply in every direction.
The LaMEM input file can be downloaded [here](../scripts/LaPalma.dat).
Make sure you are in the same directory as the *.dat file & execute the following command

````julia
Grid    = ReadLaMEM_InputFile("LaPalma.dat")
````

````
LaMEM Grid:
  nel         : (64, 64, 32)
  marker/cell : (3, 3, 3)
  markers     : (192, 192, 192)
  x           ϵ [-50.0 : 50.0]
  y           ϵ [-40.0 : 40.0]
  z           ϵ [-80.0 : 15.0]

````

The `LaMEM_grid` structure contains the number of elements in every direction and the number of markers in every cell.
It also contains `Grid.X`, `Grid.Y`, `Grid.Z`, which are the coordinates of each of the markers in the 3 directions.

In a next step we need to give each of these points a `Phase` number (which is an integer, that indicates the type of the rock that point has), as well as the temperature (in Celsius).

````julia
Phases  = ones(Int32,size(Grid.X))*2;
Temp    = ones(Float64,size(Grid.X));
nothing #hide
````

In this example we set the temperature based on the depth, assuming a constant geotherm. Depth is given in kilometers

````julia
Geotherm =  30;
Temp     =  -Grid.Z.*Geotherm;
nothing #hide
````

Since we are in La Palma, we assume that temperatures are not less than 20 degrees. Moreover, in the mantle, we have the mantle adiabat (1350 C)

````julia
Temp[Temp.<20]    .=  20;
Temp[Temp.>1350]  .=  1350;
nothing #hide
````

Because of lacking data, we have pretty much no idea where the Moho is in La Palma.
Somewhat arbitrary, we assume it to be at 40 km depth, and set the rocks below that to "mantle":

````julia
ind = findall( Grid.Z .< -40);
Phases[ind] .= 3;
nothing #hide
````

Everything above the free surface is assumed to be "air"

````julia
ind     =   AboveSurface(Grid, Topo_cart);
Phases[ind] .= 0;
nothing #hide
````

And all "air" points that are below sea-level becomes "water"

````julia
ind = findall( (Phases.==0) .& (Grid.Z .< 0));
Phases[ind] .= 1;
nothing #hide
````

You can interpret the seismicity in different ways. If you assume that there are fully molten magma chambers, there should be no earthquakes within the magma chamber (as stresses are liklely to be very small).
If, however, this is a partially molten mush, extraction of melt of that mush will change the fluid pressure and may thus release build-up stresses.
We will assume the second option in our model setup.

Looking at the seismicity, there is a swarm at around 35 km depth

````julia
AddSphere!(Phases,Temp,Grid, cen=(0,0,-35), radius=5, phase=ConstantPhase(5), T=ConstantTemp(1200));
nothing #hide
````

A shallower one exists as well

````julia
AddEllipsoid!(Phases,Temp,Grid, cen=(-1,0,-11), axes=(3,3,8), StrikeAngle=225, DipAngle=45, phase=ConstantPhase(5), T=ConstantTemp(1200));
nothing #hide
````

And there might be a mid-crustal one

````julia
AddEllipsoid!(Phases,Temp,Grid, cen=(-0,0,-23), axes=(8,8,2), StrikeAngle=0, DipAngle=0, phase=ConstantPhase(5), T=ConstantTemp(1200));
nothing #hide
````

We can generate a 3D model setup, which must include the `Phases` and `Temp` arrays

````julia
Model3D = CartData(Grid, (Phases=Phases,Temp=Temp))
Write_Paraview(Model3D,"Model3D")
````

````
Saved file: Model3D.vts

````

The model setup looks like this
![LaMEM_ModelSetup_LaPalma](../assets/img/LaMEM_ModelSetup_LaPalma.png)

We can create a LaMEM marker file from the `Model3D` setup and the (cartesian) topography

````julia
Save_LaMEMTopography(Topo_cart, "Topography.txt")
Save_LaMEMMarkersParallel(Model3D)
````

````
Written LaMEM topography file: Topography.txt
Writing LaMEM marker file -> ./markers/mdb.00000000.dat

````

Next, you can use this to run the LaMEM model and visualize the model results in Paraview as well.
If you are interested in doing this, have a look at the LaMEM [wiki](https://bitbucket.org/bkaus/lamem/wiki/Home) pages.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*
