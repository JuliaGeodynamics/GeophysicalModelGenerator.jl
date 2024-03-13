# test setting geometries in the different grid types
using Test, GeophysicalModelGenerator, GeoParams

# GeoData
Lon3D,Lat3D,Depth3D =   LonLatDepthGrid(1.0:1:10.0, 11.0:1:20.0, (-20:1:-10)*km);
Data                =   zeros(size(Lon3D));
Temp                =   ones(Float64, size(Data))*1350;
Phases              =   zeros(Int32,   size(Data));
Grid                =   GeoData(Lon3D,Lat3D,Depth3D,(DataFieldName=Data,))   

AddBox!(Phases,Temp,Grid, xlim=(2,4), zlim=(-15,-10), phase=ConstantPhase(3), DipAngle=10, T=LinearTemp(Tbot=1350, Ttop=200))
@test sum(Temp[1,1,:]) ≈ 14850.0

AddEllipsoid!(Phases,Temp,Grid, cen=(4,15,-17), axes=(1,2,3), StrikeAngle=90, DipAngle=45, phase=ConstantPhase(2), T=ConstantTemp(600))
@test sum(Temp[1,1,:]) ≈ 14850.0



# CartData 
X,Y,Z               =   XYZGrid(1.0:1:10.0, 11.0:1:20.0, -20:1:-10);
Data                =   zeros(size(X));
Temp                =   ones(Float64, size(Data))*1350;
Phases              =   zeros(Int32,   size(Data));
Grid                =   CartData(X,Y,Z,(DataFieldName=Data,))   

AddBox!(Phases,Temp,Grid, xlim=(2,4), zlim=(-15,-10), phase=ConstantPhase(3), DipAngle=10, T=LinearTemp(Tbot=1350, Ttop=200))
@test sum(Temp[1,1,:]) ≈ 14850.0

AddEllipsoid!(Phases,Temp,Grid, cen=(4,15,-17), axes=(1,2,3), StrikeAngle=90, DipAngle=45, phase=ConstantPhase(2), T=ConstantTemp(600))
@test sum(Temp[1,1,:]) ≈ 14850.0

# CartGrid
Grid                =   CreateCartGrid(size=(10,20,30),x=(0.,10), y=(0.,10), z=(2.,10))
Temp                =   ones(Float64, Grid.N...)*1350;
Phases              =   zeros(Int32,  Grid.N...);

AddBox!(Phases,Temp,Grid, xlim=(2,4), zlim=(4,8), phase=ConstantPhase(3), DipAngle=10, T=LinearTemp(Tbot=1350, Ttop=200))
@test maximum(Phases) == 3

addStripes!(Phases, Grid,stripAxes=(1,1,1),stripeWidth=0.2,stripeSpacing=1,Origin=nothing, StrikeAngle=0, DipAngle=10,phase = ConstantPhase(3),stripePhase = ConstantPhase(4))
@test maximum(Phases) == 4

# Create a CartData structure from it
Data = CartData(Grid, (T=Temp, Phases=Phases))

@test NumValue(Data.x[3,3,2]) ≈ 2.2222222222222223

# Doing the same for vertical cross-sections
Grid2D              =   CreateCartGrid(size=(10,30),x=(0.,10), z=(2.,10))
Temp2D              =   ones(Float64, Grid2D.N...)*1350;
Phases2D            =   zeros(Int32,  Grid2D.N...);

Data2D = CartData(Grid2D, (T=Temp2D, Phases=Phases2D))

@test NumValue(Data.x[3,1,2]) ≈ 2.2222222222222223


# LithosphericPhases
LP                  =   LithosphericPhases(Layers=[5 10 6], Phases=[0 1 2 3], Tlab=nothing);
X,Y,Z               =   XYZGrid(-5:1:5,-5:1:5,-20:1:5);
Phases              =   zeros(Int32,   size(X));
Temp                =   zeros(Int32,   size(X));
Phases              =   Compute_Phase(Phases, Temp, X, Y, Z, LP);

@test Phases[1,1,end]   == 3
@test Phases[1,1,7]     == 1

Phases              =   Compute_Phase(Phases, Temp, X, Y, Z, LP, Ztop=5);

@test Phases[1,1,end-4] == 0
@test Phases[1,1,5]     == 2

LP                  =   LithosphericPhases(Layers=[0.5 1.0 1.0], Phases=[0 1 2], Tlab=nothing);
Grid                =   ReadLaMEM_InputFile("test_files/SaltModels.dat");
Phases              =   zeros(Int32,   size(Grid.X));
Temp                =   zeros(Int32,   size(Grid.X));
Phases              =   Compute_Phase(Phases, Temp, Grid, LP);

@test Phases[1,1,25]    == 1
@test Phases[1,1,73]    == 0

# Create Grid & nondimensionalize it
CharDim     =   GEO_units();
Grid        =   CreateCartGrid(size=(10,20,30),x=(0.0km,10km), y=(0.0km, 10km), z=(-10.0km, 2.0km), CharDim=CharDim)
@test Grid.Δ[2] ≈ 0.0005263157894736842


# test 1D-explicit thermal solver for AddBox -----------
nel         =   96
Grid        =   CreateCartGrid(size=(nel,nel,nel),x=(-200.,200.), y=(-200.,200.), z=(-200.,0))
Temp        =   zeros(Float64, Grid.N...);
Phases      =   zeros(Int64,  Grid.N...);

# 1) horizontally layer lithosphere; UpperCrust,LowerCrust,Mantle
AddBox!(Phases,Temp,Grid, xlim=(-100,100), zlim=(-100,0), Origin=(0.0,0.0,0.0),
    phase=LithosphericPhases(Layers=[20 15 65], Phases = [1 2 3], Tlab=nothing), 
    DipAngle=0.0, T=LithosphericTemp(nz=201))

@test sum(Temp[Int64(nel/2),Int64(nel/2),:]) ≈ 36131.638045729735

# 2) inclined lithosphere; UpperCrust,LowerCrust,Mantle
Temp    =   zeros(Float64, Grid.N...);
Phases  =   zeros(Int64,  Grid.N...);

AddBox!(Phases,Temp,Grid, xlim=(-100,100), zlim=(-100,0), Origin=(0.0,0.0,0.0),
    phase=LithosphericPhases(Layers=[20 15 65], Phases = [1 2 3], Tlab=nothing), 
    DipAngle=30.0, T=LithosphericTemp(nz=201))

@test sum(Temp[Int64(nel/2),Int64(nel/2),:]) ≈ 41912.18172533137

# 3) inclined lithosphere with respect to the default origin; UpperCrust,LowerCrust,Mantle
Temp    =   zeros(Float64, Grid.N...);
Phases  =   zeros(Int64,  Grid.N...);

AddBox!(Phases,Temp,Grid, xlim=(-100,100), zlim=(-100,0),
    phase=LithosphericPhases(Layers=[20 15 65], Phases = [1 2 3], Tlab=nothing), 
    DipAngle=30.0, T=LithosphericTemp(nz=201))

@test sum(Temp[Int64(nel/2),Int64(nel/2),:]) ≈ 41316.11499878003

# 4) inclined lithosphere with only two layers
Temp    =   zeros(Float64, Grid.N...);
Phases  =   zeros(Int64,  Grid.N...);

ρM=3.0e3            # Density [ kg/m^3 ]
CpM=1.0e3           # Specific heat capacity [ J/kg/K ]
kM=2.3              # Thermal conductivity [ W/m/K ]
HM=0.0              # Radiogenic heat source per mass [H] = W/kg; [H] = [Q/rho]
ρUC=2.7e3           # Density [ kg/m^3 ]
CpUC=1.0e3          # Specific heat capacity [ J/kg/K ]
kUC=3.0             # Thermal conductivity [ W/m/K ]
HUC=617.0e-12       # Radiogenic heat source per mass [H] = W/kg; [H] = [Q/rho]

rheology = (
        # Name              = "UpperCrust",
        SetMaterialParams(;
            Phase               =   1,
            Density             =   ConstantDensity(; ρ=ρUC),
            HeatCapacity        =   ConstantHeatCapacity(; Cp=CpUC),
            Conductivity        =   ConstantConductivity(; k=kUC),
            RadioactiveHeat     =   ConstantRadioactiveHeat(; H_r=HUC*ρUC),     # [H] = W/m^3
        ),
        # Name              = "LithosphericMantle",
        SetMaterialParams(;
            Phase               =   2,
            Density             =   ConstantDensity(; ρ=ρM),
            HeatCapacity        =   ConstantHeatCapacity(; Cp=CpM),
            Conductivity        =   ConstantConductivity(; k=kM),
            RadioactiveHeat     =   ConstantRadioactiveHeat(; H_r=HM*ρM),       # [H] = W/m^3
        ),
    );

AddBox!(Phases,Temp,Grid, xlim=(-100,100), zlim=(-100,0),
    phase=LithosphericPhases(Layers=[20 80], Phases = [1 2], Tlab=nothing), 
    DipAngle=30.0, T=LithosphericTemp(rheology=rheology,nz=201))

@test sum(Temp[Int64(nel/2),Int64(nel/2),:]) ≈ 40513.969831615716

# 5) using flux lower boundary conditions
Temp    =   zeros(Float64, Grid.N...);
Phases  =   zeros(Int64,  Grid.N...);

AddBox!(Phases,Temp,Grid, xlim=(-100,100), zlim=(-100,0),
    phase=LithosphericPhases(Layers=[20 15 65], Phases = [1 2 3], Tlab=nothing), 
    DipAngle=30.0, T=LithosphericTemp(lbound="flux",nz=201))

@test sum(Temp[Int64(nel/2),Int64(nel/2),:]) ≈ 37359.648604722104


# Test the McKenzie thermal structure

# Create CartGrid struct
x        = LinRange(0.0,1200.0,64);
y        = LinRange(0.0,1200.0,64);
z        = LinRange(-660,50,64);
Cart     = CartData(XYZGrid(x, y, z));

# initialize phase and temperature matrix
Phase   = ones(Int32,size(Cart));
Temp    = ones(Float64,size(Cart))*1350;

# Create thermal structures
TsHC = HalfspaceCoolingTemp(Tsurface=20.0, Tmantle=1350, Age=120, Adiabat=0.4)
TsMK = McKenzie_subducting_slab(Tsurface = 20.0, Tmantle = 1350.0, v_cm_yr = 4.0, Adiabat = 0.0)

@test TsMK.v_cm_yr == 4.0
@test TsMK.it == 36

# Add a box with a McKenzie thermal structure

# horizontal 
Temp    = ones(Float64,size(Cart))*1350;
AddBox!(Phase, Temp, Cart; xlim=(0.0,600.0),ylim=(0.0,600.0), zlim=(-80.0, 0.0), phase = ConstantPhase(5), T=TsMK);
@test sum(Temp)  ≈ 3.518172093383281e8

# inclined slab
Temp    = ones(Float64,size(Cart))*1350;
AddBox!(Phase, Temp, Cart; xlim=(0.0,600.0),ylim=(0.0,600.0), zlim=(-80.0,0),StrikeAngle=0, DipAngle=45, phase = ConstantPhase(5), T=TsMK);
@test sum(Temp)  ≈ 3.5125017626287365e8



# horizontal slab, constant T
T_slab  = LinearWeightedTemperature(0,1,600.0,:X,ConstantTemp(1000), ConstantTemp(2000));
Temp    = ones(Float64,size(Cart))*1350;
AddBox!(Phase, Temp, Cart; xlim=(0.0,600.0),ylim=(0.0,600.0), zlim=(-80.0, 0.0), phase = ConstantPhase(5), T=T_slab);
@test sum(Temp)  ≈ 3.549127111111111e8

# horizontal slab, halfspace and McKenzie
T_slab = LinearWeightedTemperature(crit_dist=600, F1=TsHC, F2=TsMK);
Temp    = ones(Float64,size(Cart))*1350;
AddBox!(Phase, Temp, Cart; xlim=(0.0,600.0),ylim=(0.0,600.0), zlim=(-80.0, 0.0), phase = ConstantPhase(5), T=T_slab);
@test sum(Temp)  ≈ 3.499457641038468e8


Data_Final =   AddField(Cart,"Temp",Temp)


# test polygon structure

x        = LinRange(0.0,1200.0,64);
y        = LinRange(0.0,1200.0,64);
z        = LinRange(-660,50,64);
Cart     = CartData(XYZGrid(x, y, z));

# initialize phase and temperature matrix
Phase   = ones(Int32,(length(x),length(y),length(z)));
Temp    = ones(Float64,(length(x),length(y),length(z)))*1350;

AddBox!(Phase, Temp, Cart; xlim=(0.0,600.0),ylim=(0.0,600.0), zlim=(-80.0, 0.0), phase = ConstantPhase(5), T=T=ConstantTemp(120.0));

# add accretionary prism 
addPolygon!(Phase, Temp, Cart; xlim=[500.0, 200.0, 500.0],ylim=[100.0,400.0], zlim=[0.0,0.0,-60.0], phase = ConstantPhase(8), T=LinearTemp(Ttop=20, Tbot=30))

@test maximum(Phase) == 8
@test minimum(Temp) == 21.40845070422536
@test sum(Phase) == 292736

# Test the Bending slab geometry

# Create CartGrid struct
x        = LinRange(0.0,1200.0,128);
y        = LinRange(0.0,1200.0,128);
z        = LinRange(-660,50,128);
Cart     = CartData(XYZGrid(x, y, z));
X,Y,Z    = XYZGrid(x, y, z);

# initialize phase and temperature matrix
Phase   = ones(Int32,size(Cart));
Temp    = fill(1350.0,size(Cart));

t1 = Trench(Start = (400.0,400.0), End = (800.0,800.0),θ_max = 45, direction = 1.0, n_seg = 50, Length = 600.0, Thickness = 80.0, Lb = 500.0,d_decoupling = 100.0, type_bending =:Ribe)
@test t1.θ_max == 45.0
@test t1.Thickness == 80.0
@test t1.Length == 600.0
@test t1.Lb == 500.0

phase = LithosphericPhases(Layers=[5 7 88], Phases = [2 3 4], Tlab=nothing)
TsHC = HalfspaceCoolingTemp(Tsurface=20.0, Tmantle=1350, Age=30, Adiabat=0.4)
temp = TsHC;

addSlab!(Phase,Temp,Cart, t1, phase=phase, T = TsHC)
@test Temp[84,84,110]  ≈ 1045.1322688510577
@test extrema(Phase) == (1, 4)

# with weak zone
t1 = Trench(Start = (400.0,400.0), End = (800.0,800.0),θ_max = 45, direction = 1.0, n_seg = 50, Length = 600.0, Thickness = 80.0, Lb = 500.0,d_decoupling = 100.0, WeakzoneThickness=10, WeakzonePhase=9)
Phase   = ones(Int32,size(Cart));
Temp    = fill(1350.0,size(Cart));
addSlab!(Phase,Temp,Cart, t1, phase=phase, T = TsHC)
@test extrema(Phase) == (1, 9)

#Data_Final      =   CartData(X,Y,Z,(Phase=Phase,Temp=Temp)) 
#Write_Paraview(Data_Final, "Data_Final");

Phase = ones(Int32,size(Cart));
Temp = fill(1350.0,size(Cart));
TsMK = McKenzie_subducting_slab(Tsurface = 20.0, Tmantle = 1350.0, v_cm_yr = 4.0, Adiabat = 0.0)
temp = TsMK 

Phase = ones(Int32,size(Cart));
Temp = fill(1350.0,size(Cart));
TsHC = HalfspaceCoolingTemp(Tsurface=20.0, Tmantle=1350, Age=120, Adiabat=0.4)
TsMK = McKenzie_subducting_slab(Tsurface = 20.0, Tmantle = 1350.0, v_cm_yr = 4.0, Adiabat = 0.0)
T_slab = LinearWeightedTemperature(crit_dist=600, F1=TsHC, F2=TsMK);
phase = LithosphericPhases(Layers=[5 7 88], Phases = [2 3 4], Tlab=nothing)
t1 = Trench(Start = (400.0,400.0), End = (800.0,800.0),θ_max = 90.0, direction = 1.0, n_seg = 50, Length = 600.0, Thickness = 80.0, Lb = 500.0,d_decoupling = 100.0, type_bending =:Ribe,   WeakzoneThickness=10, WeakzonePhase=9)

addSlab!(Phase,Temp,Cart, t1, phase=phase, T = T_slab)
@test Temp[84,84,110]  ≈ 624.6682008876219

Data_Final      =   CartData(X,Y,Z,(Phase=Phase,Temp=Temp)) 

# 2D slab:
nx,nz = 512,128
x = range(-1000,1000, nx);
z = range(-660,0,    nz);
Grid2D = CartData(XYZGrid(x,0,z))
Phases = zeros(Int64, nx, 1, nz);
Temp = fill(1350.0, nx, 1, nz);
AddBox!(Phases, Temp, Grid2D; xlim=(-800,0.0), zlim=(-80.0, 0.0), phase = ConstantPhase(1),  T=HalfspaceCoolingTemp(Age=40));    

trench = Trench(Start=(0.0,-100.0), End=(0.0,100.0), Thickness=80.0, θ_max=30.0, Length=300, Lb=150);
addSlab!(Phases, Temp, Grid2D, trench, phase = ConstantPhase(2), T=HalfspaceCoolingTemp(Age=40));

T_slab = LinearWeightedTemperature( F1=HalfspaceCoolingTemp(Age=40), F2=McKenzie_subducting_slab(Tsurface=0,v_cm_yr=4, Adiabat = 0.0), crit_dist=600)
addSlab!(Phases, Temp, Grid2D, trench, phase = ConstantPhase(2), T=T_slab);

@test sum(Temp) ≈ 8.571402268095453e7
@test extrema(Phases) == (0, 2)

# Add them to the `CartData` dataset:
#Grid2D = CartData(Grid2D.x.val, Grid2D.y.val, Grid2D.z.val ,(;Phases, Temp))
#Write_Paraview(Grid2D,"Grid2D_SubductionCurvedMechanical");


# More sophisticated 2D example with overriding plate
nx,nz = 512,128
x = range(-1000,1000, nx);
z = range(-660,0,    nz);
Grid2D = CartData(XYZGrid(x,0,z))
Phases = zeros(Int64, nx, 1, nz);
Temp = fill(1350.0, nx, 1, nz);
lith = LithosphericPhases(Layers=[15 20 55], Phases=[3 4 5], Tlab=1250)

# Lets add the overriding plate. Note that we add this twice with a different thickness to properly represent the transition around the trench
AddBox!(Phases, Temp, Grid2D; xlim=(200,1000), zlim=(-150.0, 0.0), phase = lith, T=HalfspaceCoolingTemp(Age=80));
AddBox!(Phases, Temp, Grid2D; xlim=(0,200), zlim=(-60.0, 0.0), phase = lith, T=HalfspaceCoolingTemp(Age=80));

# The horizontal part of the oceanic plate is as before
v_spread_cm_yr = 3      #spreading velocity
lith = LithosphericPhases(Layers=[15 55], Phases=[1 2], Tlab=1250)
AddBox!(Phases, Temp, Grid2D; xlim=(-800,0.0), zlim=(-150.0, 0.0), phase = lith, T=SpreadingRateTemp(SpreadingVel=v_spread_cm_yr));

# Yet, now we add a trench as well. 
AgeTrench_Myrs = 800e3/(v_spread_cm_yr/1e2)/1e6    #plate age @ trench

# We want to add a smooth transition from a halfspace cooling thermal profile to a slab that is heated by the surrounding mantle below a decoupling depth `d_decoupling`.
T_slab = LinearWeightedTemperature( F1=HalfspaceCoolingTemp(Age=AgeTrench_Myrs), F2=McKenzie_subducting_slab(Tsurface=0,v_cm_yr=v_spread_cm_yr, Adiabat = 0.0), crit_dist=600)

# # in this case, we have a more reasonable slab thickness: 
trench = Trench(Start=(0.0,-100.0), End=(0.0,100.0), Thickness=90.0, θ_max=30.0, Length=600, Lb=200, 
                 WeakzoneThickness=15, WeakzonePhase=6, d_decoupling=125);
addSlab!(Phases, Temp, Grid2D, trench, phase = lith, T=T_slab);

# Lithosphere-asthenosphere boundary:
ind = findall(Temp .> 1250 .&& (Phases.==2 .|| Phases.==5));
Phases[ind] .= 0;

@test sum(Temp) ≈ 8.292000736425713e7  
@test extrema(Phases) == (0, 6)
#Grid2D = CartData(Grid2D.x.val,Grid2D.y.val,Grid2D.z.val, (;Phases, Temp))
#Write_Paraview(Grid2D,"Grid2D_SubductionCurvedOverriding");

