using GeophysicalModelGenerator, Test, Statistics

XYZ                 =   xyz_grid(1.0:1:10.0, 11.0:1:21.0, -23:1:-10);
Dat                 =   zeros(size(XYZ[1]));
Rho                 =   ones(Float64, size(Dat))*3000;
Phases              =   zeros(Int32,   size(Dat));
Sxx                 =   XYZ[3]*10;
Stress              =   (Sxx,Sxx,Sxx,Sxx,Sxx,Sxx,Sxx,Sxx,Sxx)
Data                =   CartData(XYZ...,(Rho=Rho,Sxx=Sxx))   
Data_tuple          =   CartData(XYZ...,(Rho=Rho,Sxx=Sxx, Stress=Stress))   

fname_asagi = write_ASAGI("test", Data)
@test fname_asagi == "test_ASAGI.nc"

# Read back file:
Data_ASAGI = read_ASAGI(fname_asagi)
@test sum(Data_ASAGI.fields.Sxx - Data.fields.Sxx) == 0
@test eltype(Data_ASAGI.fields.Rho[10]) == Float64

# Read back SeisSol file:
Data_SeisSol = read_ASAGI("test_files/tpv34_rhomulambda-inner.nc")
@test mean(Data_SeisSol.fields.rho) ≈ 2635.4805f0
@test eltype(Data_SeisSol.fields.rho[10]) == Float32

# test that specifying specific field works
fname_asagi = write_ASAGI("test", Data, (:Sxx,))
Data_ASAGI2 = read_ASAGI(fname_asagi)
@test sum(Data_ASAGI2.fields.Sxx - Data.fields.Sxx) == 0

# test that it errors if we use a tuple with non-scalar fields
@test_throws "Field Stress is not an Array but instead a NTuple{9, Array{Float64, 3}}; only Arrays are supported" write_ASAGI("test", Data_tuple)

# Cleanup
foreach(rm, filter(endswith(".nc"), readdir()))
