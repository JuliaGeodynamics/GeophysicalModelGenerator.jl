using GeophysicalModelGenerator, Test

XYZ                 =   xyz_grid(1.0:1:10.0, 11.0:1:21.0, -23:1:-10);
Dat                 =   zeros(size(XYZ[1]));
Rho                 =   ones(Float64, size(Dat))*3000;
Phases              =   zeros(Int32,   size(Dat));
Sxx                 =   XYZ[3]*10;
Data                =   CartData(XYZ...,(Rho=Rho,Sxx=Sxx))   

fname_asagi = write_ASAGI("test", Data)
@test fname_asagi == "test_ASAGI.nc"

# Read back file:
Data_ASAGI = read_ASAGI(fname_asagi)
@test sum(Data_ASAGI.fields.Sxx - Data.fields.Sxx) == 0
