using GeophysicalModelGenerator

XYZ                 =   xyz_grid(1.0:1:10.0, 11.0:1:21.0, -23:1:-10);
Dat                 =   zeros(size(XYZ[1]));
Rho                 =   ones(Float64, size(Dat))*1350;
Phases              =   zeros(Int32,   size(Dat));
Sxx                 =   XYZ[2]*10;
Data                =   CartData(XYZ...,(Rho=Rho,Sxx=Sxx))   


write_ASAGI("test", Data)
