using GeophysicalModelGenerator, Test, Statistics, NCDatasets

XYZ = xyz_grid(1.0:1:10.0, 11.0:1:21.0, -23:1:-10);
Dat = zeros(size(XYZ[1]));
Rho = ones(Float64, size(Dat)) * 3000;
Phases = zeros(Int32, size(Dat));
Sxx = XYZ[3] * 10;
Stress = (Sxx, Sxx, Sxx, Sxx, Sxx, Sxx, Sxx, Sxx, Sxx)
Data = CartData(XYZ..., (Rho = Rho, Sxx = Sxx))
Data_tuple = CartData(XYZ..., (Rho = Rho, Sxx = Sxx, Stress = Stress))

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
fname_asagi = write_ASAGI("test", Data, fields = (:Sxx,))
Data_ASAGI2 = read_ASAGI(fname_asagi)
@test sum(Data_ASAGI2.fields.Sxx - Data.fields.Sxx) == 0

# test that converting to meters works
fname_asagi = write_ASAGI("test3", Data, fields = (:Sxx,), km_to_m = true)
Data_ASAGI3 = read_ASAGI(fname_asagi)
@test Data_ASAGI3.x.val[1] ≈ 1000.0

# test that it errors if we use a tuple with non-scalar fields
@test_throws "Field Stress is not an Array but instead a NTuple{9, Array{Float64, 3}}; only Arrays are supported" write_ASAGI("test", Data_tuple)

# integer fields keep their integer type
Data_int = CartData(XYZ..., (Rho = Rho, Phases = Phases))
Data_ASAGI4 = read_ASAGI(write_ASAGI("test4", Data_int))
@test eltype(Data_ASAGI4.fields.Phases[10]) == Int32

# a requested field that does not exist is reported by name
@test_throws "Field Nope is not in the data set" write_ASAGI("test", Data, fields = (:Nope,))

# Cleanup
foreach(rm, filter(endswith(".nc"), readdir()))

# fill values in NetCDF input become NaN (tomography) or an error (topography), never numbers
ncdir = mktempdir()
tomo_nc = joinpath(ncdir, "tomo.nc")
NCDataset(tomo_nc, "c") do ds
    defVar(ds, "longitude", [0.0, 1.0], ("longitude",))
    defVar(ds, "latitude", [0.0, 1.0], ("latitude",))
    defVar(ds, "depth", [0.0, 10.0], ("depth",))
    vs = Union{Missing, Float32}[4.0f0 for i in 1:2, j in 1:2, k in 1:2]
    vs[1, 1, 1] = missing
    defVar(ds, "vs", vs, ("longitude", "latitude", "depth"); fillvalue = -999.0f0)
end
Tomo = tomo_2_GeoData(tomo_nc)
@test isnan(Tomo.fields.vel[1, 1, 1])
@test Tomo.fields.vel[2, 2, 2] == 4
@test eltype(Tomo.fields.vel) == Float32

topo_nc = joinpath(ncdir, "topo.grd")
NCDataset(topo_nc, "c") do ds
    defVar(ds, "lon", [0.0, 1.0], ("lon",))
    defVar(ds, "lat", [0.0, 1.0], ("lat",))
    defVar(ds, "z", Union{Missing, Float32}[missing 1; 2 3], ("lon", "lat"); fillvalue = -999.0f0)
end
@test_throws "missing elevations" GeophysicalModelGenerator.read_topo_netcdf(topo_nc)
