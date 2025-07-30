using Test
using GeophysicalModelGenerator, GMT

Topo = import_topo(lon = [8, 9], lat = [50, 51])
@test sum(Topo.depth.val) ≈ 1076.7045 rtol = 1e-1

Topo = import_topo([8, 9, 50, 51]);
@test sum(Topo.depth.val) ≈ 1076.7045 rtol = 1e-1

test_fwd = import_GeoTIFF("test_files/length_fwd.tif", fieldname = :forward)
@test  maximum(test_fwd.fields.forward) ≈ 33.17775km

test2 = import_GeoTIFF("test_files/UTM2GTIF.TIF")
@test   test2.fields.layer1[20, 20] == 233.0
