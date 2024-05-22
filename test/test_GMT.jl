using Test
using GeophysicalModelGenerator, GMT

Topo = import_topo(lat=[30,31], lon=[50, 51] )
@test sum(Topo.depth.val) ≈ 2777.5705

Topo = import_topo([50,51, 30,31]);
@test sum(Topo.depth.val) ≈ 2777.5705

test_fwd =  import_GeoTIFF("test_files/length_fwd.tif", fieldname=:forward)
@test  maximum(test_fwd.fields.forward) ≈ 33.17775km

test2 =  import_GeoTIFF("test_files/UTM2GTIF.TIF")
@test   test2.fields.layer1[20,20] == 105.0
