using Test
using GeophysicalModelGenerator, GMT

@testset "GMT" begin
    Topo = ImportTopo(lat=[30,31], lon=[50, 51] )
    @test sum(Topo.depth.val) ≈ 2773.3734999999997

    Topo = ImportTopo([50,51, 30,31]);
    @test sum(Topo.depth.val) ≈ 2773.3734999999997

    test_fwd =  ImportGeoTIFF("test_files/length_fwd.tif", fieldname=:forward)
    @test  maximum(test_fwd.fields.forward) ≈ 33.17775km

    test2 =  ImportGeoTIFF("test_files/UTM2GTIF.TIF")
    @test   test2.fields.layer1[20,20] == 105.0
end