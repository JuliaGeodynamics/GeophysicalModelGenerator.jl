# this tests creating paraview output from the data
using Test
using GeophysicalModelGenerator

@testset "DislocationModels" begin

    X,Y,Z = XYZGrid(-10:.02:10,-8:.02:8, 0:0);
    X0 = 3; Y0 = 1; depth = 5; L = 3; W = 2; plunge = 10; dip = 45;
    strike = 30; rake = 20; slip = 1; opening = 0.15; nu = 0.25;
    RefPoint = :Pc;
    ue,un,uv = RDdispSurf(X,Y,X0,Y0,depth,L,W,plunge,dip,strike,rake,slip,opening,nu,RefPoint);
    @test  sum(ue) ≈ 1758.97552383994
    @test  sum(un) ≈ 2396.5708120183926
    @test  sum(uv) ≈ 2467.229453347903
    
end