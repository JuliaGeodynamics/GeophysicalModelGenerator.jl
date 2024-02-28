using Test

pkg_dir = pkgdir(GeophysicalModelGenerator)

# test saving to file
Lon3D,Lat3D,Depth3D = LonLatDepthGrid(1.0:3:10.0, 11.0:4:20.0, (-20:5:-10)*km);
Data_set    =   GeophysicalModelGenerator.GeoData(Lon3D,Lat3D,Depth3D,(DataFieldName=Depth3D,))   
@test save_GMG(joinpath(pkg_dir,"test"),Data_set) == nothing


# loading from local file
data_local = load_GMG(joinpath(pkg_dir,"test"))
@test data_local.depth.val[20]==-15.0

# loading from local file
url  = "https://seafile.rlp.net/f/10f867e410bb4d95b3fe/?dl=1"
data_remote = load_GMG(url)
@test  data_remote.fields.MohoDepth[20] ≈ -17.99km


# loading remote data 
url  = "https://seafile.rlp.net/f/10f867e410bb4d95b3fe/?dl=1"
data_remote = download_data(url, "temp1.dat")
@test  data_remote[end-8:end] == "temp1.dat"

