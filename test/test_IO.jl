using Test, Downloads

pkg_dir = pkgdir(GeophysicalModelGenerator)

# Guards the downloads: skipped on push/PR when the file server is unreachable, required on the scheduled run.
url = "https://seafile.rlp.net/f/10f867e410bb4d95b3fe/?dl=1"
remote_ok = try
    Downloads.request(url; method = "HEAD", timeout = 20)
    true
catch e
    e isa Downloads.RequestError || rethrow()
    get(ENV, "GITHUB_EVENT_NAME", "") == "schedule"
end

dir = mktempdir()

# test saving to file
Lon3D, Lat3D, Depth3D = lonlatdepth_grid(1.0:3:10.0, 11.0:4:20.0, (-20:5:-10) * km);
Data_set = GeophysicalModelGenerator.GeoData(Lon3D, Lat3D, Depth3D, (DataFieldName = Depth3D,))
@test save_GMG("test", Data_set; dir) == nothing

# loading from local file
data_local = load_GMG("test", dir)
@test data_local.depth.val[20] == -15.0

# a local path containing "http" is not a URL
mkpath(joinpath(dir, "http_data"))
save_GMG("x", Data_set; dir = joinpath(dir, "http_data"))
@test load_GMG("x", joinpath(dir, "http_data")).depth.val[20] == -15.0

# URL download is loaded from a temporary file, which is removed again, and rejects HTML served as data
tmp = mktempdir()
withenv("TMPDIR" => tmp, "TMP" => tmp, "TEMP" => tmp) do
    @test load_GMG("file://" * joinpath(dir, "test.jld2")).depth.val[20] == -15.0
end
@test isempty(readdir(tmp))
# a binary download (a PNG header, not valid UTF-8) is accepted
png = joinpath(dir, "image.png")
write(png, UInt8[0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a, 0xff, 0xfe, 0x00, 0x00, 0x00, 0x0d, 0x49, 0x48])
@test read(download_data("file://" * png, "copy.png"; dir, maxattempts = 1)) == read(png)
write(joinpath(dir, "bad.jld2"), "<html>not found</html>")
@test_throws "Could not download" load_GMG("file://" * joinpath(dir, "bad.jld2"); maxattempts = 1)

if remote_ok
    data_remote = load_GMG(url)
    @test data_remote.fields.MohoDepth[20] ≈ -17.99km

    data_remote = download_data(url, "temp1.dat"; dir)
    @test data_remote[(end - 8):end] == "temp1.dat"
else
    @test_skip "file server unreachable: remote download tests skipped"
end
