using Test

dir = "test_files"
outname1 = movie_from_images(dir = dir)
@test outname1 == "test_animation.mp4"

outname2 = movie_from_images(dir = dir, copy_to_current_dir = false)
@test outname2 == "test_animation.mp4"

outname3 = movie_from_images(dir = dir, copy_to_current_dir = false, framerate = 20, outfile = "test_anim2")
@test outname3 == "test_anim2.mp4"

outname4 = movie_from_images(dir = dir, copy_to_current_dir = false, framerate = 20, outfile = "test_anim3", type = :mov_hires)
@test outname4 == "test_anim3.mov"

rm(outname1)
rm(joinpath(dir, outname2))
rm(joinpath(dir, outname3))
rm(joinpath(dir, outname4))

# default dir=pwd() with copy_to_current_dir=true must not copy the movie onto itself (#107)
mktempdir() do tmp
    for f in readdir(dir)
        endswith(f, ".png") && cp(joinpath(dir, f), joinpath(tmp, f))
    end
    cd(tmp) do
        @test movie_from_images() == "test_animation.mp4"
        @test isfile("test_animation.mp4")
    end
end
@test_throws "no images matching" movie_from_images(dir = dir, file = "nonexistent")
