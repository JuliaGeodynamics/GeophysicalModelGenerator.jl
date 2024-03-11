using Test

dir = "test_files"
outname1 = movie_from_images(dir=dir)
@test outname1=="test_animation.mp4"

outname2 = movie_from_images(dir=dir, copy_to_current_dir=false)
@test outname2=="test_animation.mp4"

outname3 = movie_from_images(dir=dir, copy_to_current_dir=false,  framerate=20, outfile="test_anim2")
@test outname3=="test_anim2.mp4"

outname4 = movie_from_images(dir=dir, copy_to_current_dir=false,  framerate=20, outfile="test_anim3", type=:mov_hires)
@test outname4=="test_anim3.mov"

rm(outname1)
rm(joinpath(dir,outname2))
rm(joinpath(dir,outname3))
rm(joinpath(dir,outname4))


