# This takes various figures from the following paper
#
# Lippitsch, R., et al. 2003. Upper mantle structure beneath the Alpine orogen from high-resolution teleseismic tomography. J. Geophys. Res. 108, 2376. https://doi.org/10.1029/2002JB002016
#
# and transfers it into paraview format
#
# For convenience we created a zip file with all cross-sections and mapviews here:
#
# 

using GeophysicalModelGenerator

# Process cross-sections of Figure 13. Note that we estimated some of the lon/lat locations
data_Fig13a         =   Screenshot_To_GeoData("Lippitsch_Fig13a.png",( 3.4, 46.0, -400.0), (16.9, 42.5, 0.0))
Write_Paraview(data_Fig13a, "Lippitsch_Fig13a") 

data_Fig13b         =   Screenshot_To_GeoData("Lippitsch_Fig13b.png",( 5.7, 52.0, -400.0), (11.6, 43.5, 0.0))
Write_Paraview(data_Fig13b, "Lippitsch_Fig13b") 

data_Fig13c         =   Screenshot_To_GeoData("Lippitsch_Fig13c.png",(17.7, 51.0, -400.0), (10.7, 43.5, 0.0))
Write_Paraview(data_Fig13c, "Lippitsch_Fig13c") 

# Mapview images
Corner_LowerLeft    =   ( 3.5, 43.0 , -150.0)
Corner_UpperRight   =   (15.5, 50.0 , -150.0)
Corner_LowerRight   =   (15.5, 43.0 , -150.0)
Corner_UpperLeft    =   (3.5 , 50.0 , -150.0)
data_Fig13_map      =   Screenshot_To_GeoData("Fig13_mapview.png",Corner_LowerLeft, Corner_UpperRight, Corner_LowerRight=Corner_LowerRight,Corner_UpperLeft=Corner_UpperLeft)
Write_Paraview(data_Fig13_map, "Lippitsch_Fig13_mapview") 

Depth                   =   -90.0;
Corner_LowerLeft        =   (Corner_LowerLeft[1],   Corner_LowerLeft[2],    Depth)
Corner_UpperRight       =   (Corner_UpperRight[1],  Corner_UpperRight[2],   Depth)
Corner_LowerRight       =   (Corner_LowerRight[1],  Corner_LowerRight[2],   Depth)
Corner_UpperLeft        =   (Corner_UpperLeft[1],   Corner_UpperLeft[2],    Depth)
data_Fig12_90km    =   Screenshot_To_GeoData("Fig12_90km.png",Corner_LowerLeft, Corner_UpperRight, Corner_LowerRight=Corner_LowerRight,Corner_UpperLeft=Corner_UpperLeft)
Write_Paraview(data_Fig12_90km, "Lippitsch_Fig12_90km") 

Depth                   = -180.0;
Corner_LowerLeft        =   (Corner_LowerLeft[1],   Corner_LowerLeft[2],    Depth)
Corner_UpperRight       =   (Corner_UpperRight[1],  Corner_UpperRight[2],   Depth)
Corner_LowerRight       =   (Corner_LowerRight[1],  Corner_LowerRight[2],   Depth)
Corner_UpperLeft        =   (Corner_UpperLeft[1],   Corner_UpperLeft[2],    Depth)
data_Fig12_180km   =   Screenshot_To_GeoData("Fig12_180km.png",Corner_LowerLeft, Corner_UpperRight, Corner_LowerRight=Corner_LowerRight,Corner_UpperLeft=Corner_UpperLeft)
Write_Paraview(data_Fig12_180km, "Lippitsch_Fig12_180km") 

Depth                   = -300.0;
Corner_LowerLeft        =   (Corner_LowerLeft[1],   Corner_LowerLeft[2],    Depth)
Corner_UpperRight       =   (Corner_UpperRight[1],  Corner_UpperRight[2],   Depth)
Corner_LowerRight       =   (Corner_LowerRight[1],  Corner_LowerRight[2],   Depth)
Corner_UpperLeft        =   (Corner_UpperLeft[1],   Corner_UpperLeft[2],    Depth)
data_Fig12_300km   =   Screenshot_To_GeoData("Fig12_300km.png",Corner_LowerLeft, Corner_UpperRight, Corner_LowerRight=Corner_LowerRight,Corner_UpperLeft=Corner_UpperLeft)
Write_Paraview(data_Fig12_300km, "Lippitsch_Fig12_300km") 

Depth                   = -400.0;
Corner_LowerLeft        =   (Corner_LowerLeft[1],   Corner_LowerLeft[2],    Depth)
Corner_UpperRight       =   (Corner_UpperRight[1],  Corner_UpperRight[2],   Depth)
Corner_LowerRight       =   (Corner_LowerRight[1],  Corner_LowerRight[2],   Depth)
Corner_UpperLeft        =   (Corner_UpperLeft[1],   Corner_UpperLeft[2],    Depth)
data_Fig12_400km   =   Screenshot_To_GeoData("Fig12_400km.png",Corner_LowerLeft, Corner_UpperRight, Corner_LowerRight=Corner_LowerRight,Corner_UpperLeft=Corner_UpperLeft)
Write_Paraview(data_Fig12_400km, "Lippitsch_Fig12_400km") 
