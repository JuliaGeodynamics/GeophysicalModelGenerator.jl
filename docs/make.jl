using GeophysicalModelGenerator
using Documenter

#DocMeta.setdocmeta!(GeophysicalModelGenerator, :DocTestSetup, :(using GeophysicalModelGenerator); recursive=true)

makedocs(;
    modules=[GeophysicalModelGenerator],
    authors="Marcel Thielmann, Boris Kaus",
    sitename="GeophysicalModelGenerator.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
    ),
    pages=[
        "Home" => "index.md",
        "Tutorials" => Any[
            "Overview" =>  "man/tutorials.md",
            "3D seismic tomography from ASCII" =>  "man/tutorial_load3DSeismicData.md",
            "3D seismic tomography from netCDF" =>  "man/tutorial_loadregular3DSeismicData_netCDF.md",
            "Visualize Moho topography" =>  "man/tutorial_MohoTopo.md",
            "Create GMT-based topography" =>  "man/tutorial_GMT_Topography.md",
            "Coastlines" =>  "man/tutorial_Coastlines.md",
            "Import screenshots" =>  "man/tutorial_Screenshot_To_Paraview.md",
            "Interpolate irregular 3D seismic tomography" =>  "man/tutorial_loadirregular3DSeismicData.md",
            "ETOPO1 Topography and geological maps" =>  "man/tutorial_GMT_Topography_GeologicalMap.md",
            "ISC earthquake data" =>  "man/tutorial_ISC_data.md",
            "Plot GPS vectors" =>  "man/tutorial_GPS.md",
            "Read UTM data" =>  "man/tutorial_UTM.md"    
        ],
        "User Guide" => Any[
            "Installation" =>  "man/installation.md",
            "Data Structures" =>  "man/datastructures.md",
            "Data Import" =>  "man/dataimport.md",
            "Paraview output" => "man/paraview_output.md",
            "Tools" => "man/tools.md",
            "Gravity code" => "man/gravity_code.md",
            "LaMEM" => "man/lamem.md"
        ],
        "List of functions"  => "man/listfunctions.md"
    ],
)

deploydocs(;
    repo="github.com/JuliaGeodynamics/GeophysicalModelGenerator.jl.git",
    branch = "gh-pages",
    target = "build",
    devbranch = "main",
    devurl = "dev",
    forcepush=true,
    push_preview = true
)

