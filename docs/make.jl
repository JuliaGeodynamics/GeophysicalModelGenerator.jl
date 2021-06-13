using GeophysicalModelGenerator
using Documenter

DocMeta.setdocmeta!(GeophysicalModelGenerator, :DocTestSetup, :(using GeophysicalModelGenerator); recursive=true)

makedocs(;
    modules=[GeophysicalModelGenerator],
    authors="Marcel Thielmann, Boris Kaus",
    repo="https://github.com/JuliaGeodynamics/GeophysicalModelGenerator/{commit}{path}#{line}",
    sitename="GeophysicalModelGenerator.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://github.com/JuliaGeodynamics/GeophysicalModelGenerator.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Tutorial" => Any[
            "Load 3D seismic tomography from CSV" =>  "man/tutorial_load3DSeismicData.md",
        ],
        "User Guide" => Any[
            "Data Structures" =>  "man/datastructures.md",
            "Paraview output" => "man/paraview_output.md"
        ],
        "List of functions"  => "man/listfunctions.md"
    ],
)

deploydocs(;
    repo="github.com/JuliaGeodynamics/GeophysicalModelGenerator.jl",
    devbranch = "main"
)
