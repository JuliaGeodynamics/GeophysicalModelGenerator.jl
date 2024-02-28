using GeophysicalModelGenerator
using .GeophysicalModelGenerator
using Documenter

#DocMeta.setdocmeta!(GeophysicalModelGenerator, :DocTestSetup, :(using GeophysicalModelGenerator); recursive=true)

# Get GeophysicalModelGenerator.jl root directory
GMG_root_dir = dirname(@__DIR__)

# Copy list of authors to not need to synchronize it manually
authors_text = read(joinpath(GMG_root_dir, "AUTHORS.md"), String)
authors_text = replace(authors_text, "in the [LICENSE.md](LICENSE.md) file" => "under [License](@ref)")
write(joinpath(@__DIR__, "src", "man", "authors.md"), authors_text)
#Contributing
contributing = read(joinpath(GMG_root_dir, "CONTRIBUTING.md"), String)
write(joinpath(@__DIR__, "src", "man", "contributing.md"), contributing)

# Copy some files from the repository root directory to the docs and modify them
# as necessary
# Based on: https://github.com/ranocha/SummationByPartsOperators.jl/blob/0206a74140d5c6eb9921ca5021cb7bf2da1a306d/docs/make.jl#L27-L41
open(joinpath(@__DIR__, "src", "man", "code_of_conduct.md"), "w") do io
  # Point to source license file
  println(io, """
  ```@meta
  EditURL = "https://github.com/JuliaGeodynamics/GeophysicalModelGenerator.jl/blob/main/CODE_OF_CONDUCT.md"
  ```
  """)
  # Write the modified contents
  println(io, "# [Code of Conduct](@id code-of-conduct)")
  println(io, "")
  for line in eachline(joinpath(dirname(@__DIR__), "CODE_OF_CONDUCT.md"))
    line = replace(line, "[AUTHORS.md](AUTHORS.md)" => "[Authors](@ref)")
    println(io, "> ", line)
  end
end
open(joinpath(@__DIR__, "src", "man", "contributing.md"), "w") do io
    # Point to source license file
    println(io, """
    ```@meta
    EditURL = "https://github.com/JuliaGeodynamics/GeophysicalModelGenerator.jl/blob/main/CONTRIBUTING.md"
    ```
    """)
    # Write the modified contents
    println(io, "# [Contributing](@id contributing)")
    println(io, "")
    for line in eachline(joinpath(dirname(@__DIR__), "CONTRIBUTING.md"))
      line = replace(line, "[LICENSE.md](LICENSE.md)" => "[License](@ref)")
      line = replace(line, "[AUTHORS.md](AUTHORS.md)" => "[Authors](@ref)")
      println(io, "> ", line)
    end
  end

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
            "Read UTM data" =>  "man/tutorial_UTM.md",
            "VoteMaps" =>  "man/Tutorial_Votemaps.md",
            "Kilometer-scale volcano" =>  "man/tutorial_local_Flegrei.md",
            "Generating LaMEM model" =>  "man/LaPalma_example.md",
            "Create movies" =>  "man/tutorial_time_Seismicity.md"
        ],
        "User Guide" => Any[
            "Installation" =>  "man/installation.md",
            "Data Structures" =>  "man/datastructures.md",
            "Data Import" =>  "man/dataimport.md",
            "Projection" =>  "man/projection.md",
            "Paraview output" => "man/paraview_output.md",
            "Tools" => "man/tools.md",
            "Visualisation" => "man/visualise.md",
            "Gravity code" => "man/gravity_code.md",
            "LaMEM" => "man/lamem.md",
            "Profile Processing" => "man/profile_processing.md"
        ],
        "List of functions"  => "man/listfunctions.md",
        "Authors" => "man/authors.md",
        "Contributing" => "man/contributing.md",
        "Code of Conduct" => "man/code_of_conduct.md",
        "License" => "man/license.md",
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
