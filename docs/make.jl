using Documenter

push!(LOAD_PATH, dirname(@__DIR__))

using GeophysicalModelGenerator

# Importing these activates package extensions
import GLMakie, GMT

#DocMeta.setdocmeta!(GeophysicalModelGenerator, :DocTestSetup, :(using GeophysicalModelGenerator); recursive=true)

# Get GeophysicalModelGenerator.jl root directory
GMG_root_dir = dirname(@__DIR__)

license = read(joinpath(GMG_root_dir, "LICENSE.md"), String)
write(joinpath(@__DIR__, "src", "man", "license.md"), license)
# Copy list of authors to not need to synchronize it manually
authors_text = read(joinpath(GMG_root_dir, "AUTHORS.md"), String)
# authors_text = replace(authors_text, "in the [LICENSE.md](LICENSE.md) file" => "under [License](@ref)")
write(joinpath(@__DIR__, "src", "man", "authors.md"), authors_text)

# Copy some files from the repository root directory to the docs and modify them
# as necessary
# Based on: https://github.com/ranocha/SummationByPartsOperators.jl/blob/0206a74140d5c6eb9921ca5021cb7bf2da1a306d/docs/make.jl#L27-L41
open(joinpath(@__DIR__, "src", "man", "license.md"), "w") do io
  # Point to source license file
  println(io, """
  ```@meta
  EditURL = "https://github.com/JuliaGeodynamics/GeophysicalModelGenerator.jl/blob/main/LICENSE.md"
  ```
  """)
  # Write the modified contents
  println(io, "# [License](@id license)")
  println(io, "")
  for line in eachline(joinpath(dirname(@__DIR__), "LICENSE.md"))
    line = replace(line, "[AUTHORS.md](AUTHORS.md)" => "[Authors](@ref)")
    println(io, "> ", line)
  end
end

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
    authors="Boris Kaus, Marcel Thielmann",
    sitename="GeophysicalModelGenerator.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
    ),
    pages=[
        "Home" => "index.md",
        "Tutorials" => Any[
            "Overview" =>  "man/tutorials.md",
            "1 - 3D seismic tomography from ASCII" =>  "man/tutorial_load3DSeismicData.md",
            "2 - 3D seismic tomography from netCDF" =>  "man/tutorial_loadregular3DSeismicData_netCDF.md",
            "3 - Visualize Moho topography" =>  "man/Tutorial_MohoTopo_Spada.md",
            "4 - Create GMT-based topography" =>  "man/tutorial_GMT_Topography.md",
            "5 - Coastlines" =>  "man/tutorial_Coastlines.md",
            "6 - Import screenshots" =>  "man/tutorial_Screenshot_To_Paraview.md",
            "7 - Interpolate irregular 3D seismic tomography" =>  "man/tutorial_loadirregular3DSeismicData.md",
            "8 - ETOPO1 Topography and geological maps" =>  "man/tutorial_GMT_Topography_GeologicalMap.md",
            "9 - ISC earthquake data" =>  "man/tutorial_ISC_data.md",
            "10 - Plot GPS vectors" =>  "man/tutorial_GPS.md",
            "11 - Read UTM data" =>  "man/tutorial_UTM.md",
            "12 - VoteMaps" =>  "man/Tutorial_Votemaps.md",
            "13 - Campi Flegrei" =>  "man/tutorial_local_Flegrei.md",
            "14 - La Palma volcano Model" =>  "man/Tutorial_LaPalma.md",
            "15 - Create movies" =>  "man/tutorial_time_Seismicity.md",
            "16 - Fault Density Map" =>  "man/Tutorial_FaultDensity.md",
            "17 - Alpine data integration" =>  "man/Tutorial_AlpineData.md",
            "18 - Jura tutorial" =>  "man/Tutorial_Jura.md"
        ],
        "User Guide" => Any[
            "Installation" =>  "man/installation.md",
            "Data Structures" =>  "man/datastructures.md",
            "Data Import" =>  "man/dataimport.md",
            "Projection" =>  "man/projection.md",
            "Surfaces" =>  "man/surfaces.md",
            "Paraview output" => "man/paraview_output.md",
            "Paraview collection" => "man/paraview_collection.md",
            "Tools" => "man/tools.md",
            "Visualisation" => "man/visualise.md",
            "Gravity code" => "man/gravity_code.md",
            "Geodynamic setups" => "man/geodynamic_setups.md",
            "LaMEM" => "man/lamem.md",
            "Profile Processing" => "man/profile_processing.md",
            "Movies" => "man/movies.md"
        ],
        "List of functions"  => "man/listfunctions.md",
        "Authors" => "man/authors.md",
        "Contributing" => "man/contributing.md",
        "Code of Conduct" => "man/code_of_conduct.md",
        "License" => "man/license.md"
    ],
    pagesonly=true,
    warnonly=true
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
