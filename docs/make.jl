using GeophysicalModelGenerator
using Documenter

DocMeta.setdocmeta!(GeophysicalModelGenerator, :DocTestSetup, :(using GeophysicalModelGenerator); recursive=true)

makedocs(;
    modules=[GeophysicalModelGenerator],
    authors="Marcel Thielmann, Boris Kaus",
    repo="https://github.com/mthielma/GeophysicalModelGenerator.jl/blob/{commit}{path}#{line}",
    sitename="GeophysicalModelGenerator.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mthielma.github.io/GeophysicalModelGenerator.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mthielma/GeophysicalModelGenerator.jl",
)
