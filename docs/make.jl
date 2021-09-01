using TMI
using Documenter

DocMeta.setdocmeta!(TMI, :DocTestSetup, :(using TMI); recursive=true)

makedocs(;
    modules=[TMI],
    authors="G Jake Gebbie",
    repo="https://github.com/ggebbie/TMI.jl/blob/{commit}{path}#{line}",
    sitename="TMI.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ggebbie.github.io/TMI.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ggebbie/TMI.jl",
)
