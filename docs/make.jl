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
        "Top-level functions" => "top.md",
        "Home" => "index.md",
        "Configuration" => "config.md",
        "Grid" => "grid.md",
        "Boundary Conditions" => "boundaries.md",
        "Plots" => "plots.md",
        "Utilities" => "utils.md",
        "Older functions" => "legacy.md",
    ],
)

deploydocs(;
           repo="github.com/ggebbie/TMI.jl",
           devbranch="main",
)
