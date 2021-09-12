using LDScore
using Documenter

DocMeta.setdocmeta!(LDScore, :DocTestSetup, :(using LDScore); recursive=true)

makedocs(;
    modules=[LDScore],
    authors="Harvey Wang",
    repo="https://github.com/harvey2phase/LDScore.jl/blob/{commit}{path}#{line}",
    sitename="LDScore.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://harvey2phase.github.io/LDScore.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/harvey2phase/LDScore.jl",
)
