using LDScoreJulia
using Documenter

DocMeta.setdocmeta!(LDScoreJulia, :DocTestSetup, :(using LDScoreJulia); recursive=true)

makedocs(;
    modules=[LDScoreJulia],
    authors="Harvey Wang",
    repo="https://github.com/harvey2phase/LDScoreJulia.jl/blob/{commit}{path}#{line}",
    sitename="LDScoreJulia.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://harvey2phase.github.io/LDScoreJulia.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/harvey2phase/LDScoreJulia.jl",
)
