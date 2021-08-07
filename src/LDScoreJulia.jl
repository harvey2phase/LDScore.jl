module LDScoreJulia
    include("LDScore/regression.jl")
    include("LDScore/sumstats.jl")

    export hsq_weights
    export aggregate
    export estimate_h2
end
