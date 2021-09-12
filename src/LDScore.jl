module LDScore
    include("ldsc.jl")
    include("LDScore/LDScoreRegression.jl")
    include("LDScore/Hsq.jl")
    include("LDScore/regression.jl")
    include("LDScore/sumstats.jl")
    include("LDScore/parse.jl")

    export hsq_weights
    export aggregate
    export estimate_h2
    export _read_sumstats
    export parse_sumstats
    export args
    export _update_func
end
