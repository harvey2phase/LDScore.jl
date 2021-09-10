module LDScoreJulia
    include("ldsc.jl")
    include("LDScore/regression.jl")
    include("LDScore/sumstats.jl")

    export hsq_weights
    export aggregate
    export estimate_h2
    export _read_sumstats
    export parse_sumstats
    export args
end
