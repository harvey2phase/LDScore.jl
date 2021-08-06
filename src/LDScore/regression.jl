using Statistics

include("abstract_regression.jl")

mutable struct Hsq <: LD_Score_Regression
    y
    x
    w
    N
    M
    n_blocks
    intercept
    slow
    step1_ii
    old_weights
    __null_intercept__

    constrain_intercept
    n_annot
    intercept_se
    twostep_filtered
    function Hsq(
        y, x, w, N, M, n_blocks, intercept, slow, step1_ii, old_weights,
    )
        hsq = new(
            y, x, w, N, M, n_blocks, intercept, slow, step1_ii, old_weights, 1,
        )
        return hsq
    end
end
