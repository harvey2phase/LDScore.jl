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
    function Hsq(
        y, x, w, N, M, n_blocks, intercept, slow, step1_ii, old_weights,
    )
        hsq = new(
            y, x, w, N, M, n_blocks, intercept, slow, step1_ii, old_weights, 1,
        )
        return hsq
    end
end


function make_ld_score_regression(
    ld_score_regression::LD_Score_Regression,
    y, x, w, N, M, n_blocks, intercept, slow, step1_ii, old_weights,
)
    #= TODO
    for i in [y, x, w, M, N]:
        try:
            if len(i.shape) != 2:
                raise TypeError('Arguments must be 2D arrays.')
        except AttributeError:
            raise TypeError('Arguments must be arrays.')
    =#


    n_snp, ld_score_regression.n_annot = size(x)
    #= TODO
    if any(i.shape != (n_snp, 1) for i in [y, w, N]):
        raise ValueError(
            'N, weights and response (z1z2 or chisq) must have shape (n_snp, 1).')
    if M.shape != (1, self.n_annot):
        raise ValueError('M must have shape (1, n_annot).')
    =#

    M_tot = sum(M)
    # shape should be [n_snp, 1]
    x_tot = sum(x, dims = 2)'
    constrain_intercept = intercept != nothing

    tot_agg = aggregate(ld_score_regression, y, x_tot, N, M_tot, intercept)
    initial_w = _update_weights(
        ld_score_regression,
        x_tot, w, N, M_tot, tot_agg, intercept,
    )
    Nbar = mean(N)

end
