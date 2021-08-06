abstract type LD_Score_Regression end

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

    constrain_intercept
    n_annot
    __null_intercept__
    function Hsq(
        y, x, w, N, M, n_blocks, intercept, slow, step1_ii, old_weights,
    )
        hsq = new(y, x, w, N, M, n_blocks, intercept, slow, step1_ii, old_weights)
        hsq.__null_intercept__ = 1
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

    print(add(ld_score_regression, 10))

    ld_score_regression
end

function aggregate(
    ld_score_regression::LD_Score_Regression, y, x_tot, N, M_tot, intercept,
)
    if intercept == nothing intecept = ld_score_regression.__null_intercept__ end
end
