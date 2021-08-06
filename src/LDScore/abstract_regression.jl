include("LD_Score_Regression.jl")

function aggregate(
    reg::LD_Score_Regression,
    y, x_tot, N, M_tot, intercept,
)
    if intercept == nothing intecept = reg.__null_intercept__ end
end

function make_ld_score_regression(
    reg::LD_Score_Regression,
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

    n_snp, reg.n_annot = size(x)
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
    reg.constrain_intercept = intercept != nothing
    reg.intercept = intercept
    reg.n_blocks = n_blocks

    tot_agg = aggregate(reg, y, x_tot, N, M_tot, intercept)
    initial_w = _update_weights(
        reg, x_tot, w, N, M_tot, tot_agg, intercept,
    )
    Nbar = mean(N)
    x = N .* x / Nbar

    if !reg.constrain_intercept
        x, x_tot = append_intercept(x), append_intercept(x_tot)
        yp = y
    else
        yp = y - intercept
        intercept_se = "NA"
    end

end

function append_intercept(x)
    print(x)
    #=
    n_row = x.shape[0]
    intercept = np.ones((n_row, 1))
    x_new = np.concatenate((x, intercept), axis=1)
    return x_new
    =#
end
