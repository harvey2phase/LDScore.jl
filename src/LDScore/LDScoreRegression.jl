include("dev_tools.jl")

"""
    LDScoreRegression

An abstract type for basic LDSC regression functionalities.
Subtypes and necessary functions should be implemented as needed;
existing subtypes include `Hsq`.

# Examples
```
mutable struct Regression <: LDScoreRegression
    ...

    function Regression(...)
        new(...)
    end
end

function update_weights(reg::Regression, ...) end
```
See `Hsq.jl` for a case study example.
"""
abstract type LDScoreRegression end


function aggregate(reg::LDScoreRegression, y, x, N, M; intercept = nothing)
    if intercept == nothing
        intercept = reg.__null_intercept__
    end

    num = M * (mean(y) - intercept)
    denom = mean(x .* N)
    return num / denom
end # function aggregate


function ld_score_regression(
    reg::LDScoreRegression,
    y, x, w, N, M, n_blocks;
    intercept = nothing, slow = false, step1_ii = nothing, old_weights = false,
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

    M_tot = sum(M)
    # shape should be [n_snp, 1]
    x_tot = sum(x, dims = 2)
    reg.constrain_intercept = intercept != nothing
    reg.intercept = intercept
    reg.n_blocks = n_blocks

    tot_agg = aggregate(reg, y, x_tot, N, M_tot; intercept=intercept)
    initial_w = update_weights(
        reg, x_tot, w, N, M_tot, tot_agg, intercept,
    )
    Nbar = mean(N)
    x = N .* x / Nbar

    if !reg.constrain_intercept
        x, x_tot = append_intercept(x), append_intercept(x_tot)
        yp = y
    else
        yp = y .- intercept
        intercept_se = "NA"
    end

    reg.twostep_filtered = nothing
    if step1_ii != nothing && reg.constrain_intercept
        throw(ErrorException(
            "twostep is not compatible with constrain_intercept."
        ))
    elseif step1_ii != nothing && reg.n_annot > 1
        throw(ErrorException(
            "twostep not compatible with partitioned LD Score yet."
        ))
    elseif step1_ii != nothing
        n1 = sum(step1_ii)
        reg.twostep_filtered = n_snp - n1
        s = dropdims(step1_ii; dims=2)
        x1 = zeros(0)
        for (i, j) in enumerate(s)
            if s[i] == 1
                append!(x1, x[i])
            end
        end
        yp1, w1, N1, initial_w1 = map(
            (a) -> reshape(a[step1_ii], (n1, 1)),
            (yp, w, N, initial_w),
        )
        update_func1 = reg._update_func(a, x1, w1, N1, M_tot, Nbar, ii=step1_ii)
        #= TODO
        step1_jknife = IRWLS(
            x1, yp1, update_func1, n_blocks, slow=slow, w=initial_w1)
        step1_int, _ = self._intercept(step1_jknife)
        yp = yp - step1_int
        x = remove_intercept(x)
        x_tot = remove_intercept(x_tot)
        update_func2 = lambda a: self._update_func(
            a, x_tot, w, N, M_tot, Nbar, step1_int)
        s = update_separators(step1_jknife.separators, step1_ii)
        step2_jknife = IRWLS(
            x, yp, update_func2, n_blocks, slow=slow, w=initial_w, separators=s)
        c = np.sum(np.multiply(initial_w, x)) / \
            np.sum(np.multiply(initial_w, np.square(x)))
        jknife = self._combine_twostep_jknives(
            step1_jknife, step2_jknife, M_tot, c, Nbar)
        =#
    elseif old_weights
        #= TODO
        initial_w = np.sqrt(initial_w)
        x = IRWLS._weight(x, initial_w)
        y = IRWLS._weight(yp, initial_w)
        jknife = jk.LstsqJackknifeFast(x, y, n_blocks)
        =#
    else
        #= TODO
        update_func = lambda a: self._update_func(
            a, x_tot, w, N, M_tot, Nbar, intercept)
        jknife = IRWLS(
            x, yp, update_func, n_blocks, slow=slow, w=initial_w)
        =#
    end

    # TODO finish implementing this function


    return reg
end # function ld_score_regression


# TODO need to find jackknife library or implement it
#function _delete_vals_tot(reg::LDScoreRegression, jknife, Nbar, M) end
#function _delete_vals_part(reg::LDScoreRegression, jknife, Nbar, M) end
#function _coef(reg::LDScoreRegression, jknife, Nbar) end
#function _cat(reg::LDScoreRegression, jknife, M, Nbar, coef, coef_cov) end
#function _tot(reg::LDScoreRegression, cat, cat_cov) end
#function _prop(reg::LDScoreRegression, jknife, M, Nbar, cat, tot) end
#function _enrichment(reg::LDScoreRegression, M, M_tot, cat, tot) end
#function _intercept(reg::LDScoreRegression, jknife) end
#function _combine_twostep_jknives(
#    reg::LDScoreRegression, step1_jknife, step2_jknife, M_tot, c; Nbar=1,
#) end
#function _delete_vals_tot(reg::LDScoreRegression) end
