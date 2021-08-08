"""
To use this abstract type and avaliable functions,
create a subtype and implement the necessary functions.

Example:
    mutable struct Regression <: LD_Score_Regression
        ...

        function Regression(...)
            new(...)
        end
    end

    function update_weights(reg::Regression, ...)
    end

See `Hsq.jl` for a case study example.
"""

abstract type LD_Score_Regression end

# TODO need to find jackknife library or implement it
#function _delete_vals_tot(reg::LD_Score_Regression, jknife, Nbar, M) end
#function _delete_vals_part(reg::LD_Score_Regression, jknife, Nbar, M) end
#function _coef(reg::LD_Score_Regression, jknife, Nbar) end
#function _cat(reg::LD_Score_Regression, jknife, M, Nbar, coef, coef_cov) end
#function _tot(reg::LD_Score_Regression, cat, cat_cov) end
#function _prop(reg::LD_Score_Regression, jknife, M, Nbar, cat, tot) end
#function _enrichment(reg::LD_Score_Regression, M, M_tot, cat, tot) end
#function _intercept(reg::LD_Score_Regression, jknife) end
#function _combine_twostep_jknives(
#    reg::LD_Score_Regression, step1_jknife, step2_jknife, M_tot, c; Nbar=1,
#) end
#function _delete_vals_tot(reg::LD_Score_Regression) end


function aggregate(reg::LD_Score_Regression, y, x, N, M; intercept=nothing)
    if intercept == nothing
        intercept = reg.__null_intercept__
    end

    num = M * (mean(y) - intercept)
    denom = mean(x .* N)
    return num / denom
end


function ld_score_regression(
    reg::LD_Score_Regression,
    y, x, w, N, M, n_blocks, intercept, slow, step1_ii, old_weights,
)
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
        yp = y - intercept
        intercept_se = "NA"
    end

    reg.twostep_filtered = nothing
    if step1_ii != nothing && reg.constrain_intercept
        throw(ErrorException(
            "twostep is not compatible with constrain_intercept."
        ))
    elseif !(step1_ii == nothing) && reg.n_annot > 1
        throw(ErrorException(
            "twostep not compatible with partitioned LD Score yet."
        ))
    elseif !(step1_ii == nothing)
        n1 = sum(step1_ii)
        reg.twostep_filtered = n_snp - n1
        s = dropdims(step1_ii; dims=2)
        x1 = zeros(0)
        for (i, j) in enumerate(s)
            if s[i] == 1
                append!(x1, x[i])
            end
        end

    # TODO finish implementing this function
    end
end
