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

    """
    Use these default values:
        - n_blocks = 200
        - intercept = nothing
        - slow = false
        - twostep = nothing
        - old_weights = false
    """
    function Hsq(
        y, x, w, N, M; n_blocks = 200, slow = false, old_weights = false,
        intercept = nothing, step1_ii = nothing,
    )
        hsq = new(
            y, x, w, N, M, n_blocks, intercept, slow, step1_ii, old_weights, 1,
        )
        return ld_score_regression(
            hsq,
            y, x, w, N, M, n_blocks, intercept, slow, step1_ii, old_weights,
        )
    end
end


function hsq_weights(ld, w_ld, N, M, hsq; intercept=nothing)
    if intercept == nothing intercept = 1.0 end
    M = M[1][1]
    hsq = max(maximum(hsq), 0.0)
    hsq = min(hsq, 1.0)
    ld = fmax(ld, 1.0)
    w_ld = fmax(w_ld, 1.0)
    c = hsq .* N ./ M

    het_w = 1.0 ./ (2 .* (c .* ld .+ intercept).^2)
    oc_w = 1.0 ./ w_ld
    w = het_w .* oc_w
    return w
end


function update_weights(reg::Hsq, ld, w_ld, N, M, hsq, intercept)
    if intercept == nothing intercept = reg.__null_intercept__ end

    return hsq_weights(ld, w_ld, N, M, hsq; intercept=intercept)
end

function overlap_output(
    hsq::Hsq,
    category_names, overlap_matrix, M_annot, M_tot, print_coefficients,
)
    overlap_matrix_prop = zeros([hsq.n_annot, hsq.n_annot])
    for i in range(hsq.n_annot)
        overlap_matrix_prop[i, :] = overlap_matrix[i, :] / M_annot
    end

    # TODO missing code here

    #=
    df = DataFrames.DataFrame({
        "Category": category_names,
        "Prop._SNPs": one_d_convert(prop_M_overlap),
        "Prop._h2": one_d_convert(prop_hsq_overlap),
        "Prop._h2_std_error": one_d_convert(prop_hsq_overlap_se),
        "Enrichment": one_d_convert(enrichment),
        "Enrichment_std_error": one_d_convert(enrichment_se),
        "Enrichment_p":diff_p,
        "Coefficient": one_d_convert(hsq.coef),
        "Coefficient_std_error": hsq.coef_se,
        "Coefficient_z-score": one_d_convert(hsq.coef) / one_d_convert(hsq.coef_se)
    })
    if print_coefficients
        return df[
            !,
            [
                "Category", "Prop._SNPs", "Prop._h2", "Prop._h2_std_error",
                "Enrichment","Enrichment_std_error", "Enrichment_p",
                "Coefficient", "Coefficient_std_error","Coefficient_z-score",
            ],
        ]
    else
        return df[
            !,
            [
                "Category", "Prop._SNPs", "Prop._h2", "Prop._h2_std_error",
                "Enrichment","Enrichment_std_error", "Enrichment_p",
            ],
        ]
    end
    =#
end
