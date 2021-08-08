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
        y, x, w, N, M, n_blocks, intercept, slow, step1_ii, old_weights,
    )
        hsq = new(
            y, x, w, N, M, n_blocks, intercept, slow, step1_ii, old_weights, 1,
        )
        return hsq
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

    prop_hsq_overlap = np.dot(
        overlap_matrix_prop, hsq.prop.T).reshape((1, hsq.n_annot))
    prop_hsq_overlap_var = np.diag(
        np.dot(np.dot(overlap_matrix_prop, hsq.prop_cov), overlap_matrix_prop.T))
    prop_hsq_overlap_se = np.sqrt(
        np.maximum(0, prop_hsq_overlap_var)).reshape((1, hsq.n_annot))
    one_d_convert = lambda x: np.array(x).reshape(np.prod(x.shape))
    prop_M_overlap = M_annot / M_tot
    enrichment = prop_hsq_overlap / prop_M_overlap
    enrichment_se = prop_hsq_overlap_se / prop_M_overlap
    overlap_matrix_diff = np.zeros([hsq.n_annot,hsq.n_annot])
    for i in range(hsq.n_annot):
        if not M_tot == M_annot[0,i]:
            overlap_matrix_diff[i, :] = overlap_matrix[i,:]/M_annot[0,i] - \
                (M_annot - overlap_matrix[i,:]) / (M_tot-M_annot[0,i])

    diff_est = np.dot(overlap_matrix_diff,hsq.coef)
    diff_cov = np.dot(np.dot(overlap_matrix_diff,hsq.coef_cov),overlap_matrix_diff.T)
    diff_se = np.sqrt(np.diag(diff_cov))
    diff_p = ["NA" if diff_se[i]==0 else 2*tdist.sf(abs(diff_est[i]/diff_se[i]),hsq.n_blocks) \
        for i in range(hsq.n_annot)]

    df = pd.DataFrame({
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
    if print_coefficients:
        df = df[["Category", "Prop._SNPs", "Prop._h2", "Prop._h2_std_error",
                "Enrichment","Enrichment_std_error", "Enrichment_p",
                 "Coefficient", "Coefficient_std_error","Coefficient_z-score"]]
    else:
        df = df[["Category", "Prop._SNPs", "Prop._h2", "Prop._h2_std_error",
                "Enrichment","Enrichment_std_error", "Enrichment_p"]]
    return df
