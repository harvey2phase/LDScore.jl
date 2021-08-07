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

# TODO add a test (works though)
function hsq_weights(ld, w_ld, N, M, hsq, intercept)
    if intercept == nothing intercept = 1 end
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

function _update_weights(
    reg::Hsq,
    ld, w_ld, N, M, hsq, intercept,
)
    if intercept == nothing intercept = reg.__null_intercept__ end

    return hsq_weights(ld, w_ld, N, M, hsq, intercept)
end
