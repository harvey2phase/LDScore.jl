abstract type LD_Score_Regression end

function aggregate(
    ld_score_regression,
    y, x_tot, N, M_tot, intercept,
)
    if intercept == nothing intecept = ld_score_regression.__null_intercept__ end
end

function _update_weights(
    ld_score_regression,
    x_tot, w, N, M_tot, tot_agg, intercept,
)
end
