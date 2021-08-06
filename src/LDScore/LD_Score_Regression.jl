abstract type LD_Score_Regression end

function _delete_vals_tot(
    reg::LD_Score_Regression,
    )
end

#=
function _delete_vals_tot(reg::LD_Score_Regression, jknife, Nbar, M):
    # Get delete values for total h2 or gencov
    n_annot = reg.n_annot
    tot_delete_vals = jknife.delete_values[
        :, 0:n_annot]  # shape (n_blocks, n_annot)
    # shape (n_blocks, 1)
    tot_delete_vals = np.dot(tot_delete_vals, M.T) / Nbar
    return tot_delete_vals
end

function _delete_vals_part(reg::LD_Score_Regression, jknife, Nbar, M):
    # Get delete values for partitioned h2 or gencov.
    n_annot = reg.n_annot
    return jknife.delete_values[:, 0:n_annot] / Nbar
end

function _coef(reg::LD_Score_Regression, jknife, Nbar):
    # Get coefficient estimates + cov from the jackknife.#
    n_annot = reg.n_annot
    coef = jknife.est[0, 0:n_annot] / Nbar
    coef_cov = jknife.jknife_cov[0:n_annot, 0:n_annot] / Nbar ** 2
    coef_se = np.sqrt(np.diag(coef_cov))
    return coef, coef_cov, coef_se
end

function _cat(reg::LD_Score_Regression, jknife, M, Nbar, coef, coef_cov):
    # Convert coefficients to per-category h2 or gencov.
    cat = np.multiply(M, coef)
    cat_cov = np.multiply(np.dot(M.T, M), coef_cov)
    cat_se = np.sqrt(np.diag(cat_cov))
    return cat, cat_cov, cat_se
end

function _tot(reg::LD_Score_Regression, cat, cat_cov):
    # Convert per-category h2 to total h2 or gencov.
    tot = np.sum(cat)
    tot_cov = np.sum(cat_cov)
    tot_se = np.sqrt(tot_cov)
    return tot, tot_cov, tot_se
end

function _prop(reg::LD_Score_Regression, jknife, M, Nbar, cat, tot):
    # Convert total h2 and per-category h2 to per-category proportion h2 or gencov.
    n_annot = reg.n_annot
    n_blocks = jknife.delete_values.shape[0]
    numer_delete_vals = np.multiply(
        M, jknife.delete_values[:, 0:n_annot]) / Nbar  # (n_blocks, n_annot)
    denom_delete_vals = np.sum(
        numer_delete_vals, axis=1).reshape((n_blocks, 1))
    denom_delete_vals = np.dot(denom_delete_vals, np.ones((1, n_annot)))
    prop = jk.RatioJackknife(
        cat / tot, numer_delete_vals, denom_delete_vals)
    return prop.est, prop.jknife_cov, prop.jknife_se
end

function _enrichment(reg::LD_Score_Regression, M, M_tot, cat, tot):
    # Compute proportion of SNPs per-category enrichment for h2 or gencov.
    M_prop = M / M_tot
    enrichment = np.divide(cat, M) / (tot / M_tot)
    return enrichment, M_prop
end

function _intercept(reg::LD_Score_Regression, jknife):
    # Extract intercept and intercept SE from block jackknife.
    n_annot = reg.n_annot
    intercept = jknife.est[0, n_annot]
    intercept_se = jknife.jknife_se[0, n_annot]
    return intercept, intercept_se
end

function _combine_twostep_jknives(reg::LD_Score_Regression, step1_jknife, step2_jknife, M_tot, c, Nbar=1):
    # Combine free intercept and constrained intercept jackknives for --two-step.
    n_blocks, n_annot = step1_jknife.delete_values.shape
    n_annot -= 1
    if n_annot > 2:
        raise ValueError(
            'twostep not yet implemented for partitioned LD Score.')

    step1_int, _ = reg._intercept(step1_jknife)
    est = np.hstack(
        (step2_jknife.est, np.array(step1_int).reshape((1, 1))))
    delete_values = np.zeros((n_blocks, n_annot + 1))
    delete_values[:, n_annot] = step1_jknife.delete_values[:, n_annot]
    delete_values[:, 0:n_annot] = step2_jknife.delete_values -\
        c * (step1_jknife.delete_values[:, n_annot] -
             step1_int).reshape((n_blocks, n_annot))  # check this
    pseudovalues = jk.Jackknife.delete_values_to_pseudovalues(
        delete_values, est)
    jknife_est, jknife_var, jknife_se, jknife_cov = jk.Jackknife.jknife(
        pseudovalues)
    jknife = namedtuple('jknife',
                        ['est', 'jknife_se', 'jknife_est', 'jknife_var', 'jknife_cov', 'delete_values'])
    return jknife(est, jknife_se, jknife_est, jknife_var, jknife_cov, delete_values)
end
=#

function _update_weights(
    reg::LD_Score_Regression,
    x_tot, w, N, M_tot, tot_agg, intercept,
)
    if intercept == nothing intercept = reg.__null_intercept__ end
end
