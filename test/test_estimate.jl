#using Pkg; Pkg.add("Revise")
#using Revise

include("../src/LDScoreJulia.jl")

using Main
using Main.LDScoreJulia
using Test

approx(x, y, eps) = abs(x - y) <= eps

function vector_approx(x, y, eps)
    if size(x) != size(y)
        return false
    end
    for (i, j) in enumerate(x)
        if !approx(x[i], y[i], eps)
            return false
        end
    end
    return true
end

# Test case from original ldsc repo
LDScoreJulia.estimate_h2(
    Dict([
        ("ref_id", "LDScoreJulia/test/test_ldscore/oneld_onefile"),
        ("w_ld", "LDScoreJulia/test/test_ldscore/w"),
        ("h2", "LDScoreJulia/test/test_sumstats/1"),
        ("out", "LDScoreJulia/test/test_out/1"),

        # default parser arguments
        ("n_blocks", 200),
        ("χ²_max", nothing),
        ("two_step", nothing),
        ("intercept_h2", nothing),
    ])
)

def _overlap_output(self, category_names, overlap_matrix, M_annot, M_tot, print_coefficients):
    '''LD Score regression summary for overlapping categories.'''
    overlap_matrix_prop = np.zeros([self.n_annot,self.n_annot])
    for i in range(self.n_annot):
        overlap_matrix_prop[i, :] = overlap_matrix[i, :] / M_annot

    prop_hsq_overlap = np.dot(
        overlap_matrix_prop, self.prop.T).reshape((1, self.n_annot))
    prop_hsq_overlap_var = np.diag(
        np.dot(np.dot(overlap_matrix_prop, self.prop_cov), overlap_matrix_prop.T))
    prop_hsq_overlap_se = np.sqrt(
        np.maximum(0, prop_hsq_overlap_var)).reshape((1, self.n_annot))
    one_d_convert = lambda x: np.array(x).reshape(np.prod(x.shape))
    prop_M_overlap = M_annot / M_tot
    enrichment = prop_hsq_overlap / prop_M_overlap
    enrichment_se = prop_hsq_overlap_se / prop_M_overlap
    overlap_matrix_diff = np.zeros([self.n_annot,self.n_annot])
    for i in range(self.n_annot):
        if not M_tot == M_annot[0,i]:
            overlap_matrix_diff[i, :] = overlap_matrix[i,:]/M_annot[0,i] - \
                (M_annot - overlap_matrix[i,:]) / (M_tot-M_annot[0,i])

    diff_est = np.dot(overlap_matrix_diff,self.coef)
    diff_cov = np.dot(np.dot(overlap_matrix_diff,self.coef_cov),overlap_matrix_diff.T)
    diff_se = np.sqrt(np.diag(diff_cov))
    diff_p = ['NA' if diff_se[i]==0 else 2*tdist.sf(abs(diff_est[i]/diff_se[i]),self.n_blocks) \
        for i in range(self.n_annot)]

    df = pd.DataFrame({
        'Category': category_names,
        'Prop._SNPs': one_d_convert(prop_M_overlap),
        'Prop._h2': one_d_convert(prop_hsq_overlap),
        'Prop._h2_std_error': one_d_convert(prop_hsq_overlap_se),
        'Enrichment': one_d_convert(enrichment),
        'Enrichment_std_error': one_d_convert(enrichment_se),
        'Enrichment_p':diff_p,
        'Coefficient': one_d_convert(self.coef),
        'Coefficient_std_error': self.coef_se,
        'Coefficient_z-score': one_d_convert(self.coef) / one_d_convert(self.coef_se)
    })
    if print_coefficients:
        df = df[['Category', 'Prop._SNPs', 'Prop._h2', 'Prop._h2_std_error',
                'Enrichment','Enrichment_std_error', 'Enrichment_p',
                 'Coefficient', 'Coefficient_std_error','Coefficient_z-score']]
    else:
        df = df[['Category', 'Prop._SNPs', 'Prop._h2', 'Prop._h2_std_error',
                'Enrichment','Enrichment_std_error', 'Enrichment_p']]
    return df
