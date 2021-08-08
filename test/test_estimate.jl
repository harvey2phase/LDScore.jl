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
x = LDScoreJulia.estimate_h2(
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
        ("overlap_annot", false),
    ])
)
print(x.tot)
