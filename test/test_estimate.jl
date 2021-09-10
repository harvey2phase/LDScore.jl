"""
Tests for the function `estimate_h2`.
Tests all fail because `ld_score_regression` is only half-implemented.
"""

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

print(LDScoreJulia.args)
LDScoreJulia.args["ref-ld"] = "test/test_ldscore/oneld_onefile"
LDScoreJulia.args["w-ld"] = "test/test_ldscore/w"
LDScoreJulia.args["h²"] = "test/test_sumstats/1"
LDScoreJulia.args["out"] = "test/test_out/1"

# Test case from original ldsc repo
x = LDScoreJulia.estimate_h2()

# TODO finish implementing `ld_score_regression`
