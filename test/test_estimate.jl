"""
Tests for the function `estimate_h2`.
Tests all fail because `ld_score_regression` is only half-implemented.
"""

include("../src/LDScore.jl")

using Main
using Main.LDScore
using Test

approx(x, y, eps) = abs(x - y) <= eps

function vector_approx(x, y, eps)
    if size(x) != size(y)
        return false
    end
    for (i, j) in enumerate(x)
        if !approx(x[i], y[i], eps) return false end
    end
    return true
end

LDScore.args["ref-ld"] = "test/test_ldscore/oneld_onefile"
LDScore.args["w-ld"] = "test/test_ldscore/w"
LDScore.args["h²"] = "test/test_sumstats/1"
LDScore.args["out"] = "test/test_out/1"

# Test case from original ldsc repo
x = LDScore.estimate_h2()

# TODO finish implementing `ld_score_regression`
