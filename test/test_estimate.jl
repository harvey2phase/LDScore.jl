"""
Tests for the `estimate_h2` function.
(Tests all fail because `ld_score_regression` is only half-implemented.)
"""

include("../src/LDScore.jl")
include("test_tools.jl")

using Main
using Main.LDScore
using Test


LDScore.args["ref-ld"] = "test/test_ldscore/oneld_onefile"
LDScore.args["w-ld"] = "test/test_ldscore/w"
LDScore.args["h²"] = "test/test_sumstats/1"
LDScore.args["out"] = "test/test_out/1"

# Test case from original ldsc repo
x = LDScore.estimate_h2()

# TODO finish implementing `ld_score_regression`
