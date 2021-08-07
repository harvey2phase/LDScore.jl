using Main
using Main.LDScoreJulia
using Test

println("==================== NEW ===========================")

# TODO begin testset ===========================================================

# column vectors d = 6 * 1
y = [[0.2061] [0.2601] [4.3514] [6.1703] [5.0221] [2.418 ]]'
x = [[6.1388] [6.4785] [8.88] [4.8064] [3.6219] [3.7244]]'
w = [[6.1388] [6.4785] [8.88] [4.8064] [3.6219] [3.7244]]'
N = [[592497.] [591333.] [696881.] [617319.] [687344.] [685637.]]'
M = [[1173569.]]
n_blocks = 2
intercept = nothing
slow = false
step1_ii = [[true] [true] [true] [true] [true] [true]]'
old_weights = false

synthetic_hsq = LDScoreJulia.Hsq(
    y, x, w, N, M, n_blocks, intercept, slow, step1_ii, old_weights,
)

synthetic_hsq.old_weights == false
LDScoreJulia.ld_score_regression(
    synthetic_hsq,
    y, x, w, N, M, n_blocks, intercept, slow, step1_ii, old_weights,
)

@test typeof(synthetic_hsq) == LDScoreJulia.Hsq
@test synthetic_hsq.__null_intercept__ == 1
@test synthetic_hsq.y == y
@test synthetic_hsq.x == x
@test synthetic_hsq.w == w
@test synthetic_hsq.N == N
@test synthetic_hsq.M == M
@test synthetic_hsq.n_blocks == n_blocks
@test synthetic_hsq.intercept == intercept
@test synthetic_hsq.slow == slow
@test synthetic_hsq.step1_ii == step1_ii

# end test set -----------------------------------------------------------------

#=
# TODO begin test set ==========================================================

chisq = ones((4, 1)) .* 4
ld = ones((4, 1))
w_ld = ones((4, 1))
N = 9 .* ones((4, 1))
M = matrix((7))
hsq = LDScoreJulia.Hsq(chisq, ld, w_ld, N, M, n_blocks=3, intercept=1)

@testset "Hsq weights" begin

end

# end test set -----------------------------------------------------------------
=#
