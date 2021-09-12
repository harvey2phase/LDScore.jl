""" Basic tests for the `Hsq` struct and associated function. """

include("../src/LDScore.jl")
include("test_tools.jl")


using Main
using Main.LDScore
using Test

"Accuracy of approximatly equal; the lower the more accurate"
eps = 1.0e-6

#@testset "Ensure basic Hsq struct works" begin
    # column vectors d = 6 * 1
    y = [[0.2061] [0.2601] [4.3514] [6.1703] [5.0221] [2.418]]'
    x = [[6.1388] [6.4785] [8.88] [4.8064] [3.6219] [3.7244]]'
    w = [[6.1388] [6.4785] [8.88] [4.8064] [3.6219] [3.7244]]'
    N = [[592497.0] [591333.0] [696881.0] [617319.0] [687344.0] [685637.0]]'
    M = [[1173569.0]]
    n_blocks = 2
    intercept = nothing
    slow = false
    old_weights = false

    hsq = LDScore.Hsq(
        y, x, w, N, M, n_blocks = n_blocks, intercept = intercept, slow = slow,
        old_weights = old_weights,
    )

    hsq.old_weights == false

    @test typeof(hsq) == LDScore.Hsq
    @test hsq.__null_intercept__ == 1
    @test hsq.y == y
    @test hsq.x == x
    @test hsq.w == w
    @test hsq.N == N
    @test hsq.M == M
    @test hsq.n_blocks == n_blocks
    @test hsq.intercept == intercept
    @test hsq.slow == slow
#end


#@testset "Hsq functions" begin
    χ² = ones((4, 1)) .* 4
    ld = ones((4, 1))
    w_ld = ones((4, 1))
    N = 9 .* ones((4, 1))
    M = [[7]]
    hsq = LDScore.Hsq(
        χ², ld, w_ld, N, M, n_blocks = 3, intercept = 1,
        slow = false, old_weights = false, twostep = nothing,
    )

    @testset "weights" begin
        x = [[0.09570313] [0.09570313] [0.09570313] [0.09570313]]'
        a = LDScore.hsq_weights(ld, w_ld, N, M, 1)
        b = LDScore.hsq_weights(ld, w_ld, N, M, 2)
        @test vector_approx(a, x, eps)
        @test vector_approx(a, b, eps)

        x = [[0.5] [0.5] [0.5] [0.5]]'
        a = LDScore.hsq_weights(ld, w_ld, N, M, 0)
        b = LDScore.hsq_weights(ld, w_ld, N, M, -1)
        @test vector_approx(a, x, eps)
        @test vector_approx(a, b, eps)
    end

    @testset "aggregate" begin
        χ² = ones((10, 1)) .* 3 ./ 2
        ld = ones((10, 1)) .* 100
        N = ones((10, 1)) .* 100000
        M = 1e7
        agg = LDScore.aggregate(hsq, χ², ld, N, M)
        approx(agg, 0.5, eps)
        agg = LDScore.aggregate(hsq, χ², ld, N, M; intercept = 1.5)
        approx(agg, 0, eps)
    end
#end
