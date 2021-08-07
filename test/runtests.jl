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

@testset "Hsq struct" begin
    # column vectors d = 6 * 1
    y = [[0.2061] [0.2601] [4.3514] [6.1703] [5.0221] [2.418]]'
    x = [[6.1388] [6.4785] [8.88] [4.8064] [3.6219] [3.7244]]'
    w = [[6.1388] [6.4785] [8.88] [4.8064] [3.6219] [3.7244]]'
    N = [[592497.0] [591333.0] [696881.0] [617319.0] [687344.0] [685637.0]]'
    M = [[1173569.0]]
    n_blocks = 2
    intercept = nothing
    slow = false
    step1_ii = [[true] [true] [true] [true] [true] [true]]'
    old_weights = false

    hsq = LDScoreJulia.Hsq(
        y, x, w, N, M, n_blocks, intercept, slow, step1_ii, old_weights,
    )

    hsq.old_weights == false
    LDScoreJulia.ld_score_regression(
        hsq, y, x, w, N, M, n_blocks, intercept, slow, step1_ii, old_weights,
    )

    @test typeof(hsq) == LDScoreJulia.Hsq
    @test hsq.__null_intercept__ == 1
    @test hsq.y == y
    @test hsq.x == x
    @test hsq.w == w
    @test hsq.N == N
    @test hsq.M == M
    @test hsq.n_blocks == n_blocks
    @test hsq.intercept == intercept
    @test hsq.slow == slow
    @test hsq.step1_ii == step1_ii
end


@testset "Hsq functions" begin
    eps = 1.0e-6

    chisq = ones((4, 1)) .* 4
    ld = ones((4, 1))
    w_ld = ones((4, 1))
    N = 9 .* ones((4, 1))
    M = [[7]]
    n_blocks = 3
    intercept = 1
    slow = false
    twostep = nothing
    old_weights = false
    hsq = LDScoreJulia.Hsq(
        chisq, ld, w_ld, N, M, n_blocks, intercept, slow, twostep, old_weights,
    )

    @testset "weights" begin
        intercept = nothing

        x = [[0.09570313] [0.09570313] [0.09570313] [0.09570313]]'
        a = hsq_weights(ld, w_ld, N, M, 1, intercept)
        b = hsq_weights(ld, w_ld, N, M, 2, intercept)
        @test vector_approx(a, x, eps)
        @test vector_approx(a, b, eps)

        x = [[0.5] [0.5] [0.5] [0.5]]'
        a = hsq_weights(ld, w_ld, N, M, 0, intercept)
        b = hsq_weights(ld, w_ld, N, M, -1, intercept)
        @test vector_approx(a, x, eps)
        @test vector_approx(a, b, eps)
    end

    @testset "aggregate" begin
        chisq = ones((10, 1)) .* 3 ./ 2
        ld = ones((10, 1)) .* 100
        N = ones((10, 1)) .* 100000
        M = 1e7
        agg = aggregate(hsq, chisq, ld, N, M)
        approx(agg, 0.5, eps)
        agg = aggregate(hsq, chisq, ld, N, M; intercept=1.5)
        approx(agg, 0, eps)
    end
end

function f(x; y=3)
    x + y
end

f(3; y=9)
