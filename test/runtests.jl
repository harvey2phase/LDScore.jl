using Main
using Main.LDScoreJulia
using Test

@testset "LDScoreJulia.jl" begin
    # Write your tests here.
end

# synthetic dataset
# column vectors d = 6 * 1
y = [[0.2061] [0.2601] [4.3514] [6.1703] [5.0221] [2.418 ]]'
x = [[6.1388] [6.4785] [8.88] [4.8064] [3.6219] [3.7244]]'
w = [[6.1388] [6.4785] [8.88] [4.8064] [3.6219] [3.7244]]'
N = [[592497.] [591333.] [696881.] [617319.] [687344.] [685637.]]'
M = [[1173569.]]
n_blocks = 200
intercept = nothing
slow = false
step1_ii = [[true] [true] [true] [true] [true] [true]]
old_weights = false
synthetic_hsq = Hsq(
    y, x, w, N, M, n_blocks, intercept, slow, step1_ii, old_weights,
)

@test typeof(synthetic_hsq) == Hsq
