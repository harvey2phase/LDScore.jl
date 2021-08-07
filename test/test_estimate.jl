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


estimate_h2(
    Dict([
        ("h2", "height.sumstats.gz"),
        ("ref-ld-chr", "eur_w_ld_chr/"),
        ("w-ld-chr", "eur_w_ld_chr/"),
        ("out", "height_h2"),
    ])
)
