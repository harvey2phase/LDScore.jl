using Statistics

include("LD_Score_Regression.jl")
include("Hsq.jl")

# TODO improve error-handeling (the Julian way)

function fmax(col_vec, val)
    new = zeros(0)
    for a in col_vec append!(new, max(a, val)) end
    return reshape(new, (size(new)...,1))
end


function append_intercept(x)
    n_row = size(x)[1]
    intercept = ones(n_row, 1)
    x_new = [x; intercept]
    return x_new
end
