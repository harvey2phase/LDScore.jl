using Statistics

include("LD_Score_Regression.jl")
include("Hsq.jl")

# TODO improve error-handeling to the Julian way

function fmax(col_vec, val)
    new = zeros(0)
    for a in col_vec append!(new, max(a, val)) end
    return reshape(new, (size(new)...,1))
end

# TODO probably move functions with reg::Hsq to Hsq.jl

# TODO add a test (works though)
function weights(
    reg::Hsq,
    ld, w_ld, N, M, hsq, intercept,
)
    if intercept == nothing intercept = 1 end
    hsq = max(maximum(hsq), 0.0)
    hsq = min(hsq, 1.0)
    ld = fmax(ld, 1.0)
    w_ld = fmax(w_ld, 1.0)
    c = hsq .* N ./ M

    het_w = 1.0 ./ (2 .* (c .* ld .+ intercept).^2)
    oc_w = 1.0 ./ w_ld
    w = het_w .* oc_w
    println("w")
    println(w)
    return w
end

function _update_weights(
    reg::Hsq,
    ld, w_ld, N, M, hsq, intercept,
)
    if intercept == nothing intercept = reg.__null_intercept__ end

    return weights(reg, ld, w_ld, N, M, hsq, intercept)
end

function aggregate(
    reg::LD_Score_Regression,
    y, x, N, M, intercept,
)
    if intercept == nothing
        intercept = reg.__null_intercept__
    end

    num = M * (mean(y) - intercept)
    println("num")
    println(num)
    denom = mean(x .* N)
    return num / denom
end


function make_ld_score_regression(
    reg,
    y, x, w, N, M, n_blocks, intercept, slow, step1_ii, old_weights,
)
    #= TODO
    for i in [y, x, w, M, N]:
        try:
            if len(i.shape) != 2:
                raise TypeError('Arguments must be 2D arrays.')
        except AttributeError:
            raise TypeError('Arguments must be arrays.')
    =#

    n_snp, reg.n_annot = size(x)
    #= TODO
    if any(i.shape != (n_snp, 1) for i in [y, w, N]):
        raise ValueError(
            'N, weights and response (z1z2 or chisq) must have shape (n_snp, 1).')
    if M.shape != (1, self.n_annot):
        raise ValueError('M must have shape (1, n_annot).')
    =#

    M_tot = sum(M)
    # shape should be [n_snp, 1]
    x_tot = sum(x, dims = 2)
    reg.constrain_intercept = intercept != nothing
    reg.intercept = intercept
    reg.n_blocks = n_blocks

    tot_agg = aggregate(reg, y, x_tot, N, M_tot, intercept)
    initial_w = _update_weights(
        reg, x_tot, w, N, M_tot, tot_agg, intercept,
    )
    Nbar = mean(N)
    x = N .* x / Nbar
    println("X AFTER A BIT")
    println(x)

    if !reg.constrain_intercept
        x, x_tot = append_intercept(x), append_intercept(x_tot)
        println("x_tot")
        println(x_tot)
        yp = y
    else
        yp = y - intercept
        intercept_se = "NA"
    end

    reg.twostep_filtered = nothing
    if !(step1_ii == nothing) && reg.constrain_intercept
        throw(ErrorException("twostep is not compatible with constrain_intercept."))
    elseif !(step1_ii == nothing) && reg.n_annot > 1
        throw(ErrorException("twostep not compatible with partitioned LD Score yet."))
    elseif !(step1_ii == nothing)
        n1 = sum(step1_ii)
        reg.twostep_filtered = n_snp - n1
        s = dropdims(step1_ii; dims=2)
        x1 = zeros(0)
        for (i, j) in enumerate(s)
            if s[i] == 1
                append!(x1, x[i])
            end
        end
        #=
        yp1, w1, N1, initial_w1 = map(
            lambda a: a[step1_ii].reshape((n1, 1)), (yp, w, N, initial_w))
            update_func1 = lambda a: self._update_func(
            a, x1, w1, N1, M_tot, Nbar, ii=step1_ii)
            step1_jknife = IRWLS(
            x1, yp1, update_func1, n_blocks, slow=slow, w=initial_w1)
            step1_int, _ = self._intercept(step1_jknife)
            yp = yp - step1_int
            x = remove_intercept(x)
            x_tot = remove_intercept(x_tot)
            update_func2 = lambda a: self._update_func(
            a, x_tot, w, N, M_tot, Nbar, step1_int)
            s = update_separators(step1_jknife.separators, step1_ii)
            step2_jknife = IRWLS(
            x, yp, update_func2, n_blocks, slow=slow, w=initial_w, separators=s)
            c = np.sum(np.multiply(initial_w, x)) / \
            np.sum(np.multiply(initial_w, np.square(x)))
            jknife = self._combine_twostep_jknives(
            step1_jknife, step2_jknife, M_tot, c, Nbar)
        =#
    end
end

function append_intercept(x)
    n_row = size(x)[1]
    intercept = ones(n_row, 1)
    x_new = [x; intercept]
    return x_new
end
