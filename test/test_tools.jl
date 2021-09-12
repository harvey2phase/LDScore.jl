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
