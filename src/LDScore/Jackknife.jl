using Statistics

"""
Base class for jackknife objects. Input involves x,y, so this base class is tailored
for statistics computed from independent and dependent variables (e.g., regressions).
The __delete_vals_to_pseudovalues__ and __jknife__ methods will still be useful for other
sorts of statistics, but the __init__ method will need to be overriden.

# Attributes
- `n_blocks::Integer`: Number of jackknife blocks.
- `p::Integer`: Dimensionality of the independent varianble.
- `N::Integer`: Number of datapoints (equal to `size(x)[1]`).

# Functions
- `jknife(pseudovalues)`:
    Computes jackknife estimate and variance from the jackknife pseudovalues.
- `delete_vals_to_pseudovalues(delete_vals, est)`:
    Converts delete values and the whole-data estimate to pseudovalues.
- `get_separators()`:
    Returns (approximately) evenly-spaced jackknife block boundaries.
"""
abstract type Jackknife end


"""
"Constructor" for abstract type `Jackknife`.

# Parameters
- `x::Matric`: Independent variable with shape `(n, p)`.
- `y::Matrix`: Dependent variable with shape `(n, 1)`.
- `n_blocks::Integer`: Number of jackknife blocks.
"""
function jackknife(x, y; n_blocks = nothing, separators = nothing)
    self.N, self.p = _check_shape(x, y)
    if separators != nothing
        if max(separators) != self.N
            error("Max(separators) must be equal to number of data points.")
        end
        if min(separators) != 0
            error("Max(separators) must be equal to 0.")
        end
        self.separators = sort(separators)
        self.n_blocks = length(separators) - 1
    elseif n_blocks != nothing
        self.n_blocks = n_blocks
        self.separators = self.get_separators(self.N, self.n_blocks)
    else
        error("Must specify either n_blocks are separators.")
    end

    if self.n_blocks > self.N
        error("More blocks than data points.")
    end
end # function jackknife


"""
Converts pseudovalues to jackknife estimate and variance.

# Arguments
- `pseudovalues::Matrix`: pf floats with shape `(n_blocks, p)`.

# Returns
- `jknife_est::Matrix`: Jackknifed estimate with shape `(1, p)`.
- `jknife_var::Matrix`: Variance of jackknifed estimate with shape `(1, p)`.
- `jknife_se::Matrix`: Standard error of jackknifed estimate, equal to `sqrt(jknife_var)` with shape `(1, p)`.
- `jknife_cov::Matrix`: Covariance matrix of jackknifed estimate with shape `(p, p)`.
"""
function jknife(cls, pseudovalues)
    n_blocks = size(pseudovalues)[1]
    jknife_cov = atleast_2d(cov(pseudovalues.T, ddof = 1) / n_blocks)
    jknife_var = atleast_2d(diag(jknife_cov))
    jknife_se = atleast_2d(sqrt(jknife_var))
    jknife_est = atleast_2d(mean(pseudovalues, axis=0))
    return (jknife_est, jknife_var, jknife_se, jknife_cov)
end # function jknife

"""
Converts whole-data estimate and delete values to pseudovalues.
Raises error if `size(est) != (1, size(delete_values)[2])`

# Arguments
- `delete_values::Matrix`: Delete values with shape (of the matrix) `(n_blocks, p)`.
- `est::Matrix`: Whole-data estimate with shape `(1, p)`.

# Returns
- `pseudovalues::Matrix`: Psuedovalues with shape `(n_blocks, p)`.
"""
function delete_values_to_pseudovalues(cls, delete_values, est)
    n_blocks, p = size(delete_values)
    if size(est) != (1, p)
        error("Different number of parameters in delete_values than in est.")
    end

    return n_blocks * est - (n_blocks - 1) * delete_values
end


"Define evenly-spaced block boundaries."
function get_separators(cls, N, n_blocks)
    return convert(Vector{Int}, floor.(range(0, N, n_blocks + 1)))
end


"Check that `x` and `y` have the correct shapes (for regression jackknives)."
function _check_shape(x, y)
    if length(size(x)) != 2 || length(size(y)) != 2
        error("x and y must be 2D arrays.")
    end
    if size(x)[1] != size(y)[1]
        error("Number of datapoints in x != number of datapoints in y.")
    end
    if size(y)[2] != 1
        error("y must have shape (n_snp, 1)")
    end

    n, p = size(x)
    if p > n
        error("More dimensions than datapoints.")
    end

    return (n, p)
end # function _check_shape


function _check_shape_block(xty_block_values, xtx_block_values)
    """Check that xty_block_values and xtx_block_values have correct shapes."""
    if size(xtx_block_values)[1:3] != size(xty_block_values)
        error("Shape of xty_block_values must equal shape of first two dimensions of xty_block_values.")
    end
    if length(size(xtx_block_values)) < 3
        error("xtx_block_values must be a 3D array.")
    end
    if size(xtx_block_values)[2] != size(xtx_block_values)[3]
        error("Last two axes of xtx_block_values must have same dimension.")
    end

    return size(xtx_block_values)[1:3]
end # function _check_shape_block
