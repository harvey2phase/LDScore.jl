using Statistics

"""
Base class for jackknife objects. Input involves x,y, so this base class is tailored
for statistics computed from independent and dependent variables (e.g., regressions).
The __delete_vals_to_pseudovalues__ and __jknife__ methods will still be useful for other
sorts of statistics, but the __init__ method will need to be overriden.


# Attributes
n_blocks : int
    Number of jackknife blocks
p : int
    Dimensionality of the independent varianble
N : int
    Number of datapoints (equal to x.shape[0])

# Methods
-------
jknife(pseudovalues):
    Computes jackknife estimate and variance from the jackknife pseudovalues.
delete_vals_to_pseudovalues(delete_vals, est):
    Converts delete values and the whole-data estimate to pseudovalues.
get_separators():
    Returns (approximately) evenly-spaced jackknife block boundaries.
"""
abstract type Jackknife end


"""
"Constructor" for abstract type `Jackknife`.

# Parameters
- x : np.matrix with shape (n, p)
    Independent variable.
- y : np.matrix with shape (n, 1)
    Dependent variable.
- n_blocks : int
    Number of jackknife blocks
- *args, **kwargs :
    Arguments for inheriting jackknives.
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
- pseudovalues: np.matrix pf floats with shape (n_blocks, p)

# Returns
- jknife_est: np.matrix with shape (1, p)
    Jackknifed estimate.
- jknife_var: np.matrix with shape (1, p)
    Variance of jackknifed estimate.
- jknife_se: np.matrix with shape (1, p)
    Standard error of jackknifed estimate, equal to sqrt(jknife_var).
- jknife_cov: np.matrix with shape (p, p)
    Covariance matrix of jackknifed estimate.
"""
function jknife(cls, pseudovalues)
    n_blocks = size(pseudovalues)[1]
    jknife_cov = atleast_2d(np.cov(pseudovalues.T, ddof=1) / n_blocks)
    jknife_var = np.atleast_2d(np.diag(jknife_cov))
    jknife_se = np.atleast_2d(np.sqrt(jknife_var))
    jknife_est = np.atleast_2d(np.mean(pseudovalues, axis=0))
    return (jknife_est, jknife_var, jknife_se, jknife_cov)
end # function jknife

"""
Converts whole-data estimate and delete values to pseudovalues.

# Arguments
- delete_values: np.matrix with shape (n_blocks, p). Delete values.
- est: np.matrix with shape (1, p). Whole-data estimate.

# Returns
- pseudovalues: np.matrix with shape (n_blocks, p) Psuedovalues.

# Raises
- ValueError: If est.shape != (1, delete_values.shape[1])
"""
function delete_values_to_pseudovalues(cls, delete_values, est)
    n_blocks, p = delete_values.shape
    if est.shape != (1, p)
        error("Different number of parameters in delete_values than in est.")
    end

    return n_blocks * est - (n_blocks - 1) * delete_values
end


"Define evenly-spaced block boundaries."
function get_separators(cls, N, n_blocks)
    return np.floor(np.linspace(0, N, n_blocks + 1)).astype(int)
end


"Check that x and y have the correct shapes (for regression jackknives)."
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
    if xtx_block_values.shape[0:2] != xty_block_values.shape
        error("Shape of xty_block_values must equal shape of first two dimensions of xty_block_values.")
    end
    if len(xtx_block_values.shape) < 3
        error("xtx_block_values must be a 3D array.")
    end
    if xtx_block_values.shape[1] != xtx_block_values.shape[2]
        error("Last two axes of xtx_block_values must have same dimension.")
    end

    return xtx_block_values.shape[0:2]
end # function _check_shape_block
