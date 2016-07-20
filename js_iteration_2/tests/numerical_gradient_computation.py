import numpy as np

NUMERICAL_GRADIENT_TOLERANCE = 0.0001  # 0.00001   # 0.001
assert NUMERICAL_GRADIENT_TOLERANCE > 0.0000059
# assert NUMERICAL_GRADIENT_TOLERANCE > 0.00000001


def make_vector3(x, y, z):
    xyz = np.array([x, y, z])
    assert not np.any(np.isnan(xyz))
    assert not np.any(np.isinf(xyz))

    if issubclass(type(x), np.ndarray):
        return np.array((float(x[0]), float(y[1]), float(z[1])))

    return np.array((float(x), float(y), float(z)))


def repeat_vect3(N, v3):
    _x = v3
    xa = np.tile(np.expand_dims(_x, axis=0), (N, 1))
    assert xa.shape[0] == N
    return xa


def weights(k, x0, xs):
    """Calculate weights for the finite difference approximation.
    Arguments:
    k - The k-th derivative will be approximated.
    x0 - The point at which to approximate the derivative.
    xs - The grid points at which the function's value is known.
    This uses the algorithm described in:
    B. Fornberg, "Calculation of weights in finite difference formulas",
    SIAM Review 40 (1998), pp. 685-691.
    """
    # Size of the system.
    n = xs.size
    assert k < n, "need more grid points to calculate this derivative"
    # Measure points relative to x0.
    xs = xs - x0
    # Weight matrix to calculate successive finite difference
    # approximations.
    w = np.zeros((k+1, n))
    # Before starting, we want to pre-compute certain reusable
    # quantities.
    product_ratio = np.ones(n)
    for j in range(n):
        for i in range(j):
            product_ratio[j] *= xs[j] - xs[i]
    product_ratio[1:] = product_ratio[:n-1] / product_ratio[1:]
    # Each iteration of this loop will produce weights for the
    # derivatives that use one point more than the previous iteration.
    for j in range(n - k):
        # 0-th derivative approximation.
        if j == 0:
            # The approximation to the 0th derivative given only one
            # function value is trivially to just use that value.
            w[0, 0] = 1
        else:
            w[0, j] = - xs[j-1] * w[0, j-1] * product_ratio[j]
            for i in range(1, j+1):
                w[0, j-i] = xs[j] * w[0, j-i] / (xs[j] - xs[j-i])
        for m in range(1, k+1):
            # Generate weights for each derivative using the
            # previous one.
            # m is the derivative we are currently working on,
            # and l is the number of points used in this round,
            # minus one.
            l = j+m
            w[m, l] = (m*w[m-1, l-1] - xs[l-1]*w[m, l-1]) * product_ratio[l]
            for i in range(1, l+1):
                w[m, l-i] = (xs[l] * w[m, l-i] - m*w[m-1, l-i]) \
                            / (xs[l] - xs[l-i])
    return w[k, :]


def numerical_gradient(pos0, delta_t=0.01/10.0/10.0, order=5):

    m = order  # sample points: -m,...,-1,0,1,2,...,+m

    sample_points = range(-m, m+1)
    n = m*2+1

    x0 = 0
    findiff_weights = weights(k=1, x0=x0, xs=np.array(sample_points) * delta_t)

    pos = repeat_vect3(1, pos0)
    pos3 = np.tile(pos, (3*n, 1))

    assert not issubclass(pos.dtype.type, np.integer)

    dx = repeat_vect3(1, make_vector3(1, 0, 0))
    dy = repeat_vect3(1, make_vector3(0, 1, 0))
    dz = repeat_vect3(1, make_vector3(0, 0, 1))
    dxyz = [dx, dy, dz]

    ci = 0
    for d in range(3):
        for i in sample_points:
            dd = dxyz[d] #we need to change the value of pos3
            pos3[ci, :] = pos3[ci, :] + (dd * delta_t * float(i))
            ci += 1
    v = np.array([-0.8, -0.8, -0.8, -0.8, -0.8, -0.64, -0.48, -0.48, -0.48, -0.48, -0.48, -0.48, -0.48, -0.48, -0.48, -0.48, -0.48, -0.48, -0.48, -0.48, -0.48, -0.48, -0.48, -0.48, -0.48, -0.8, -0.64, -0.64, -0.64, -0.64, -0.64, -0.32, -0.32])
    # v = (-0.8)*np.ones(pos3.shape[0])
    # v = iobj.implicitFunction(pos3) #replace by the text file

    v3 = np.reshape(v, (3, n), order='C')

    """ Calculating the numerical derivative using finite difference (convolution with weights) """
    # convolusion
    grad_cnv = np.dot(v3, findiff_weights)
    print(grad_cnv)
    # Detecting sharp edges (non-smooth points, i.e. corners and edges and ridges)
    if np.max(np.abs(grad_cnv)) > 100:
        pass

    return grad_cnv.reshape(1, 3)


if __name__ == '__main__':
    """ numerical """
    x = make_vector3(0.2, 0, 0)
    g2 = numerical_gradient(x)
