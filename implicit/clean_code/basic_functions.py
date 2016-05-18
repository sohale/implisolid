import numpy as np
import sys

from implicit_config import TOLERANCE
from implicit_config import VERBOSE


def make_inverse(m):
    assert not issubclass(m.dtype.type, np.integer)
    assert m.shape == (4, 4), "Matrix must be 4x4"
    invm = np.linalg.inv(m)
    assert np.allclose(np.dot(invm, m), np.eye(4), atol=TOLERANCE), "Matrix inversion failed: Matrix is singular or bad conditioned"
    assert np.allclose(np.dot(m, invm), np.eye(4), atol=TOLERANCE), "Matrix inversion failed: Matrix is singular or bad conditioned"
    error = np.sum(np.abs(np.dot(invm, m) - np.eye(4)))
    if VERBOSE:
        print("Error of the inverse matrix: %2.20f" % error)
    v0001 = np.reshape(np.array((0, 0, 0, 1)), (4))
    assert np.allclose(invm[3, :], v0001, atol=TOLERANCE), "Last row of the inverse matrix should be 0,0,0,1 "
    return invm


def check_matrix4(m):
    assert not issubclass(m.dtype.type, np.integer)
    assert m.shape == (4, 4), "Matrix must be 4x4"
    assert np.allclose(m[3, :], np.array((0, 0, 0, 1)), atol=0.00000000001), "Last row of any Matrix4 must be 0,0,0,1"
    assert not np.any(np.isnan(m.ravel()))
    assert not np.any(np.isinf(m.ravel()))


def check_vector4(p):
    assert not issubclass(p.dtype.type, np.integer)
    assert p.shape == (4,), "Vector must be a numpy array of (4) elements"
    assert p[3] == 1.0, "4th element of every Vector must be 1.0"
    assert not np.any(np.isnan(p.ravel()))
    assert not np.any(np.isinf(p.ravel()))


def check_vector3(p):
    assert not issubclass(p.dtype.type, np.integer)
    assert p.shape == (3,), "Vector must be a numpy array of (3) elements"
    assert not np.any(np.isnan(p.ravel()))
    assert not np.any(np.isinf(p.ravel()))


def check_vector2(p):
    assert not issubclass(p.dtype.type, np.integer)
    assert p.shape == (2,), "Vector must be a numpy array of (3) elements"
    assert not np.any(np.isnan(p.ravel()))
    assert not np.any(np.isinf(p.ravel()))


def check_vector4_vectorized(p):
    assert not issubclass(p.dtype.type, np.integer)
    assert p.ndim == 2
    assert p.shape[1:] == (4,), "Vector must be a numpy array of (Nx4) elements"
    e = np.sum(np.abs(p[:, 3] - 1))
    if e > 0.0:
        print("EERROR:", e)
    assert np.allclose(p[:, 3], 1, 0.00000000000001), "4th element of every Vector must be 1.0"
    assert not np.any(np.isnan(p.ravel()))
    assert not np.any(np.isinf(p.ravel()))


def check_vector3_vectorized(p):
    assert not issubclass(p.dtype.type, np.integer)
    assert p.ndim == 2
    assert p.shape[1:] == (3,), "Vector must be a numpy array of (Nx3) elements"
    assert not np.any(np.isnan(p.ravel()))
    assert not np.any(np.isinf(p.ravel()))


def check_scalar_vectorized(v, N=None):
    assert not issubclass(v.dtype.type, np.integer)
    n = v.shape[0]
    if v.ndim == 2:
        assert v.shape[1] == 1
        assert v.shape == (n, 1), "values must be a numpy array of (N,) or (Nx1) elements"
    if N is not None:
        assert v.shape[0] == N
    assert not np.any(np.isnan(v.ravel()))
    assert not np.any(np.isinf(v.ravel()))


def make_vector4(x, y, z):
    xyz = np.array([x, y, z])
    assert not np.any(np.isnan(xyz))
    assert not np.any(np.isinf(xyz))

    if issubclass(type(x), np.ndarray):
        return np.array((float(x[0]), float(y[1]), float(z[1]), 1.0))

    return np.array((float(x), float(y), float(z), 1.0))


def make_vector3(x, y, z):
    xyz = np.array([x, y, z])
    assert not np.any(np.isnan(xyz))
    assert not np.any(np.isinf(xyz))

    if issubclass(type(x), np.ndarray):
        return np.array((float(x[0]), float(y[1]), float(z[1])))

    return np.array((float(x), float(y), float(z)))


def normalize_vector(v, snapToZero=False):
    assert not np.any(np.isnan(v.ravel()))
    assert not np.any(np.isinf(v.ravel()))

    r = v.copy()
    # r[:] = np.sign(r[:]) * np.abs(r[:]) ** POW
    r = r / np.sqrt(np.dot(r, r))
    assert (r[0]*r[0] + r[1]*r[1] + r[2]*r[2] - 1) < 0.00000000001
    if snapToZero:
        for i in range(0, 3):
            if np.abs(r[i]) < 0.0000001:
                r[i] = 0
    r[3] = 1
    return r


def normalize_vector3(v, snapToZero=False):
    assert not np.any(np.isnan(v.ravel()))
    assert not np.any(np.isinf(v.ravel()))

    r = v.copy()
    r = r / np.sqrt(np.dot(r, r))
    assert (r[0]*r[0] + r[1]*r[1] + r[2]*r[2] - 1) < 0.00000000001
    if snapToZero:
        for i in range(0, 3):
            if np.abs(r[i]) < 0.0000001:
                r[i] = 0
    return r


def make_random_vector(norm, POW, type="rand"):
    if type == "rand":
        r = np.random.rand(3)*2 - 1
    elif type == "randn":
        r = np.random.randn(3)
    else:
        raise ValueError("nknown random distribution")

    r[:] = np.sign(r[:]) * np.abs(r[:]) ** POW
    r = r / np.sqrt(np.dot(r, r))
    assert (r[0]*r[0] + r[1]*r[1] + r[2]*r[2] - 1) < 0.00000000001
    r = r * norm
    for i in range(0, 3):
        if np.abs(r[i]) < 0.0000001:
            r[i] = 0
    return np.array((r[0], r[1], r[2], 1))


def make_random_vector3(norm, POW, type="rand"):
    if type == "rand":
        r = np.random.rand(3)*2 - 1
    elif type == "randn":
        r = np.random.randn(3)
    else:
        raise ValueError("nknown random distribution")

    r[:] = np.sign(r[:]) * np.abs(r[:]) ** POW
    r = r / np.sqrt(np.dot(r, r))
    assert (r[0]*r[0] + r[1]*r[1] + r[2]*r[2] - 1) < 0.00000000001
    r = r * norm
    for i in range(0, 3):
        if np.abs(r[i]) < 0.0000001:
            r[i] = 0
    return np.array((r[0], r[1], r[2]))


def make_random_vector_vectorized(N, norm, POW, type="rand", normalize=True):
    r = np.ones((N, 4))

    # not tested:
    if type == "rand":
        r[:, 0:3] = np.random.rand(N, 3)*2 - 1
    elif type == "randn":
        r[:, 0:3] = np.random.randn(N, 3)
    else:
        raise ValueError("nknown random distribution")

    r[:, 0:3] = np.sign(r[:, 0:3]) * np.abs(r[:, 0:3]) ** POW
    n3 = np.sqrt(np.sum(r[:, 0:3] * r[:, 0:3], axis=1, keepdims=True))
    if normalize:
        r[:, 0:3] = r[:, 0:3] / np.tile(n3, (1, 3))
        s_1 = np.sum(r[:, 0:3] * r[:, 0:3], axis=1) - 1
        assert np.all(np.abs(s_1) < 0.00000000001)
    r[:, 0:3] = r[:, 0:3] * norm
    r[:, 3] = 1
    check_vector4_vectorized(r)
    return r


def make_random_vector3_vectorized(N, norm, POW, type="rand", normalize=True):
    r = np.ones((N, 3))
    # not tested:
    if type == "rand":
        r = np.random.rand(N, 3)*2 - 1
    elif type == "randn":
        r = np.random.randn(N, 3)
    else:
        raise ValueError("nknown random distribution")

    r = np.sign(r) * np.abs(r) ** POW
    n3 = np.sqrt(np.sum(r * r, axis=1, keepdims=True))
    if normalize:
        r = r / np.tile(n3, (1, 3))
        s_1 = np.sum(r * r, axis=1) - 1
        assert np.all(np.abs(s_1) < 0.00000000001)
    r = r * norm
    check_vector3_vectorized(r)
    return r


def normalize_vector4_vectorized(v, zero_normal="leave_zero_norms"):
    """ returns vectors of either length 1 or zero. """
    N = v.shape[0]
    assert not issubclass(v.dtype.type, np.integer)
    assert not np.any(np.isnan(v))
    assert not np.any(np.isinf(v))

    norms = np.sqrt(np.sum(v[:, 0:3]*v[:, 0:3], axis=1, keepdims=True))
    denominator = np.tile(norms, (1, 4))
    if zero_normal == "leave_zero_norms":
        zeros_i = np.abs(norms.ravel()) < 0.00000001
        non_zero_i = np.logical_not(zeros_i)
        if not np.any(zeros_i):
            c = 1.0 / denominator
        else:
            c = np.ones(denominator.shape)
            c[non_zero_i, :] = 1.0 / denominator[non_zero_i, :]
    else:
        pass
    assert not np.any(np.isnan(c))
    assert not np.any(np.isinf(c))
    r = v * c
    assert r.shape[0] == N
    df = np.sum(r[:, 0:3] * r[:, 0:3], axis=1)
    e1a = np.all(np.abs(df[non_zero_i]-1.0) < 0.00000000001)
    e0a = np.all(np.abs(df[zeros_i]) < 0.00000000001)

    if not (e1a and e0a):
        print("r:", r)
        print("v:", v)
        print("c:", v)
        print("denom: ", denominator)
        print(norms)
        print(denominator)
        print(np.sum(r[:, 0:3] * r[:, 0:3], axis=1))

    assert e1a and e0a
    r[:, 3] = 1
    return r


def normalize_vector3_vectorized(v, zero_normal="leave_zero_norms"):
    """ returns vectors of either length 1 or zero. """
    N = v.shape[0]
    assert not issubclass(v.dtype.type, np.integer)
    assert not np.any(np.isnan(v))
    assert not np.any(np.isinf(v))
    norms = np.sqrt(np.sum(v * v, axis=1, keepdims=True))
    denominator = np.tile(norms, (1, 3))
    if zero_normal == "leave_zero_norms":
        zeros_i = np.abs(norms.ravel()) < 0.00000001
        non_zero_i = np.logical_not(zeros_i)
        if not np.any(zeros_i):
            c = 1.0 / denominator
        else:
            c = np.ones(denominator.shape)
            c[non_zero_i, :] = 1.0 / denominator[non_zero_i, :]
    else:
        pass
    assert not np.any(np.isnan(c))
    assert not np.any(np.isinf(c))
    r = v * c
    assert r.shape[0] == N
    df = np.sum(r * r, axis=1)
    e1a = np.all(np.abs(df[non_zero_i]-1.0) < 0.00000000001)
    e0a = np.all(np.abs(df[zeros_i]) < 0.00000000001)

    if not (e1a and e0a):
        print("r:", r)
        print("v:", v)
        print("c:", v)
        print("denom: ", denominator)
        print(norms)
        print(denominator)
        print(np.sum(r*r, axis=1))

    assert e1a and e0a
    return r


def repeat_vect4(N, v4):
    check_vector4(v4)
    _x = v4
    xa = np.tile(np.expand_dims(_x, axis=0), (N, 1))
    assert xa.shape[0] == N
    return xa


def repeat_vect3(N, v3):
    check_vector3(v3)
    _x = v3
    xa = np.tile(np.expand_dims(_x, axis=0), (N, 1))
    assert xa.shape[0] == N
    return xa


def almost_equal1(a, b, TOLERANCE):
    """ Scalar version """
    assert not issubclass(a.dtype.type, np.integer)
    np.isscalar(a)
    np.isscalar(b)
    return np.sum(np.abs(a - b)) < TOLERANCE


def is_python3():

    v = sys.version_info.major
    return v == 3
