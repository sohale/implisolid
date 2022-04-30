import numpy as np

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
    assert not np.any( np.isnan(m.ravel()) )
    assert not np.any( np.isinf(m.ravel()) )

def check_matrix4_vectorized(m):
    assert not issubclass(m.dtype.type, np.integer)
    assert m.shape[1, :] == (4, 4), "Matrix must be Nx4x4"
    for i in range(m.shape[0]):
        assert np.allclose(m[i, 3, :], np.array((0, 0, 0, 1)), atol=0.00000000001), "Last row of any Matrix4 must be 0,0,0,1"
    assert not np.any( np.isnan(m.ravel()) )
    assert not np.any( np.isinf(m.ravel()) )


def check_matrix3_vectorized(h):
    #print(h.ndim)
    #print(h.ndim == 3)
    #print(h.shape)
    #print(h.shape[1:])
    #print(h.shape[1:] == (3,3) )
    assert h.ndim == 3, "not 3d"
    #assert(h.shape[1] == 3)
    #assert(h.shape[2] == 3)
    assert h.shape[1:] == (3, 3), "not :x3x3"
    assert not np.any( np.isnan(h.ravel()) )
    assert not np.any( np.isinf(h.ravel()) )


def check_vector4(p):
    assert not issubclass(p.dtype.type, np.integer)
    assert p.shape == (4,), "Vector must be a numpy array of (4) elements"
    assert p[3] == 1.0, "4th element of every Vector must be 1.0"
    assert not np.any( np.isnan(p.ravel()) )
    assert not np.any( np.isinf(p.ravel()) )


def check_vector4_vectorized(pa):
    #assert not issubclass(np.dtype('int8').type, np.integer)
    assert not issubclass(pa.dtype.type, np.integer)
    assert pa.ndim == 2
    assert pa.shape[1:] == (4,), "Vector must be a numpy array of (Nx4) elements"
    e = np.sum(np.abs(pa[:, 3]-1))
    if e > 0.0:
        print("EERROR:", e)
    assert np.allclose(pa[:, 3], 1, 0.00000000000001), "4th element of every Vector must be 1.0"
    assert not np.any( np.isnan(pa.ravel()) )
    assert not np.any( np.isinf(pa.ravel()) )

def check_scalar_vectorized(va, N=None):
    #assert va.ndim == 2 #dont force 2 dim. can accesp .shape==(100,)
    assert not issubclass(va.dtype.type, np.integer)
    n = va.shape[0]
    if va.ndim == 2:
        assert va.shape[1] == 1
        assert va.shape == (n, 1), "values must be a numpy array of (N,) or (Nx1) elements"
    if not N is None:
        assert va.shape[0] == N
    assert not np.any( np.isnan(va.ravel()) )
    assert not np.any( np.isinf(va.ravel()) )

def make_vector4_numpy(v):
    assert issubclass(type(v), np.ndarray)
    assert v.size == 3 or v.size == 4

    assert not np.any( np.isnan(v.ravel()) )
    assert not np.any( np.isinf(v.ravel()) )

    v = v.ravel()
    return np.array((float(v[0]), float(v[1]), float(v[2]), 1.0))


def make_vector4(x, y, z):
    xyz = np.array([x,y,z])
    assert not np.any( np.isnan(xyz) )
    assert not np.any( np.isinf(xyz) )

    if issubclass(type(x), np.ndarray):
        return np.array((float(x[0]), float(y[1]), float(z[1]), 1.0))

    return np.array((float(x), float(y), float(z), 1.0))

def make_vector4_vectorized(x, y, z):
    return make_vector4(x, y, z).reshape((1,4))


def almost_equal4(a, b, TOLERANCE):
    assert not np.any( np.isnan(a.ravel()) )
    assert not np.any( np.isinf(b.ravel()) )
    assert not issubclass(a.dtype.type, np.integer)
    check_vector4(a)
    check_vector4(b)
    return np.sum(np.abs(a - b)) < TOLERANCE


def almost_equal1(a, b, TOLERANCE):
    """ Scalar version """
    assert not issubclass(a.dtype.type, np.integer)
    np.isscalar(a)
    np.isscalar(b)
    return np.sum(np.abs(a - b)) < TOLERANCE


def almost_equal4_vectorized(a, b, TOLERANCE):
    assert a.ndim == 2
    assert b.ndim == 2
    check_vector4_vectorized(a)
    check_vector4_vectorized(a)
    return np.sum(np.abs(a[:, :] - b[:, :])) < TOLERANCE


def check_matrix3(m):
    #print(m.shape)
    assert m.shape == (3, 3)


def make_random_vector(norm, POW, type="rand"):
    if type == "rand":
        r = np.random.rand(3)*2 - 1
    elif type == "randn":
        r = np.random.randn(3)
    else:
        raise Error("nknown random distribution")

    r[:] = np.sign(r[:]) * np.abs(r[:]) ** POW
    r = r / np.sqrt(np.dot(r, r))
    assert (r[0]*r[0] + r[1]*r[1] + r[2]*r[2] - 1) < 0.00000000001
    r = r * norm
    for i in range(0, 3):
        if np.abs(r[i]) < 0.0000001:
            r[i] = 0
    return np.array((r[0], r[1], r[2], 1))


def make_random_vector_vectorized(N, norm, POW, type="rand", normalize=True):
    r = np.ones((N, 4))

    #not tested:
    if type == "rand":
        r[:, 0:3] = np.random.rand(N, 3)*2 - 1
    elif type == "randn":
        r[:, 0:3] = np.random.randn(N, 3)
    else:
        raise Error("nknown random distribution")

    r[:, 0:3] = np.sign(r[:, 0:3]) * np.abs(r[:, 0:3]) ** POW
    #r[:,0:3] = r[0:3] / np.tile( np.sqrt( np.sum(r[:,0:3] * r[:,0:3], axis=1, keepdims=True) ) , (1,3) )
    n3 = np.sqrt(np.sum(r[:, 0:3] * r[:, 0:3], axis=1, keepdims=True))
    if normalize:
        r[:, 0:3] = r[:, 0:3] / np.tile(n3, (1, 3))
        s_1 = np.sum(r[:, 0:3] * r[:, 0:3], axis=1) - 1
        assert np.all(np.abs(s_1) < 0.00000000001)
    r[:, 0:3] = r[:, 0:3] * norm
    r[:, 3] = 1
    check_vector4_vectorized(r)
    return r


def normalize_vector(v, snapToZero=False):
    assert not np.any( np.isnan(v.ravel()) )
    assert not np.any( np.isinf(v.ravel()) )

    r = v.copy()
    #r[:] = np.sign(r[:]) * np.abs(r[:]) ** POW
    r = r / np.sqrt(np.dot(r, r))
    assert (r[0]*r[0] + r[1]*r[1] + r[2]*r[2] - 1) < 0.00000000001
    if snapToZero:
        for i in range(0, 3):
            if np.abs(r[i]) < 0.0000001:
                r[i] = 0
    r[3] = 1
    return r

#todo: http://floating-point-gui.de/errors/comparison/

#todo: write tests for this
def normalize_vector4_vectorized(v, zero_normal="leave_zero_norms"):
    """ returns vectors of either length 1 or zero. """
    N = v.shape[0]
    assert not issubclass(v.dtype.type, np.integer)
    assert not np.any( np.isnan(v) )
    assert not np.any( np.isinf(v) )

    # norms = np.linalg.norm(v[:,0:3], axis = 1, keepdims=True, ord=2)
    norms = np.sqrt(np.sum(v[:,0:3] * v[:,0:3], axis=1, keepdims=True))
    denominator = np.tile(norms, (1, 4))
    if zero_normal=="leave_zero_norms":
        zeros_i = np.abs(norms.ravel()) < 0.00000001
        non_zero_i = np.logical_not(zeros_i)
        if not np.any(zeros_i):
            c = 1.0 / denominator
        else:
            c = np.ones(denominator.shape)
            c[non_zero_i,:] = 1.0 / denominator[non_zero_i,:]
    else:
        pass
    assert not np.any( np.isnan(c) )
    assert not np.any( np.isinf(c) )
    r = v * c
    assert r.shape[0] == N
    df = np.sum(r[:,0:3] * r[:,0:3], axis=1)
    #print(df.shape)
    #print(df)
    #print(non_zero_i)
    e1a = np.all(np.abs(df[non_zero_i]-1.0) < 0.00000000001)
    e0a = np.all(np.abs(df[zeros_i]) < 0.00000000001)

    if not (e1a and e0a):
        print("r:", r)
        print("v:", v)
        print("c:", v)
        print("denom: ", denominator)
        print(norms)
        print(denominator)
        print(np.sum(r[:,0:3] * r[:,0:3], axis=1))

    #print (e1a)
    #print (e0a)
    assert e1a and e0a  # np.all(np.logical_or(e1a, e0a))
    r[:, 3] = 1
    return r


def repeat_vect4(N, v4):
    check_vector4(v4)
    _x = v4
    xa = np.tile(np.expand_dims(_x, axis=0), (N, 1))
    assert xa.shape[0] == N
    return xa


import sys
def is_python3():
    #import sys
    v = sys.version_info.major
    return v == 3
