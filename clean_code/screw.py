import numpy as np
from implicit import ImplicitFunction
from basic_functions import check_vector3_vectorized, check_vector3, repeat_vect3, check_vector2, make_vector3

import vector3
from implicit_config import config


def numerical_gradient(iobj, pos0, delta_t=0.01/10.0/10.0, order=5):

    check_vector3(pos0)
    assert issubclass(type(iobj), vector3.ImplicitFunction)
    m = order

    _VERBOSE = False
    import finite_diff_weights

    sample_points = range(-m, m+1)
    n = m*2+1

    x0 = 0
    findiff_weights = finite_diff_weights.weights(k=1, x0=x0, xs=np.array(sample_points) * delta_t)

    pos0_3 = repeat_vect3(1, pos0)
    pos = np.tile(pos0_3, (3*n, 1))
    assert not issubclass(pos.dtype.type, np.integer)

    dx = repeat_vect3(1, make_vector3(1, 0, 0))
    dy = repeat_vect3(1, make_vector3(0, 1, 0))
    dz = repeat_vect3(1, make_vector3(0, 0, 1))
    dxyz = [dx, dy, dz]

    ci = 0
    for d in range(3):
        for i in sample_points:
            dd = dxyz[d]
            pos[ci, :] = pos[ci, :] + (dd * delta_t * float(i))
            ci += 1

    v = iobj.implicitFunction(pos)
    v3 = np.reshape(v, (3, n), order='C')

    Lipchitz_beta = 1

    d0 = np.abs(np.diff(v3, axis=1))
    nonsmooth_ness = d0 / (np.abs(delta_t)**Lipchitz_beta)

    d = np.abs(np.diff(v3, n=1, axis=1)) / np.abs(delta_t)
    d = d - np.tile(np.mean(d, axis=1, keepdims=True), (1, d.shape[1]))
    d = np.abs(d) / np.abs(delta_t)
    d = d - np.tile(np.mean(d, axis=1, keepdims=True), (1, d.shape[1]))

    if(np.max(np.ravel(nonsmooth_ness))) > 100*10:
        print "warning: nonsmooth ",

    """ Calculating the numerical derivative using finite difference (convolution with weights) """

    grad_cnv = np.dot(v3, findiff_weights)

    if np.max(np.abs(grad_cnv)) > 100:
        pass

    if _VERBOSE:
        np.set_printoptions(formatter={'all': lambda x: ''+("%2.19f" % (x,))})

    """ Calculating the numerical derivative using 'mean of diff' """
    grad_mean = np.mean(-np.diff(v3, axis=1) / delta_t, axis=1)
    if _VERBOSE:
        print("grad_mean: ", grad_mean)
        print("grad_convolusion: ", grad_cnv)

    if False:
        g = iobj.implicitGradient(pos0_3)
    if _VERBOSE:
        print("grad_analytical: ", g)

        print("Errors:")
        print("conv error: ", g - grad_cnv)

        print("to be continued")

    return grad_cnv


class Screw(ImplicitFunction):
    def __init__(self, A, w, u, slen, r0, delta, twist_rate, phi0=0.0, phi_func=None):
        """ class name is: vectorized.Screw
            twist_rate:  mm per cycle.
            phi0: cycle (0.5 = half a cycle)

            A bug fixed. But not tested
        """
        (self.A, self.w, self.u, self.slen, self.r0, self.delta, self.twist_rate, self.phi0, self.phi_func) = \
            (A, w, u, slen, r0, delta, twist_rate, phi0, phi_func)
        assert phi_func is None
        v = np.cross(u, w)
        assert w.shape == (3,)
        assert u.shape == (3,)
        assert v.shape == (3,)
        assert A.shape == (3,)

        self.u = u[:, np.newaxis]
        self.v = v[:, np.newaxis]
        self.w = w[:, np.newaxis]
        self.A = A[:, np.newaxis]
        assert self.u.shape == (3, 1)
        assert self.v.shape == (3, 1)
        assert self.w.shape == (3, 1)
        assert self.A.shape == (3, 1)

        assert np.isscalar(r0)
        assert np.isscalar(delta)
        assert np.isscalar(twist_rate)
        assert np.isscalar(phi0)
        norm_tol = 0.00000001  # 1000 km
        assert np.abs(np.linalg.norm(w)-1.0) < norm_tol
        assert np.abs(np.linalg.norm(u)-1.0) < norm_tol
        self.UVW = np.concatenate((self.u, self.v, self.w), axis=1)
        self.UVW_inv = np.linalg.inv(self.UVW)
        assert self.slen > 0
        assert self.integrity_invariant()

    def integrity_invariant(self):
        norm_tol = 0.00000001  # 1000 km
        matrix_inv_tol = 0.000001
        sane = True
        sane = sane and np.abs(np.linalg.norm(self.w)-1.0) < norm_tol
        sane = sane and np.abs(np.linalg.norm(self.u)-1.0) < norm_tol
        sane = sane and np.abs(np.linalg.norm(self.v)-1.0) < norm_tol
        sane = sane and self.UVW.shape == (3, 3)
        sane = sane and np.allclose(np.dot(self.UVW, self.UVW_inv), np.eye(3), matrix_inv_tol)
        sane = sane and np.allclose(np.dot(self.UVW_inv, self.UVW), np.eye(3), matrix_inv_tol)
        sane = sane and self.slen > 0
        sane = sane and self.slen > config.numerical_min_length
        sane = sane and self.r0 > config.numerical_min_length
        sane = sane and np.linalg.norm(self.v[:, 0] - np.cross(self.u[:, 0], self.w[:, 0])) < norm_tol
        return sane

    def implicitFunction(self, xa):
        check_vector3_vectorized(xa)
        assert self.w.shape == (3, 1)
        assert self.A.shape == (3, 1)
        x = xa
        count = x.shape[0]
        assert x.shape == (count, 3)
        aa = np.tile(np.transpose(self.A), (count, 1))

        assert x.shape == (count, 3)
        assert aa.shape == (count, 3)
        t = np.dot(x - aa, self.w.ravel())
        assert t.shape == (count,)
        t1 = np.dot(x - aa, self.w).reshape((1, count))
        assert t1.shape == (1, count)
        p = np.transpose(aa) + np.dot(self.w, t1)   # 3x1 * 1xcount
        assert p.shape == (3, count)
        assert x.shape == (count, 3)
        ab = np.dot(self.UVW_inv, np.transpose(x) - p)
        assert ab.shape == (3, count)
        theta = np.arctan2(ab[1, :], ab[0, :])
        r = np.linalg.norm(x - np.transpose(p), ord=2, axis=1)
        assert r.shape == (count,)

        inside_ness = t / self.slen
        inside_ness = 1 - 2 * np.abs(inside_ness-0.5)

        inside_ness = (inside_ness > 0)*1.0

        pi2 = np.pi*2

        def phi(x):
            return np.abs(2*(x - np.floor(x)) - 1.0)*2.0 - 1.0
        return (-r + self.r0 + self.delta * phi(t / self.twist_rate - theta/pi2 + self.phi0))*inside_ness

    def implicitGradient(self, x):
        check_vector3_vectorized(x)
        count = x.shape[0]
        g = np.zeros((count, 3))
        for i in range(x.shape[0]):
            v = x[i, 0:3]
            g[i, :] = numerical_gradient(self, v)
        return g

    def curvature(self, x):
        check_vector2(x)
        raise NotImplementedError()

__all__ = ['Screw']
