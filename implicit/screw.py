import numpy as np
from implicit_vectorized import ImplicitFunctionVectorized
from basic_types import check_vector4_vectorized
from numerical_utils import numerical_gradient
from implicit_config import config


class Screw(ImplicitFunctionVectorized):
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
        # print(self.slen)
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
        check_vector4_vectorized(xa)
        #print(self.w.shape)
        assert self.w.shape == (3, 1)
        assert self.A.shape == (3, 1)
        x = xa[:, 0:3]
        #print(x.shape)
        count = x.shape[0]
        assert x.shape == (count, 3)
        #x - self.A
        #t = np.dot(np.transpose(self.w), x - self.A)
        aa = np.tile(np.transpose(self.A), (count, 1))

        assert x.shape == (count, 3)
        assert aa.shape == (count, 3)
        t = np.dot(x - aa, self.w.ravel())
        #print(t.shape)
        assert t.shape == (count,)
        t1 = np.dot(x - aa, self.w).reshape((1, count))
        assert t1.shape == (1, count)
        p = np.transpose(aa) + np.dot(self.w, t1)   # 3x1 * 1xcount
        assert p.shape == (3, count)
        assert x.shape == (count, 3)
        ab = np.dot(self.UVW_inv, np.transpose(x) - p)
        assert ab.shape == (3, count)
        theta = np.arctan2(ab[1, :], ab[0, :])
        # assert ab[2,:] == t
        r = np.linalg.norm(x - np.transpose(p), ord=2, axis=1)
        assert r.shape == (count,)

        # inside = t - self.len > 0
        inside_ness = t / self.slen
        # print(self.slen)

        #todo: see Cylinder
        #inside_ness[inside_ness<1] = 1
        #inside_ness[t<0] = 1 + t[t<0]
        #inside_ness = (inside_ness-1)*100+1
        #inside_ness[inside_ness>1] = 1
        inside_ness = 1 - 2 * np.abs(inside_ness-0.5)
        #inside_ness = (inside_ness-1)*100+1
        inside_ness = (inside_ness > 0)*1.0

        pi2 = np.pi*2

        def phi(x):
            return np.abs(2*(x - np.floor(x)) - 1.0)*2.0 - 1.0
        screw_ness = (-r + self.r0 + self.delta * phi(t / self.twist_rate - theta/pi2 + self.phi0))*inside_ness

        lidness0 = t
        lidness1 = self.slen - t

        m3 = np.concatenate((lidness0[:, np.newaxis], lidness1[:, np.newaxis], screw_ness[:, np.newaxis]), axis=1)
        fval = np.min(m3, axis=1)

        #return screw_ness
        return fval

    def implicitGradient(self, x):
        check_vector4_vectorized(x)
        count = x.shape[0]
        g = np.zeros((count, 4))
        for i in range(x.shape[0]):
            v = x[i, 0:4]
            # inefficient: not vectorised
            g[i, :] = numerical_gradient(self, v, is_vectorized=True)
        return g

    def curvature(self, x):
        check_vect2(x)
        raise NotImplementedError()

__all__ = ['Screw']
