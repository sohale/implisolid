import numpy as np
from implicit_vectorized import ImplicitFunctionVectorized
from basic_types import check_vector4_vectorized, make_vector4, check_vector4
#from numerical_utils import numerical_gradient
from implicit_config import config

# todo: class CutCone


class SimpleCylinder(ImplicitFunctionVectorized):

    def __init__(self, A, w, u, radius_u, radius_v, c_len):
        """ The cross section of a SimpleCylinder is always a circle.
        """
        (self.A, self.w, self.u, self.radius_u, self.radius_v, self.c_len) = \
            (A, w, u, radius_u, radius_v, c_len)

        check_vector4(w)
        check_vector4(u)
        check_vector4(A)
        assert w[3] == 1
        assert u[3] == 1
        assert A[3] == 1
        w = w[:3]
        u = u[:3]
        A = A[:3]

        v = np.cross(u[:3], w[:3])
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

        assert np.isscalar(radius_u)
        assert np.isscalar(radius_v)

        norm_tol = 0.00000001  # 1000 km
        assert np.abs(np.linalg.norm(w)-1.0) < norm_tol
        assert np.abs(np.linalg.norm(u)-1.0) < norm_tol
        self.UVW = np.concatenate( (self.u, self.v, self.w), axis=1 )
        self.UVW_inv = np.linalg.inv( self.UVW )
        assert self.c_len > 0
        assert self.integrity_invariant()

    def integrity_invariant(self):
        #config = threeD_printing_config_profile
        norm_tol = 0.00000001  # 1000 km
        matrix_inv_tol = 0.000001
        #sane = True
        sd = {"sane": True}
        def check(boolean, reason):
            sd["sane"] = sd["sane"] and boolean
            if not boolean:
                print("Error:", reason)
            pass
        check(np.abs(np.linalg.norm(self.w)-1.0) < norm_tol, "w1")
        check(np.abs(np.linalg.norm(self.u)-1.0) < norm_tol, "u1")
        check(np.abs(np.linalg.norm(self.v)-1.0) < norm_tol, "v1")
        check(self.UVW.shape == (3, 3), "33")
        check(np.allclose(np.dot(self.UVW, self.UVW_inv), np.eye(3), matrix_inv_tol), "MMinv")
        check(np.allclose(np.dot(self.UVW_inv, self.UVW), np.eye(3), matrix_inv_tol), "MinvM")
        check(self.c_len > config.numerical_min_length, "minlen")  # functions, not actual models
        check(self.c_len > 0, "clen")
        check(self.radius_u > config.numerical_min_length, "r0u")
        check(self.radius_v > config.numerical_min_length, "r0v")
        check(self.radius_u == self.radius_v, "rr")
        check(np.linalg.norm(self.v.ravel() - np.cross(self.u.ravel(), self.w.ravel())) < norm_tol, "x")
        return sd["sane"]

    def implicitFunction(self, xa, return_grad=False):
        check_vector4_vectorized(xa)
        assert self.w.shape == (3, 1)
        assert self.A.shape == (3, 1)
        x = xa[:, 0:3]
        count = x.shape[0]
        assert x.shape == (count, 3)
        aa = np.tile(np.transpose(self.A), (count, 1))

        assert x.shape == (count, 3)
        assert aa.shape == (count, 3)
        t_ = np.dot(x - aa, self.w)  # Nx1
        assert t_.shape == (count,1)
        t = t_[:, 0]  # N,
        assert t.shape == (count,)
        t_arr_1xN = t.reshape((1, count))
        assert t_arr_1xN.shape == (1, count)
        p = np.transpose(aa) + np.dot(self.w, t_arr_1xN)   # 3x1 * 1xcount
        assert p.shape == (3, count)
        assert x.shape == (count, 3)
        # alpha_beta = np.dot(self.UVW_inv, np.transpose(x) - p)  # alpha,beta  #not correct
        # assert alpha_beta.shape == (3, count)
        # theta = np.arctan2(alpha_beta[1, :], alpha_beta[0, :])
        assert self.radius_u == self.radius_v

        r = np.linalg.norm(x - np.transpose(p), ord=2, axis=1)
        assert r.shape == (count,)

        # todo: the following is incorrect

        t0 = t
        t1 = self.c_len - t
        r_ = self.radius_u - r
        # t0,t1,r_  as if  x,y,z

        print t0.shape, t1.shape, r_.shape
        m3 = np.concatenate((t0[:, np.newaxis], t1[:, np.newaxis], r_[:, np.newaxis]), axis=1)
        print m3.shape
        fval = np.min(m3, axis=1)

        if not return_grad:
            return fval 
        else:
            # c_t01 = (t - self.c_len/2.0)/(self.c_len/2.0)   # [-1,+1]
            # c_t = (1.0 + np.sign(c_t01))/2.0  # {0,1}
            # t01_grad = self.w * (1-c_t) + (-self.w) * c_t  # toward center
            # c_
            c_t0 = np.logical_and(t0 <= t1,  t0 <= r_)
            c_t1 = np.logical_and(t1 <= t0,  t1 <= r_)
            c_r = np. logical_and(r_ <= t0,  r_ <= t1)
            assert c_t0.ndim == 1
            assert c_t1.ndim == 1
            c_t0 = np.tile(c_t0[:,np.newaxis], (1,3) )
            c_t1 = np.tile(c_t1[:,np.newaxis], (1,3) )
            c_r = np.tile(c_r[:,np.newaxis], (1,3) )
            grad_t0 = np.tile(self.w[np.newaxis,:,0], (count,1))
            grad_t1 = np.tile(-self.w[np.newaxis,:,0], (count,1))
            grad_r = np.transpose(p) - x
            print c_t0.shape, "c_t0"
            print c_t1.shape, "c_t1"
            print c_r.shape, "c_r"
            print grad_r.shape, "g_r"
            print grad_t0.shape, "g_t0"
            print grad_t1.shape, "g_t1"
            #print t0.shape, "t0"
            #print t1.shape, "t1"
            #return (c_t0) * grad_t0 + (c_t1) * grad_t1 + (c_r) * grad_r
            a=  (c_t0) * grad_t0 
            b= (c_t1) * grad_t1 
            c=  (c_r) * grad_r
            g = a+b+c

            g4 = np.concatenate( (g, np.ones((count, 1), dtype=float)), axis=1)
            check_vector4_vectorized(g4)
            #in progress
            return fval, g4

    def implicitGradient(self, x):
        # return numeric_utils.generic_slow_gradient(self, x)
        f, g = self.implicitFunction(x, return_grad=True)
        print("ggggg", g.shape)
        return g

    def curvature(self, x):
        check_vect2(x)
        raise NotImplementedError()


__all__ = ['SimpleCylinder']
