import numpy as np

from basic_functions import check_vector3_vectorized
from basic_functions import check_vector3, make_vector3


class ImplicitFunction(object):
    """ Functions in this type receive numpy vectors of size Nx4 """

    def implicitFunction(self, pv):
        raise Exception()

    def implicitGradient(self, pv):
        """ Returns a vector of size N x 3 where N is the number of points.  """
        raise Exception()

    def hessianMatrix(self, pv):
        """ Returns a vector of size N x 3 x 3 where N is the number of points"""
        raise Exception()

    def integrity_invariant(self):
        return False


class SignedDistanceImplicit(object):
    """ Interface only (can be used in multiple inheritance) """
    pass


class UnitSphere(ImplicitFunction):
    def implicitFunction(self, pv):

        return 1.0 - np.sum(pv * pv, axis=1)

    def implicitGradient(self, pv):
        assert pv.ndim == 2
        assert pv.shape[1:] == (3,)
        grad = -2*pv

        return grad


class UnitCube1(ImplicitFunction, SignedDistanceImplicit):
    def __init__(self, size=1.0):
        self.p0 = []
        self.n0 = []

        def side(x, y, z):
            p0 = make_vector3(x, y, z)
            p0 = p0 / 2.0 * size
            n0 = -make_vector3(x, y, z)
            n0 = n0
            self.p0 += [p0]
            self.n0 += [n0]
            # print(self.p0[-1])
            check_vector3(self.p0[-1])
            check_vector3(self.n0[-1])

            def norm2(v):
                return v[0]*v[0]+v[1]*v[1]+v[2]*v[2]
            assert norm2(self.n0[-1]) - 1 == 0.0

        side(1, 0, 0)
        side(-1, 0, 0)
        side(0, 1, 0)
        side(0, -1, 0)
        side(0, 0, 1)
        side(0, 0, -1)

    def implicitFunction(self, p):
        check_vector3_vectorized(p)

        sides = len(self.p0)
        n = p.shape[0]
        temp = np.zeros((n, sides))
        for i in range(sides):
            p0 = self.p0[i]
            n0 = self.n0[i]
            sub = p - np.tile(p0[np.newaxis, :], (n, 1))
            vi = np.dot(sub, n0)
            temp[:, i] = vi
        va = np.amin(temp, axis=1)
        return va

    def implicitGradient(self, p):
        check_vector3_vectorized(p)

        sides = 6
        na = np.zeros((sides, 3))
        n = p.shape[0]
        temp = np.zeros((n, sides))
        for i in range(len(self.p0)):
            p0 = self.p0[i]
            n0 = self.n0[i]
            sub = p - np.tile(p0[np.newaxis, :], (n, 1))
            vi = np.dot(sub, n0)

            temp[:, i] = vi

            na[i, :] = n0

        ia = np.argmin(temp, axis=1)

        assert ia.shape == (n,)

        g = na[ia, :]

        check_vector3_vectorized(g)
        return g


__all__ = ['ImplicitFunction', 'UnitSphere', 'UnitCube1']

if __name__ == "__main__":
    pass
