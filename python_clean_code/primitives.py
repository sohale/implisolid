import numpy as np

from basic_functions import check_vector3, make_vector3
from implicit_config import INTEGRITY_TOLERANCES_NORM


class ImplicitFunctionPointwise(object):

    def implicitFunction(self, p):
        raise Exception()

    def implicitGradient(self, p):
        raise Exception()

    def hessianMatrix(self, p):
        raise Exception()

    def integrity_invariant(self):
        return False


class MultiImplicitObject(object):

    def multiImplicitFunction(self, p):
        raise Exception()

    def multiImplicitGradients(self, p):
        raise Exception()

    def multiImplicitHessians(self, p):
        raise Exception()

    def integrity_invariant(self):
        return False


class SignedDistanceImplicitPointwise(object):
    """ Interface only (can be used in multiple inheritance) """
    pass


class PrimitiveBase(ImplicitFunctionPointwise):
    pass


class UnitSphere(PrimitiveBase):

    def implicitFunction(self, p):
        check_vector3(p)
        return 1.0 - (np.dot(p, p))

    def implicitGradient(self, p):
        check_vector3(p)
        grad = -2*p
        check_vector3(grad)
        return grad

    def integrity_invariant(self):
        return True


def norm2(v):
    return v[0]*v[0]+v[1]*v[1]+v[2]*v[2]


class UnitCube1(PrimitiveBase, SignedDistanceImplicitPointwise):
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

        # assert norm2(self.n0[-1]) - 1 == 0.0

        side(1, 0, 0)
        side(-1, 0, 0)
        side(0, 1, 0)
        side(0, -1, 0)
        side(0, 0, 1)
        side(0, 0, -1)

    def integrity_invariant(self):
        integrity = True
        for i in range(6):
            integrity = integrity and np.abs(norm2(self.n0[i]) - 1) < INTEGRITY_TOLERANCES_NORM

        return integrity

    def implicitFunction(self, p):
        check_vector3(p)
        chosen_i = None
        v = +np.infty
        for i in range(len(self.p0)):
            p0 = self.p0[i]
            n0 = self.n0[i]
            vi = np.dot(p-p0, n0)
            if vi < v:
                v = vi
                chosen_i = i

        assert chosen_i is not None

        return v

    def implicitGradient(self, p):
        check_vector3(p)
        chosen_i = None
        v = +np.infty
        for i in range(len(self.p0)):
            p0 = self.p0[i]
            n0 = self.n0[i]
            vi = np.dot(p-p0, n0)
            if vi < v:
                v = vi
                chosen_i = i
                grad = n0

        assert chosen_i is not None
        check_vector3(grad)
        return grad


__all__ = ['UnitSphere', 'UnitCube1', 'ImplicitFunctionPointwise']

if __name__ == "__main__":
    pass
