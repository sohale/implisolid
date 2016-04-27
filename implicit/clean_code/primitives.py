import numpy as np

from basic_functions import check_vector4, check_matrix3


class ImplicitFunctionPointwise(object):

    def implicitFunction(self, p):
        raise VirtualException()

    def implicitGradient(self, p):
        raise VirtualException()

    def hessianMatrix(self, p):
        raise VirtualException()

    def integrity_invariant(self):
        return False


class MultiImplicitObject(object):

    def multiImplicitFunction(self, p):
        #return [0]
        raise VirtualException()

    def multiImplicitGradients(self, p):
        raise VirtualException()

    def multiImplicitHessians(self, p):
        raise VirtualException()

    def integrity_invariant(self):
        return False


class SignedDistanceImplicitPointwise(object):
    """ Interface only (can be used in multiple inheritance) """
    pass


class PrimitiveBase(ImplicitFunctionPointwise):
    pass

#UnitSphere is not a SignedDistanceImplicitPointwise

class UnitSphere(PrimitiveBase):

    def implicitFunction(self, p):
        check_vector4(p)
        return 1.0 - (np.dot(p[:3], p[:3]))

    def implicitGradient(self, p):
        check_vector4(p)
        grad = -2*p
        grad[3] = 1
        check_vector4(grad)
        return grad

    def hessianMatrix(self, p):
        check_vector4(p)
        h = np.array([[-2, 0, 0],  [0, -2, 0],  [0, 0, -2]], ndmin=2)
        check_matrix3(h)
        return h

    def integrity_invariant(self):
        return True

from basic_functions import make_vector4
from implicit_config import INTEGRITY_TOLERANCES_NORM

class UnitCube1(PrimitiveBase, SignedDistanceImplicitPointwise):
    def __init__(self, size=1.0):
        self.p0 = []
        self.n0 = []

        def side(x, y, z):
            p0 = (make_vector4(x, y, z) + 0.0)
            p0 = p0 / 2.0 * size
            # n0 points inwards
            n0 = -make_vector4(x, y, z)
            n0 = n0
            self.p0 += [p0]
            self.n0 += [n0]
            self.p0[-1][3] = 1
            self.n0[-1][3] = 1

            def norm2(v):
                return v[0]*v[0]+v[1]*v[1]+v[2]*v[2]

            #print(norm2(self.n0[-1][0:3]))
            assert norm2(self.n0[-1][0:3]) - 1 == 0.0

        #for x in range( -1,2 ):
        #    for y in range( -1,2 ):
        #        for z in range( -1,2 ):
        #           side(x,y,z)
        side(1, 0, 0)
        side(-1, 0, 0)
        side(0, 1, 0)
        side(0, -1, 0)
        side(0, 0, 1)
        side(0, 0, -1)
        #print ( self.p0 )
        #print ( self.n0 )

    def integrity_invariant(self):
        integrity = True
        for i in range(6):
            integrity = integrity and np.abs(norm2(self.n0[i][0:3]) - 1) < INTEGRITY_TOLERANCES_NORM
        #todo: Check convexity
        return integrity

    #todo: evaluate gradient and implicit function at the same time
    def implicitFunction(self, p):
        check_vector4(p)
        chosen_i = None
        v = +np.infty
        for i in range(len(self.p0)):
            p0 = self.p0[i]
            n0 = self.n0[i]
            vi = np.dot(p-p0, n0)
            if vi < v:
                v = vi
                chosen_i = i

        assert not chosen_i is None

        return v

    def implicitGradient(self, p):
        check_vector4(p)
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

        assert not chosen_i is None

        #not tested
        grad[3] = 1
        check_vector4(grad)
        return grad

    def hessianMatrix(self, p):
        check_vector4(p)
        h = np.array([[0, 0, 0],  [0, 0, 0],  [0, 0, 0]], ndmin=2)
        check_matrix3(h)
        return h


__all__ = ['UnitSphere', 'UnitCube1', 'ImplicitFunctionPointwise']

if __name__ == "__main__":
    pass
