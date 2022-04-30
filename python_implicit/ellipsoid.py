import numpy as np

from primitives import ImplicitFunctionPointwise
from primitives import UnitSphere
from basic_types import make_inverse, check_vector4, check_matrix4


class Ellipsoid(ImplicitFunctionPointwise):

    def __init__(self, m):
        check_matrix4(m)
        self.matrix = m.copy()
        self.invmatrix = make_inverse(m)
        self.sphere = UnitSphere()
        assert isinstance(self.sphere, ImplicitFunctionPointwise)

    def implicitFunction(self, p):
        check_vector4(p)
        tp = np.dot(self.invmatrix, p)
        v = self.sphere.implicitFunction(tp)
        return v

    def implicitGradient(self, p):  # -> Vector3D :
        check_vector4(p)
        tp = np.dot(self.invmatrix, p)
        g = self.sphere.implicitGradient(tp)
        g[3] = 0  # important
        v4 = np.dot(np.transpose(self.invmatrix), g)
        v4[3] = 1
        return v4

    def hessianMatrix(self, p):
        #warning: not tested
        check_vector4(p)
        tp = np.dot(self.invmatrix, p)
        h1 = self.sphere.hessianMatrix(tp)
        h = np.dot(h1, self.invmatrix)  # which one is correct?
        h = np.dot(self.invmatrix, h1)   # which one is correct?
        raise VirtualException()
        return h


__all__ = ['Ellipsoid']
