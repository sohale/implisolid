import numpy as np

from basic_types import make_inverse, is_python3
import implicit_vectorized
from basic_types import check_matrix4
from basic_types import check_vector4_vectorized, check_scalar_vectorized
from implicit_vectorized import ImplicitFunctionVectorized
from transformed import Transformable


class Ellipsoid(ImplicitFunctionVectorized):

    def __init__(self, m):
        assert m.shape == (4, 4)
        np.allclose(m[3, :], np.array((0, 0, 0, 1)), atol=0.001)

        self.matrix = m.copy()
        self.invmatrix = make_inverse(m)
        self.sphere = implicit_vectorized.UnitSphere()
        assert isinstance(self.sphere, implicit_vectorized.ImplicitFunctionVectorized)

    def implicitFunction(self, pa):
        check_vector4_vectorized(pa)
        tp = np.dot(self.invmatrix, np.transpose(pa))  # inefficient. todo: multiply from right => will be efficient
        tp = np.transpose(tp)  # inefficient.

        v = self.sphere.implicitFunction(tp)
        check_scalar_vectorized(v)
        return v

    def implicitGradient(self, pa):  # -> Vector3D :
        check_vector4_vectorized(pa)
        tp = np.dot(self.invmatrix,  np.transpose(pa))
        tp = np.transpose(tp)  # inefficient

        g = self.sphere.implicitGradient(tp)
        check_vector4_vectorized(g)  # not needed
        #print("g:  ", g)
        g[:, 3] = 0  # important
        v4 = np.dot(np.transpose(self.invmatrix),  np.transpose(g))
        v4 = np.transpose(v4)  # not efficient
        #print("v4:  ", v4)
        v4[:, 3] = 1
        check_vector4_vectorized(v4)
        return v4

    def hessianMatrix(self, p):
        check_vector4_vectorized(p)
        raise VirtualException()

    def __str__(self):
        return "Ellipsoid(vectorized)" + str(self.matrix)  # .split().join(";")


class Transformed(ImplicitFunctionVectorized, Transformable):
    """ See the @vectorized.Ellipsoid and @nonvec.Transformed classes."""
    def __init__(self, base_object, m=None):
        if is_python3():
            super().__init__(initialMatrix=m)
        else:
            #super(self.__class__, self).__init__(initialMatrix=m)
            super(Transformed, self).__init__(initialMatrix=m)
            #print(self.__class__)

        #print(type(base_object))
        assert issubclass(type(base_object), ImplicitFunctionVectorized)
        self.base_object = base_object

        if m is None:
            m = np.eye(4)
        check_matrix4(m)
        self.matrix = m.copy()
        self.invmatrix = make_inverse(m)

        assert isinstance(self.base_object, ImplicitFunctionVectorized)

    def implicitFunction(self, p):
        check_vector4_vectorized(p)

        tp = np.dot(self.invmatrix, np.transpose(p))
        tp = np.transpose(tp)
        v = self.base_object.implicitFunction(tp)
        check_scalar_vectorized(v)
        return v

    def implicitGradient(self, p):  # -> Vector3D :
        check_vector4_vectorized(p)
        tp = np.dot(self.invmatrix, np.transpose(p))
        tp = np.transpose(tp)
        g = self.base_object.implicitGradient(tp)
        check_vector4_vectorized(g)
        g[:, 3] = 0  # important

        v4 = np.dot(np.transpose(self.invmatrix), np.transpose(g))
        v4 = np.transpose(v4)  # not efficient
        v4[:, 3] = 1
        check_vector4_vectorized(v4)
        return v4

    def hessianMatrix(self, p):
        #warning: not tested
        check_vector4_vectorized(p)
        tp = np.dot(self.invmatrix, np.transpose(p))
        tp = np.transpose(tp)

        h1 = self.base_object.hessianMatrix(tp)
        h = np.dot(h1, self.invmatrix)  # which one is correct?
        h = np.dot(self.invmatrix, np.tanspose(h1))   # which one is correct?
        raise VirtualException()
        return h


__all__ = ['Ellipsoid', 'Transformed']
