import numpy as np

from basic_functions import make_inverse, is_python3
import implicit
from basic_functions import check_matrix4
from basic_functions import check_scalar_vectorized
from basic_functions import check_vector3_vectorized
from implicit import ImplicitFunction
from transformed import Transformable


class Ellipsoid(ImplicitFunction):

    def __init__(self, m):
        assert m.shape == (4, 4)
        np.allclose(m[3, :], np.array((0, 0, 0, 1)), atol=0.001)

        self.matrix = m.copy()
        self.invmatrix = make_inverse(m)
        self.sphere = implicit.UnitSphere()
        assert isinstance(self.sphere, implicit.ImplicitFunction)

    def implicitFunction(self, pa):
        check_vector3_vectorized(pa)
        pa = np.concatenate((pa, np.ones((pa.shape[0], 1))), axis=1)
        tp = np.dot(self.invmatrix, np.transpose(pa))  # inefficient. todo: multiply from right => will be efficient
        tp = np.transpose(tp)
        tp = tp[:, :3]
    #    print pa, tp
        v = self.sphere.implicitFunction(tp)
        check_scalar_vectorized(v)
#        print v
        return v

    def implicitGradient(self, pa):  # -> Vector3D :
        check_vector3_vectorized(pa)
        pa = np.concatenate((pa, np.ones((pa.shape[0], 1))), axis=1)
        tp = np.dot(self.invmatrix, np.transpose(pa))
        tp = np.transpose(tp)  # inefficient
        tp = tp[:, :3]
        g = self.sphere.implicitGradient(tp)
        check_vector3_vectorized(g)  # not needed

        g = np.concatenate((g, np.ones((g.shape[0], 1))), axis=1)
        v4 = np.dot(np.transpose(self.invmatrix), np.transpose(g))
        v4 = np.transpose(v4)  # not efficient
        v3 = v4[:, :3]
        check_vector3_vectorized(v3)
        return v3

    def hessianMatrix(self, p):
        check_vector3_vectorized(p)
        raise VirtualException()

    def __str__(self):
        return "Ellipsoid(vectorized)" + str(self.matrix)  # .split().join(";")


class Transformed(ImplicitFunction, Transformable):
    """ See the @vectorized.Ellipsoid and @nonvec.Transformed classes."""
    def __init__(self, base_object, m=None):
        if is_python3():
            super().__init__(initialMatrix=m)
        else:
            # super(self.__class__, self).__init__(initialMatrix=m)
            super(Transformed, self).__init__(initialMatrix=m)
            # print(self.__class__)

        # print(type(base_object))
        assert issubclass(type(base_object), ImplicitFunction)
        self.base_object = base_object

        if m is None:
            m = np.eye(4)
        check_matrix4(m)
        self.matrix = m.copy()
        self.invmatrix = make_inverse(m)

        assert isinstance(self.base_object, ImplicitFunction)

    def implicitFunction(self, p):
        check_vector3_vectorized(p)
        p = np.concatenate((p, np.ones((p.shape[0], 1))), axis=1)
        tp = np.dot(self.invmatrix, np.transpose(p))
        tp = np.transpose(tp)
        tp = tp[:, :3]
        v = self.base_object.implicitFunction(tp)
        check_scalar_vectorized(v)
        return v

    def implicitGradient(self, p):  # -> Vector3D :
        check_vector3_vectorized(p)
        p = np.concatenate((p, np.ones((p.shape[0], 1))), axis=1)
        tp = np.dot(self.invmatrix, np.transpose(p))
        tp = np.transpose(tp)
        tp = tp[:, :3]
        g = self.base_object.implicitGradient(tp)
        check_vector3_vectorized(g)
        # g[:, 3] = 0  # important

        g = np.concatenate((g, np.ones((g.shape[0], 1))), axis=1)

        v4 = np.dot(np.transpose(self.invmatrix), np.transpose(g))
        v4 = np.transpose(v4)  # not efficient
        v3 = v4[:, :3]
        check_vector3_vectorized(v3)
        return v3

    def hessianMatrix(self, p):
        check_vector3_vectorized(p)
        p = np.concatenate((p, np.ones((p.shape[0], 1))), axis=1)
        tp = np.dot(self.invmatrix, np.transpose(p))
        tp = np.transpose(tp)

        h1 = self.base_object.hessianMatrix(tp)
        h = np.dot(h1, self.invmatrix)  # which one is correct?
        h = np.dot(self.invmatrix, np.tanspose(h1))   # which one is correct?
        raise VirtualException()
        return h


__all__ = ['Ellipsoid', 'Transformed']
