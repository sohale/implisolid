import numpy as np
from basic_functions import check_scalar_vectorized
from basic_functions import check_vector3_vectorized

from implicit import ImplicitFunction


class CrispSubtract(ImplicitFunction):
    def __init__(self, a, b):
        assert isinstance(a, ImplicitFunction)
        assert isinstance(b, ImplicitFunction)
        self.a = a
        self.b = b

    def implicitFunction(self, p):
        check_vector3_vectorized(p)
        va = self.a.implicitFunction(p)
        vb = - self.b.implicitFunction(p)

        c = 1 - np.greater(va, vb)
        v = va * c + vb * (1-c)
        check_scalar_vectorized(v)
        return v

    def implicitGradient(self, p):
        check_vector3_vectorized(p)
        va = self.a.implicitFunction(p)
        vb = - self.b.implicitFunction(p)

        c = 1 - np.greater(va, vb)
        c = np.tile(np.expand_dims(c, axis=1), (1, 3))
        assert c.shape[1:] == (3,)

        grada = self.a.implicitGradient(p)
        gradb = -self.b.implicitGradient(p)

        grad = grada * c + gradb * (1-c)
        check_vector3_vectorized(grad)
        return grad


class CrispUnion(ImplicitFunction):
    def __init__(self, a, b):
        assert isinstance(a, ImplicitFunction)
        assert isinstance(b, ImplicitFunction)
        self.a = a
        self.b = b

    def implicitFunction(self, pa):
        check_vector3_vectorized(pa)

        va = self.a.implicitFunction(pa)
        vb = self.b.implicitFunction(pa)
        c = np.greater(va, vb)
        v = va * c + vb * (1-c)

        check_scalar_vectorized(v)

        return v

    def implicitGradient(self, p):
        check_vector3_vectorized(p)
        va = self.a.implicitFunction(p)
        vb = self.b.implicitFunction(p)

        c = np.greater(va, vb)
        c = np.tile(np.expand_dims(c, axis=1), (1, 3))
        assert c.shape[1:] == (3,)
        grada = self.a.implicitGradient(p)
        gradb = self.b.implicitGradient(p)
        grad = grada * c + gradb * (1-c)
        check_vector3_vectorized(grad)

        return grad


class CrispIntersection(ImplicitFunction):
    """ """
    def __init__(self, a, b):
        assert isinstance(a, ImplicitFunction)
        assert isinstance(b, ImplicitFunction)
        self.a = a
        self.b = b

    def implicitFunction(self, p):
        check_vector3_vectorized(p)
        va = self.a.implicitFunction(p)
        vb = self.b.implicitFunction(p)
        c = 1.0-np.greater(va, vb)
        v = va * c + vb * (1-c)
        check_scalar_vectorized(v)
        return v

    def implicitGradient(self, p):
        check_vector3_vectorized(p)
        va = self.a.implicitFunction(p)
        vb = self.b.implicitFunction(p)

        c = 1 - np.greater(va, vb)

        c = np.tile(c[:, np.newaxis], (1, 3))
        assert c.shape[1:] == (3,)

        grada = self.a.implicitGradient(p)
        gradb = self.b.implicitGradient(p)
        grad = grada * c + gradb * (1-c)
        check_vector3_vectorized(grad)
        return grad


    # todo: non-crisp np.greater, i.e. smooth transition min.

__all__ = ['CrispSubtract', 'CrispIntersection', 'CrispUnion']
