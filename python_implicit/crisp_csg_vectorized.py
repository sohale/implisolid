import numpy as np
from basic_types import check_vector4_vectorized, check_scalar_vectorized
from basic_types import check_matrix3_vectorized
#from implicit_config import TOLERANCE
#from primitives import ImplicitFunctionPointwise

#from primitives import UnitSphere
from implicit_function_vectorized import ImplicitFunctionVectorized

class CrispSubtract(ImplicitFunctionVectorized):
    def __init__(self, a, b):
        assert isinstance(a, ImplicitFunctionVectorized)
        assert isinstance(b, ImplicitFunctionVectorized)
        self.a = a
        self.b = b

    def implicitFunction(self, p):
        check_vector4_vectorized(p)
        va = self.a.implicitFunction(p)
        vb = - self.b.implicitFunction(p)

        c = 1 - np.greater(va, vb)
        v = va * c + vb * (1-c)
        check_scalar_vectorized(v)
        return v

    def implicitGradient(self, p):
        check_vector4_vectorized(p)
        va = self.a.implicitFunction(p)
        vb = - self.b.implicitFunction(p)

        c = 1 - np.greater(va, vb)
        c = np.tile(np.expand_dims(c, axis=1), (1, 4))
        assert c.shape[1:] == (4,)

        grada = self.a.implicitGradient(p)
        gradb = -self.b.implicitGradient(p)
        gradb[:, 3] = 1
        grad = grada * c + gradb * (1-c)
        check_vector4_vectorized(grad)
        return grad

    def hessianMatrix(self, p):
        check_vector4_vectorized(p)
        va = self.a.implicitFunction(p)
        vb = - self.b.implicitFunction(p)

        c = 1 - np.greater(va, vb)
        c = np.tile(c[:, np.newaxis, np.newaxis], (1, 4, 4))
        assert c.shape[1:] == (4, 4)

        ha = self.a.hessianMatrix(p)
        hb = -self.b.hessianMatrix(p)
        h = ha * c + hb * (1-c)
        check_matrix3_vectorized(h)
        return h


class CrispUnion(ImplicitFunctionVectorized):
    def __init__(self, a, b):
        assert isinstance(a, ImplicitFunctionVectorized)
        assert isinstance(b, ImplicitFunctionVectorized)
        self.a = a
        self.b = b

    def implicitFunction(self, pa):
        check_vector4_vectorized(pa)
        va = self.a.implicitFunction(pa)
        vb = self.b.implicitFunction(pa)

        c = np.greater(va, vb)
        v = va * c + vb * (1-c)
        check_scalar_vectorized(v)

        return v

    def implicitGradient(self, p):
        check_vector4_vectorized(p)
        va = self.a.implicitFunction(p)
        vb = self.b.implicitFunction(p)

        c = np.greater(va, vb)  # shape: (N,)
        c = np.tile(np.expand_dims(c, axis=1), (1, 4))
        assert c.shape[1:] == (4,)
        grada = self.a.implicitGradient(p)
        gradb = self.b.implicitGradient(p)
        grad = grada * c + gradb * (1-c)
        check_vector4_vectorized(grad)

        return grad

    def hessianMatrix(self, p):
        check_vector4_vectorized(p)
        va = self.a.implicitFunction(p)  # why not used???
        vb = self.b.implicitFunction(p)
        c = np.greater(va, vb)  # shape: (N,)
        c = np.tile(np.expand_dims(c, axis=1), (1, 4))

        ha = self.a.hessianMatrix(p)
        hb = self.b.hessianMatrix(p)
        h = ha * c + hb * (1-c)
        check_matrix3_vectorized(h)

        return h

#not tested


class CrispIntersection(ImplicitFunctionVectorized):
    """ """
    def __init__(self, a, b):
        assert isinstance(a, ImplicitFunctionVectorized)
        assert isinstance(b, ImplicitFunctionVectorized)
        self.a = a
        self.b = b

    def implicitFunction(self, p):
        check_vector4_vectorized(p)
        va = self.a.implicitFunction(p)
        vb = self.b.implicitFunction(p)
        c = 1.0-np.greater(va, vb)
        v = va * c + vb * (1-c)
        check_scalar_vectorized(v)
        return v

    def implicitGradient(self, p):
        check_vector4_vectorized(p)
        va = self.a.implicitFunction(p)
        vb = self.b.implicitFunction(p)

        c = 1 - np.greater(va, vb)
        c = np.tile(c[:, np.newaxis], (1, 4))
        assert c.shape[1:] == (4,)

        grada = self.a.implicitGradient(p)
        gradb = self.b.implicitGradient(p)
        grad = grada * c + gradb * (1-c)

        check_vector4_vectorized(grad)
        return grad

    def hessianMatrix(self, p):
        check_vector4_vectorized(p)
        va = self.a.implicitFunction(p)
        vb = self.b.implicitFunction(p)

        c = 1 - np.greater(va, vb)
        c = np.tile(c[:, np.newaxis, np.newaxis], (1, 4, 4))
        assert c.shape[1:] == (4, 4)

        ha = self.a.hessianMatrix(p)
        hb = self.b.hessianMatrix(p)
        h = ha * c + hb * (1-c)

        check_matrix3_vectorized(h)
        return h

    #todo: non-crisp np.greater, i.e. smooth transition min.

__all__ = ['CrispSubtract', 'CrispIntersection', 'CrispUnion']
