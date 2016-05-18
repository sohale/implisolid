import numpy as np
from basic_functions import check_scalar_vectorized, check_vector3_vectorized

from implicit import ImplicitFunction


class SimpleBlend (ImplicitFunction):
    """ A simplistic blend using the exp() function """

    def __init__(self, a, b, afactor=1.0, bfactor=1.0):
        assert isinstance(a, ImplicitFunction)
        assert isinstance(b, ImplicitFunction)
        self.a = a
        self.b = b
        self.afactor = afactor
        self.bfactor = bfactor
        assert self.afactor > 0
        assert self.bfactor > 0

    def implicitFunction(self, p):
        check_vector3_vectorized(p)
        va = self.a.implicitFunction(p)
        vb = self.b.implicitFunction(p)

        fa = self.afactor
        fb = self.bfactor
        v = -(1 - (fa*np.exp(va) + fb*np.exp(vb)) / (fa+fb))
        check_scalar_vectorized(v)
        return v

    def implicitGradient(self, p):
        check_vector3_vectorized(p)
        va = self.a.implicitFunction(p)
        vb = self.b.implicitFunction(p)
        ca = np.tile(va[:, np.newaxis], (1, 3))
        cb = np.tile(vb[:, np.newaxis], (1, 3))

        grada = self.a.implicitGradient(p)
        gradb = self.b.implicitGradient(p)
        # not tested
        fa = self.afactor
        fb = self.bfactor
        grad = + (fa*np.exp(ca)*grada + fb*np.exp(cb)*gradb) / (fa+fb)
        check_vector3_vectorized(grad)
        return grad


__all__ = ['SimpleBlend']
