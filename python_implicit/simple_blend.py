import numpy as np
from basic_types import check_vector4_vectorized, check_scalar_vectorized
from implicit_config import TOLERANCE
from primitives import ImplicitFunctionPointwise

from primitives import UnitSphere
from implicit_vectorized import ImplicitFunctionVectorized


class SimpleBlend (ImplicitFunctionVectorized):
    """ A simplistic blend using the exp() function """

    def __init__(self, a, b, afactor=1.0, bfactor=1.0):
        assert isinstance(a, ImplicitFunctionVectorized)
        assert isinstance(b, ImplicitFunctionVectorized)
        self.a = a
        self.b = b
        self.afactor = afactor
        self.bfactor = bfactor
        assert self.afactor > 0
        assert self.bfactor > 0

    def implicitFunction(self, p):
        check_vector4_vectorized(p)
        va = self.a.implicitFunction(p)
        vb = self.b.implicitFunction(p)
        #v = va * self.afactor  +  vb * self.bfactor + 1

        fa = self.afactor
        fb = self.bfactor
        v = -(1 - (fa*np.exp(va)   +  fb*np.exp(vb) ) / (fa+fb))
        check_scalar_vectorized(v)
        return v

    def implicitGradient(self, p):
        check_vector4_vectorized(p)
        va = self.a.implicitFunction(p)
        vb = self.b.implicitFunction(p)
        ca = np.tile( va[:,np.newaxis], (1,4) )
        cb = np.tile( vb[:,np.newaxis], (1,4) )

        grada = self.a.implicitGradient(p)
        gradb = self.b.implicitGradient(p)
        #not tested
        fa = self.afactor
        fb = self.bfactor
        grad = + (fa*np.exp(ca)*grada   +  fb*np.exp(cb)*gradb ) / (fa+fb)
        grad[:,3] = 1
        check_vector4_vectorized(grad)
        return grad

    def hessianMatrix(self, p):
        check_vector4_vectorized(p)
        ha = self.a.hessianMatrix(p)
        hb = self.b.hessianMatrix(p)
        #todo
        h = ha * self.afactor  +  hb * self.bfactor
        check_matrix3_vectorized(h)
        return h

__all__ = ['SimpleBlend']
