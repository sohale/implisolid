import numpy as np
from basic_types import check_vector4_vectorized, check_scalar_vectorized
from basic_types import check_matrix3_vectorized

from implicit_vectorized import ImplicitFunctionVectorized
from numerical_utils import numerical_gradient_slow_func


#R-functions were introduced for better CSG opeations in 1963 by Vladimir Logvinovich Rvachev
#R-disjunction: Union (or)
#R-conjunction: Intersection (and)
#todo: RFunction in its general form

# Incomplete. Not tested.

class RSubtract(ImplicitFunctionVectorized):
    """ CSG Subtract based on Rvachev functions. If alpha = 1 => we have equivalent to CrispSubtract """

    def __init__(self, a, b, alpha):
        """ Special values for alpha: alpha==1.0 => we are CrispSubtract """
        assert isinstance(a, ImplicitFunctionVectorized)
        assert isinstance(b, ImplicitFunctionVectorized)
        self.a = a
        self.b = b
        self.alpha = alpha
        self.sgn = -1  # -1 = conjunction
        self.integrity_invariant()

    def integrity_invariant(self):
        return self.alpha  and (self.sgn == -1 or self.sgn == 1)

    def implicitFunction(self, p):
        check_vector4_vectorized(p)
        va = self.a.implicitFunction(p)
        vb = - self.b.implicitFunction(p)

        v = (1.0/(1+self.alpha))*(va + vb + self.sgn * np.sqrt(va**2 + vb**2 - (2*self.alpha)*va*vb))
        check_scalar_vectorized(v)
        return v

    def implicitGradient(self, x):
        """also see the class Screw. same code."""
        return numerical_gradient_slow_func(self, x)

    def hessianMatrix(self, p):
        raise NotImplementedError()


__all__ = ['RSubtract',]
