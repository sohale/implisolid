import numpy as np
#from implicit_config import TOLERANCE

#from implicit_config import VERBOSE

#from basic_types import make_inverse, check_matrix4, make_vector4
from basic_types import check_vector4_vectorized, check_matrix3_vectorized
from basic_types import make_vector4, check_vector4  # check_matrix4
from numerical_utils import numerical_gradient_vectorized_v1

# @profile
# def memoize(f):
#     memo = {}
#     def inner(*args):
#         if str(args) not in memo:
#             memo[str(args)] = f(*args)
#         return memo[str(args)]
#     return inner


class ImplicitFunctionVectorized(object):
    """ Functions in this type receive numpy vectors of size Nx4 """

    def implicitFunction(self, pv):
        raise VirtualException()

    def implicitGradient(self, pv):
        """ Returns a vector of size N x 4 where N is the number of points. result[:,3] must be 1.0 for all. """
        #raise VirtualException()
        print("Warning: Numerical gradient is used")
        return numerical_gradient_vectorized_v1(self, pv)

    def hessianMatrix(self, pv):
        """ Returns a vector of size N x 4 x 4 where N is the number of points"""
        raise VirtualException()

    def integrity_invariant(self):
        return False

    def assert_integrity(self):
        """ This is not an assert. Unless an assert, it has to be called in production, at runtime and through exception. """
        ok = self.integrity_invariant()
        if not ok:
            # todo: provide more information
            raise ObjectError()


class SignedDistanceImplicitVectorized(object):
    """ Interface only (can be used in multiple inheritance) """
    pass
