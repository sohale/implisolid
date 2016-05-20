__all__ = []

import implicit
from implicit import *
__all__.extend(implicit.__all__)


from crisp_csg import *
import crisp_csg
__all__.extend(crisp_csg.__all__)


from basic_functions import repeat_vect3, check_vector3_vectorized
__all__.extend(['repeat_vect3', 'check_vector3_vectorized'])

from simple_blend import *
import simple_blend
__all__.extend(simple_blend.__all__)

#from simple_blend import *
from ellipsoid import *
import ellipsoid
__all__.extend(ellipsoid.__all__)



from screw import *
import screw
__all__.extend(screw.__all__)

# from rvachev_csg import *
# import rvachev_csg
# __all__.extend(rvachev_csg.__all__)

from cylinder_simple import *
import cylinder_simple
__all__.extend(cylinder_simple.__all__)


def is_implicit_type(iobj_vec):
    return issubclass(type(iobj_vec), ImplicitFunctionVectorized)

__all__.extend(["is_implicit_type"])


# from implicit_vectorized import *
# import implicit_vectorized
# __all__.extend(implicit_vectorized.__all__)


# __all__.extend(ellipsoid_vectorized.__all__)
# __all__.extend(simple_blend.__all__)


# from basic_functions import *


# print(crisp_csg_vectorized.CrispSubtract)

# __all__.extend([
#    'normalize_vector4_vectorized',
#    'make_random_vector_vectorized',
#    'almost_equal4_vectorized',
#    'check_scalar_vectorized',
#    'check_vector4_vectorized',
#    'check_matrix3_vectorized',
#    'check_matrix4_vectorized'
#    ])
