__all__ = []

from ellipsoid import *   # actual call
import ellipsoid  # access barrier
__all__.extend(ellipsoid.__all__)


from primitives import *
import primitives
__all__.extend(primitives.__all__)


from crisp_csg import *
import crisp_csg
__all__.extend(crisp_csg.__all__)

from transformed import *
import transformed
__all__.extend(transformed.__all__)

from simple_blend_nonvec import *
import simple_blend_nonvec
__all__.extend(simple_blend_nonvec.__all__)


def is_implicit_type(iobj_nonvec):
    return issubclass(type(iobj_nonvec), ImplicitFunctionPointwise)

__all__.extend(["is_implicit_type"])
