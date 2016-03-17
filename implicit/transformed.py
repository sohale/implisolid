import numpy as np

from primitives import ImplicitFunctionPointwise
from basic_types import make_inverse, check_vector4, check_matrix4, is_python3

#from lib.transformations import *
from lib import transformations as tf

class Transformable(object):
    """ #interface only. Use Builder pattern."""

    def __init__(self, initialTransformable=None, initialMatrix=None):
        """ Usage:
            t = Transformatble()  or t = Transformable(t) or t = Transformable(m) """
        assert initialTransformable is None or initialMatrix is None

        if not initialMatrix is None:
            check_matrix4(initialMatrix)
            matrix = initialMatrix
        elif not initialTransformable is None:
            assert issubclass(type(initialTransformable), Transformable)
            assert not initialTransformable is None
            matrix = initialTransformable.matrix
        else:
            matrix = np.eye(4)

        check_matrix4(matrix)

        self.matrix = matrix
        self.invmatrix = make_inverse(self.matrix)

    #def rotate_euler(x, y, z, type="EulerXYZ"):
    #    assert type == "EulerXYZ"
    #    raise
    #    return self

    def rotate(self, angle, along, units="rad"):
        if units == "rad":
            pass
        elif units == "deg":
            angle = angle / 360.0 * np.pi*2.0
        else:
            raise UsageError()  # UsageError

        check_vector4(along)
        rm = tf.rotation_matrix(angle, along[0:3])
        self.matrix = np.dot(rm , self.matrix)
        self.invmatrix = make_inverse(self.matrix)

        #print(angle /(3.1415926536*2) * 360 )
        #print(rm)
        #print(self.matrix)
        return self

    def move(self, x, y, z):
        self.matrix[:, 3] += [x, y, z, 0]
        self.invmatrix = make_inverse(self.matrix)
        return self

    def resize(self, s):
        self.matrix[:, 0:3] *= s
        self.invmatrix = make_inverse(self.matrix)
        return self

#class Transformed(ImplicitFunctionPointwise, Transformable):


class Transformed(ImplicitFunctionPointwise, Transformable):
    """ See the @Ellipsoid class. Just replace base_sphere with base_object. Note: super and self have shared members self.matric and self.invmatrix."""
    #Rotated. LinearTransformation. HomogeneousCoordinates
    def __init__(self, base_object, m=None):
        """ """
        if is_python3():
            super().__init__(initialMatrix=m)
        else:
            super(Transformed, self).__init__(initialMatrix=m)

        #assert type(baseImplicitClass) is type
        #assert issubclass(baseImplicitClass, ImplicitFunctionPointwise)
        #self.base_object = baseImplicitClass()
        assert issubclass(type(base_object), ImplicitFunctionPointwise)
        self.base_object = base_object

        if m is None:
            m = np.eye(4)
        check_matrix4(m)
        self.matrix = m.copy()
        self.invmatrix = make_inverse(m)

        assert isinstance(self.base_object, ImplicitFunctionPointwise)

    def implicitFunction(self, p):
        check_vector4(p)
        tp = np.dot(self.invmatrix, p)
        v = self.base_object.implicitFunction(tp)
        return v

    def implicitGradient(self, p):  # -> Vector3D :
        check_vector4(p)
        tp = np.dot(self.invmatrix, p)
        g = self.base_object.implicitGradient(tp)
        g[3] = 0  # important
        v4 = np.dot(np.transpose(self.invmatrix), g)
        v4[3] = 1
        return v4

    def hessianMatrix(self, p):
        #warning: not tested
        check_vector4(p)
        tp = np.dot(self.invmatrix, p)
        h1 = self.base_object.hessianMatrix(tp)
        h = np.dot(h1, self.invmatrix)  # which one is correct?
        h = np.dot(self.invmatrix, h1)   # which one is correct?
        raise VirtualException()
        return h

__all__ = ['Transformable', 'Transformed']
