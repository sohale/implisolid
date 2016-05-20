import numpy as np

from primitives import ImplicitFunctionPointwise
from basic_functions import make_inverse, check_matrix4, is_python3
from basic_functions import check_vector3

# from lib.transformations import *
import transformations as tf


class Transformable(object):
    """ #interface only. Use Builder pattern."""

    def __init__(self, initialTransformable=None, initialMatrix=None):
        """ Usage:
            t = Transformatble()  or t = Transformable(t) or t = Transformable(m) """
        assert initialTransformable is None or initialMatrix is None

        if initialMatrix is not None:
            check_matrix4(initialMatrix)
            matrix = initialMatrix
        elif initialTransformable is not None:
            assert issubclass(type(initialTransformable), Transformable)
            assert initialTransformable is not None
            matrix = initialTransformable.matrix
        else:
            matrix = np.eye(4)

        check_matrix4(matrix)

        self.matrix = matrix
        self.invmatrix = make_inverse(self.matrix)

    def rotate(self, angle, along, units="rad"):
        if units == "rad":
            pass
        elif units == "deg":
            angle = angle / 360.0 * np.pi*2.0
        else:
            raise ValueError()

        check_vector3(along)
        rm = tf.rotation_matrix(angle, along)
        self.matrix = np.dot(rm, self.matrix)
        self.invmatrix = make_inverse(self.matrix)

        return self

    def move(self, x, y, z):
        self.matrix[:, 3] += [x, y, z, 0]
        self.invmatrix = make_inverse(self.matrix)
        return self

    def resize(self, s):
        self.matrix[:, 0:3] *= s
        self.invmatrix = make_inverse(self.matrix)
        return self


class Transformed(ImplicitFunctionPointwise, Transformable):
    """ See the @Ellipsoid class. Just replace base_sphere with base_object. Note: super and self have shared members self.matric and self.invmatrix."""
    def __init__(self, base_object, m=None):
        """ """
        if is_python3():
            super().__init__(initialMatrix=m)
        else:
            super(Transformed, self).__init__(initialMatrix=m)

        assert issubclass(type(base_object), ImplicitFunctionPointwise)
        self.base_object = base_object

        if m is None:
            m = np.eye(4)
        check_matrix4(m)
        self.matrix = m.copy()
        self.invmatrix = make_inverse(m)

        assert isinstance(self.base_object, ImplicitFunctionPointwise)

    def implicitFunction(self, p):
        check_vector3(p)
        p = np.concatenate((p, np.ones((p.shape[0], 1))), axis=1)
        tp = np.dot(self.invmatrix, p)
        v = self.base_object.implicitFunction(tp)
        return v

    def implicitGradient(self, p):
        check_vector3(p)
        p = np.concatenate((p, np.ones((p.shape[0], 1))), axis=1)
        tp = np.dot(self.invmatrix, p)
        g = self.base_object.implicitGradient(tp)
        v3 = np.dot(np.transpose(self.invmatrix), g)
        return v3


__all__ = ['Transformable', 'Transformed']
