from implicit import ImplicitFunction
import numpy as np
from basic_functions import check_scalar_vectorized, check_vector3_vectorized


class Transformable1(object):
    pass


class TwistZ(ImplicitFunction, Transformable1):
    """ Twists alont the Z axis only, around Z axis.
    To change this, apply a matrix and its reverse (a separate class)
    Also see class @Transformed and @Screw."""

    def __init__(self, base_object, twist_rate):
        """twist_rate:  # cycles per mm
        """

        self.lamda = np.pi*2*twist_rate

        assert issubclass(type(base_object), ImplicitFunction)
        self.base_object = base_object

        assert isinstance(self.base_object, ImplicitFunction)

    def integrity_invariant(self):
        return True

    def implicitFunction(self, p):
        check_vector3_vectorized(p)

        N = p.shape[0]
        print "self.lamda", self.lamda
        theta = p[:, 2] * self.lamda
        print theta.shape
        assert theta.shape == (N,)
        ca = np.cos(theta)
        sa = np.sin(theta)
        print theta.shape, "theta"
        print theta

        p2 = np.concatenate((
            ca[:, np.newaxis]*p[:, 0, np.newaxis] - sa[:, np.newaxis]*p[:, 1, np.newaxis],
            sa[:, np.newaxis]*p[:, 0, np.newaxis] + ca[:, np.newaxis]*p[:, 1, np.newaxis],
            p[:, 2, np.newaxis],), axis=1)

        v = self.base_object.implicitFunction(p2)
        check_scalar_vectorized(v)
        return v

    def implicitGradient(self, p):  # -> Vector3D :

        return None


__all__ = ['Ellipsoid', 'Transformed']
