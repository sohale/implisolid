from implicit import ImplicitFunction
import numpy as np
from basic_functions import check_scalar_vectorized, check_vector3_vectorized


class Transformable1(object):
    pass


# ImplicitFuncitonVectorized3
class TwistZ(ImplicitFunction, Transformable1):
    """ Twists alont the Z axis only, around Z axis.
    To change this, apply a matrix and its reverse (a separate class)
    Also see class @Transformed and @Screw."""

    def __init__(self, base_object, twist_rate):
        """twist_rate:  # cycles per mm
        """
        # axis_x0 is z always for now

        self.lamda = np.pi*2*twist_rate
        # if is_python3():
        #    super().__init__(initialMatrix=m)
        # else:
        #    super(Transformed, self).__init__(initialMatrix=m)

        assert issubclass(type(base_object), ImplicitFunction)
        self.base_object = base_object

        # if m is None:
        #    m = np.eye(4)
        # check_matrix4(m)
        # self.matrix = m.copy()
        # self.invmatrix = make_inverse(m)

        assert isinstance(self.base_object, ImplicitFunction)

    def integrity_invariant(self):
        return True

    def implicitFunction(self, p):
        check_vector3_vectorized(p)
        # vec3.check_vector3_vectorized(p)
        N = p.shape[0]
        print "self.lamda", self.lamda
        theta = p[:, 2] * self.lamda
        print theta.shape
        assert theta.shape == (N,)
        ca = np.cos(theta)
        sa = np.sin(theta)
        print theta.shape, "theta"
        print theta

        # aa = ca[:, np.newaxis]
        # bb = p[:, 0, np.newaxis]
        # cc = - sa[:, np.newaxis]
        # dd = p[:, 1, np.newaxis]
        # print aa.shape, "aa"
        # print bb.shape, "bb"
        # print cc.shape, "cc"
        # print dd.shape, "dd"
        # x1 = ca[:, np.newaxis]*p[:, 0] - sa[:, np.newaxis]*p[:, 1]
        # x1 = aa*bb  # +cc*dd

        p2 = np.concatenate((
            ca[:, np.newaxis]*p[:, 0, np.newaxis] - sa[:, np.newaxis]*p[:, 1, np.newaxis],
            sa[:, np.newaxis]*p[:, 0, np.newaxis] + ca[:, np.newaxis]*p[:, 1, np.newaxis],
            p[:, 2, np.newaxis],), axis=1)

        v = self.base_object.implicitFunction(p2)
        check_scalar_vectorized(v)
        return v
        p2 = p.copy()
        p2[:, 0:2] = m
        print "m.shape", m.shape

        v = self.base_object.implicitFunction(p2)
        check_scalar_vectorized(v)
        return v

    def implicitGradient(self, p):  # -> Vector3D :

        return None

    def hessianMatrix(self, p):
        # warning: not tested
        # check_vector3_vectorized(p)
        # tp = np.dot(self.invmatrix, vector3.make_vector4(np.transpose(p)))
        # tp = np.transpose(tp)
        # h1 = self.base_object.hessianMatrix(tp)
        # h = np.dot(h1, self.invmatrix)  # which one is correct?
        # h = np.dot(self.invmatrix, vec3.make_v4(np.tanspose(h1)))   # which one is correct?
        raise Exception()
    #    return h1


__all__ = ['Ellipsoid', 'Transformed']
