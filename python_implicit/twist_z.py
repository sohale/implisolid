class Transformable1(object):
    pass

# does not work yet. Dont use.

from implicit_vectorized import ImplicitFunctionVectorized
import numpy as np
from basic_types import check_vector4_vectorized, make_vector4, check_vector4, check_scalar_vectorized, make_inverse

#ImplicitFuncitonVectorized3
class TwistZ(ImplicitFunctionVectorized, Transformable1):
    """ Twists alont the Z axis only, around Z axis.
    To change this, apply a matrix and its reverse (a separate class)
    Also see class @Transformed and @Screw."""


    def __init__(self, base_object, twist_rate):
        """twist_rate:  # cycles per mm
        """
        # axis_x0 is z always for now

        self.lamda = np.pi*2*twist_rate
        #if is_python3():
        #    super().__init__(initialMatrix=m)
        #else:
        #    super(Transformed, self).__init__(initialMatrix=m)

        assert issubclass(type(base_object), ImplicitFunctionVectorized)
        self.base_object = base_object

        #if m is None:
        #    m = np.eye(4)
        #check_matrix4(m)
        #self.matrix = m.copy()
        #self.invmatrix = make_inverse(m)

        #m = np.eye(4)
        #self.invmatrix = make_inverse(m)

        assert isinstance(self.base_object, ImplicitFunctionVectorized)

    def integrity_invariant(self):
        return True

    def implicitFunction(self, p):
        check_vector4_vectorized(p)
        # vec3.check_vector3_vectorized(p)
        N = p.shape[0]
        print("self.lamda", self.lamda)
        theta = p[:, 2] * self.lamda
        print(theta.shape)
        assert theta.shape == (N,)
        ca = np.cos(theta)
        sa = np.sin(theta)
        print(theta.shape, "theta")
        print(theta)

        #aa = ca[:, np.newaxis]
        #bb = p[:, 0, np.newaxis]
        #cc = - sa[:, np.newaxis]
        #dd = p[:, 1, np.newaxis]
        #print aa.shape, "aa"
        #print bb.shape, "bb"
        #print cc.shape, "cc"
        #print dd.shape, "dd"


        #x1 = ca[:, np.newaxis]*p[:, 0] - sa[:, np.newaxis]*p[:, 1]
        #x1 = aa*bb  # +cc*dd

        p2 = np.concatenate((
            ca[:, np.newaxis]*p[:, 0, np.newaxis] - sa[:, np.newaxis]*p[:, 1, np.newaxis],
            sa[:, np.newaxis]*p[:, 0, np.newaxis] + ca[:, np.newaxis]*p[:, 1, np.newaxis],
            p[:, 2, np.newaxis],  # z
            p[:, 3, np.newaxis]   # 1.
            ), axis=1)

        v = self.base_object.implicitFunction(p2)
        check_scalar_vectorized(v)
        return v

        #cs = np.concatenate((ca[:, np.newaxis], sa[:, np.newaxis]), axis=1)
        #assert cs.shape == (N, 2)
        # Nx2x2
        #cs_4 = np.concatenate( (
        #    ca[:, np.newaxis],
        #    -sa[:, np.newaxis],
        #    sa[:, np.newaxis],
        #    ca[:, np.newaxis]
        #    ), axis=1)
        assert cs_4.shape == (N, 4)
        cs_2x2 = cs_4.reshape(N, 2, 2)  #contains the reverse rotation matrix
        assert cs_2x2.shape == (N, 2, 2)
        #xy = np.dot(cs_2x2, p[:, 0:2])
        #assert xy.shape == (N, 2)
        p2 = p.copy()
        #p2[:, 0:2] = np.dot(cs_2x2, p[:,0:2])
        #m = np.tensordot(cs_2x2, p[:, 0:2], axes=(2, 1))

        #np.sum(a*b, axis=1)

        print("m.shape", m.shape)
        p2[:, 0:2] = m

        v = self.base_object.implicitFunction(p2)
        check_scalar_vectorized(v)
        return v

    #def implicitGradient(self, p):  # -> Vector3D :
        #check_vector4_vectorized(p)
        #tp = np.dot(self.invmatrix, vec3.make_v4(np.transpose(p)))
        #tp = np.transpose(tp)
        #g = self.base_object.implicitGradient(tp)
        #check_vector4_vectorized(g)
        # #g[:, 3] = 0  # important
        #v4 = np.dot(np.transpose(self.invmatrix), vec3.make_v4(np.transpose(g)))
        #v4 = np.transpose(v4)  # not efficient
        # #v4[:, 3] = 1
        #check_vector4_vectorized(v4)
        #return v4

    def implicitGradient__(self, p):  # -> Vector3D :
        check_vector4_vectorized(p)
        p = np.concatenate((p, np.ones((p.shape[0], 1))), axis=1)
        tp = np.dot(self.invmatrix, np.transpose(p))
        tp = np.transpose(tp)
        #tp = tp[:, :3]
        g = self.base_object.implicitGradient(tp)

        check_vector4_vectorized(g)

        g = np.concatenate((g, np.ones((g.shape[0], 1))), axis=1)
        v4 = np.dot(np.transpose(self.invmatrix), (np.transpose(g)))
        v4 = np.transpose(v4)
        #v3 = v4[:, :3]

        check_vector4_vectorized(v4)
        return v4

    def hessianMatrix(self, p):
        #warning: not tested
        check_vector4_vectorized(p)
        #tp = np.dot(self.invmatrix, vec3.make_v4(np.transpose(p)))
        #tp = np.transpose(tp)

        h1 = self.base_object.hessianMatrix(tp)
        #h = np.dot(h1, self.invmatrix)  # which one is correct?
        #h = np.dot(self.invmatrix, vec3.make_v4(np.tanspose(h1)))   # which one is correct?
        raise NotImplementedException()
        return h


__all__ = ['Ellipsoid', 'Transformed']
