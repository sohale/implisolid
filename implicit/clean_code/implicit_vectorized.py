import numpy as np
#from implicit_config import TOLERANCE

#from implicit_config import VERBOSE

#from basic_functions import make_inverse, check_matrix4, make_vector4
from basic_functions import check_vector3_vectorized, check_matrix3_vectorized

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
        """ Returns a vector of size N x 3 where N is the number of points.  """
        raise VirtualException()

    def hessianMatrix(self, pv):
        """ Returns a vector of size N x 3 x 3 where N is the number of points"""
        raise VirtualException()

    def integrity_invariant(self):
        return False


class SignedDistanceImplicitVectorized(object):
    """ Interface only (can be used in multiple inheritance) """
    pass


class UnitSphere(ImplicitFunctionVectorized):
    # @profile
    def implicitFunction(self, pv):

        return 1.0 - np.sum(pv * pv, axis=1)
    # @profile
    def implicitGradient(self, pv):
        assert pv.ndim == 2
        assert pv.shape[1:] == (3,)
        grad = -2*pv

        return grad

    def hessianMatrix(self, vp):
        n = vp.shape[0]
        h = np.array([[[-2, 0, 0],  [0, -2, 0],  [0, 0, -2]]], ndmin=3)   # 1x 3x3
        print(h)
        print(h.shape)
        check_matrix3_vectorized(h)
        vh = np.expand_dims(h, axis=0)  # insert one dimension: 4x4 -> N x (4x4)
        return h


class UnitCube1(ImplicitFunctionVectorized, SignedDistanceImplicitVectorized):
    def __init__(self, size=1.0):
        self.p0 = []
        self.n0 = []

        def side(x, y, z):
            p0 = (make_vector3(x, y, z) + 0.0)
            p0 = p0 / 2.0 * size
            n0 = -make_vector3(x, y, z)
            n0 = n0
            self.p0 += [p0]
            self.n0 += [n0]
            #print(self.p0[-1])
            check_vector3(self.p0[-1])
            check_vector3(self.n0[-1])

            def norm2(v):
                return v[0]*v[0]+v[1]*v[1]+v[2]*v[2]
            assert norm2(self.n0[-1]) - 1 == 0.0

        side(1, 0, 0)
        side(-1, 0, 0)
        side(0, 1, 0)
        side(0, -1, 0)
        side(0, 0, 1)
        side(0, 0, -1)

    #todo: evaluate gradient and implicit function at the same time
    # @memoize
    # @profile
    def implicitFunction(self, p):
        check_vector3_vectorized(p)

        sides = len(self.p0)
        n = p.shape[0]
        temp = np.zeros((n, sides))
        for i in range(sides):
            p0 = self.p0[i]
            n0 = self.n0[i]
            sub = p - np.tile(p0[np.newaxis, :], (n, 1))
            # assert np.allclose(sub[:, 3], 0)
            vi = np.dot(sub, n0)
            temp[:, i] = vi
        #ia = np.argmin(temp, axis=1)
        va = np.amin(temp, axis=1)
        return va

        """
        n = p.shape[0]
        #temp = np.zeros((n,6,3))
        mp0 = np.zeros((1,6,3))
        for i in range(len(self.p0)):
            p0 = self.p0[i]
            #n0 = self.n0[i]
            mp0[0,i,:] = - p0

        #temp = np.zeros((n,6,3))
        temp = np.tile(p[:,np.newaxis,:], (1,6,1))
        temp = temp - np.tile(mp0, (n,:,:))  # x-p0
        temp = np.dot( temp, n0)
        for i in range(len(self.p0)):
            p0 = self.p0[i]
            n0 = self.n0[i]

            temp[:,i,:] = p[i,0:3] - p0
            #vi = np.dot(x-p0, n0)
        pass
        """
        """
        chosen_i = None
        v = +np.infty
        for i in range(len(self.p0)):
            p0 = self.p0[i]
            n0 = self.n0[i]
            vi = np.dot( p-p0, n0 )
            if vi < v:
                v = vi
                chosen_i = i


        assert not chosen_i is None

        return v
        """
    # @memoize
    # @profile
    def implicitGradient(self, p):
        check_vector3_vectorized(p)

        sides = 6
        na = np.zeros((sides, 3))
        n = p.shape[0]
        temp = np.zeros((n, sides))
        for i in range(len(self.p0)):
            p0 = self.p0[i]
            n0 = self.n0[i
            sub = p - np.tile(p0[np.newaxis, :], (n, 1))
            vi = np.dot(sub, n0)
            #print(vi)
            temp[:, i] = vi

            na[i, :] = n0

        ia = np.argmin(temp, axis=1)
        #va = np.amin(temp, axis=1)
        assert ia.shape == (n,)

        g = na[ia, :]

        check_vector3_vectorized(g)
        return g


        """
        check_vector4(p)
        chosen_i = None
        v = +np.infty
        for i in range(len(self.p0)):
            p0 = self.p0[i]
            n0 = self.n0[i]
            vi = np.dot( p-p0, n0 )
            if vi < v:
                v = vi
                chosen_i = i
                grad = n0

        assert not chosen_i is None

        #not tested
        grad[3] = 1
        check_vector4(grad)
        return grad
        """
    # @memoize
    # @profile
    def hessianMatrix(self, p):
        check_vector3(p)
        h = np.array([[0, 0, 0],  [0, 0, 0],  [0, 0, 0]], ndmin=2)
        check_matrix3(h)
        #todo : array of matrix3s
        n = p.shape[0]
        ha = np.tile(h[np.newaxis, :, :], (n, 1, 1))
        return h


__all__ = ['ImplicitFunctionVectorized', 'UnitSphere', 'UnitCube1']

if __name__ == "__main__":
    pass
