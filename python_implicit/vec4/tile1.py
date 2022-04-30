import numpy as np
from implicit_vectorized import ImplicitFunctionVectorized
from basic_types import check_vector4_vectorized, make_vector4, check_vector4, check_scalar_vectorized
#from numerical_utils import numerical_gradient
from implicit_config import config


def check_vector3(p):
    assert not issubclass(p.dtype.type, np.int), "vector cannot be integer"
    assert issubclass(p.dtype.type, np.float), "vector must be float"
    assert p.shape == (3,), "Vector must be a numpy array of (3) elements"
    assert not np.any( np.isnan(p.ravel()) )
    assert not np.any( np.isinf(p.ravel()) )


class Tile1D(ImplicitFunctionVectorized):   # Movable, uniform_scalable, rotatable

    def __init__(self, base_object, start, direction, tilecount, back_count=0):
        """
        """

        assert issubclass(type(base_object), ImplicitFunctionVectorized)
        self.base_object = base_object
        assert isinstance(self.base_object, ImplicitFunctionVectorized)

        (self.start, self.tilecount, self.back_count) = \
            (start, tilecount, back_count)

        check_vector3(start)
        check_vector3(direction)

        assert np.linalg.norm(direction) > config.numerical_min_length

        self.dir_normalised = direction / np.linalg.norm(direction)
        self.dir_len = np.linalg.norm(direction)

        #self.UVW = np.concatenate( (self.u, self.v, self.w), axis=1 )
        #self.UVW_inv = np.linalg.inv( self.UVW )

        assert float(int(self.tilecount)+1) == self.tilecount+1
        assert self.tilecount > 0
        assert float(int(self.back_count)+1) == self.back_count+1
        assert self.back_count >= 0

        assert np.isscalar(self.dir_len)
        assert self.dir_len > 0
        assert self.integrity_invariant()
        assert self.dir_len > config.numerical_min_length

    def integrity_invariant(self):
        #config = threeD_printing_config_profile
        norm_tol = 0.00000001  # 1000 km
        matrix_inv_tol = 0.000001

        #sane = True
        sd = {"sane": True}
        def check(boolean, reason):
            sd["sane"] = sd["sane"] and boolean
            if not boolean:
                print("Error:", reason)
            pass
        check(np.abs(np.linalg.norm(self.dir_normalised)-1.0) < norm_tol, "dir_normalised")
        check(self.start.shape == (3,), "3")
        check(self.dir_normalised.shape == (3,), "3")
        check(np.isscalar(self.dir_len), "")
        #check(np.allclose(np.dot(self.UVW, self.UVW_inv), np.eye(3), matrix_inv_tol), "MMinv")
        #check(np.allclose(np.dot(self.UVW_inv, self.UVW), np.eye(3), matrix_inv_tol), "MinvM")
        check(self.dir_len > config.numerical_min_length, "minlen")  # functions, not actual models
        check(self.dir_len > 0, "dir_len +")
        check(np.floor(self.tilecount) == self.tilecount, "int")
        check( self.tilecount > 0, "point   count must be a natural number (positive integer)" )
        check(np.floor(self.back_count) == self.back_count, "int")
        check( self.back_count >= 0, "point   count must be a natural number (positive integer)" )
        check( np.isscalar(self.dir_len), "scalar" )
        check( self.dir_len > config.numerical_min_length, "vector too small")

        return sd["sane"]

    def _move(self, xa):
        x = xa[:, 0:3].copy()
        pointcount = x.shape[0]

        assert self.dir_normalised.ndim == 1
        assert x.shape[1] == self.dir_normalised.size
        x = x - self.start[np.newaxis, :]
        la = np.dot(x, self.dir_normalised)
        assert la.shape == (pointcount,)
        la = np.floor(la/self.dir_len)
        print la.shape
        la[la > self.tilecount-1-self.back_count + 0.1] = 0
        la[la < -self.back_count - 0.1] = 0
        la = la*self.dir_len

        print self.dir_normalised.shape
        x = x - np.dot(la[:, np.newaxis], self.dir_normalised[np.newaxis, :])  # inplace
        x = x + self.start[np.newaxis, :]
        assert x.shape == (pointcount, 3)
        x4 = np.concatenate((x, np.ones((pointcount, 1))), axis=1)
        check_vector4_vectorized(x4)
        return x4

    #alternative design:
    #def transform(self, xa, return_grad=False):
    #    return self._move(self, xa)

    def implicitFunction(self, xa, return_grad=False):
        check_vector4_vectorized(xa)
        x4 = self._move(xa)
        #return x4
        v = self.base_object.implicitFunction(x4)
        check_scalar_vectorized(v)
        return v

        if not return_grad:
            return v
        else:
            raise NotImplementedError()
            check_vector4_vectorized(g4)
            return v, g4

    def implicitGradient(self, x):
        check_vector4_vectorized(xa)
        x4 = self._move(xa)
        g4 = self.base_object.implicitGradient(x4)
        check_vector4_vectorized(g4)
        return g4

    def curvature(self, x):
        check_vect2(x)
        raise NotImplementedError()


__all__ = ['Tile1D']
