import numpy as np
from implicit_2d import Implicit2D

#not tested: quick draft
class UnitCircle(Implicit2D):
    def implicitFunction(self, x):
        check_vect2(x)
        print("not tested")
        return 1 - np.sum(x*x, axis=1)

    def implicitGradient(self, x):
        check_vect2(x)
        print("not tested")
        return -2*x

    def curvature(self, x):
        check_vect2(x)
        print("not tested")
        return -2*np.eye(2)


#not tested: quick draft
class BitMapImplicit2D(Implicit2D):
    def __init__(self, filename, w, h, angle):
        self.sz = np.array([w,h])
        raise NotImplementedError()

    def implicitFunction(self, x):
        check_vect2(x)
        raise NotImplementedError()

    def implicitGradient(self, x):
        check_vect2(x)
        raise NotImplementedError()

    def curvature(self, x):
        check_vect2(x)
        raise NotImplementedError()


#not tested: quick draft
class Conic(Implicit2D):
    def __init__(self):
        raise NotImplementedError()

    def implicitFunction(self, x):
        check_vect2(x)
        raise NotImplementedError()

    def implicitGradient(self, x):
        check_vect2(x)
        raise NotImplementedError()

    def curvature(self, x):
        check_vect2(x)
        raise NotImplementedError()
