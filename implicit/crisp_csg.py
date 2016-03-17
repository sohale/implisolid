from basic_types import check_vector4, check_matrix3
#from implicit_config import TOLERANCE
from primitives import ImplicitFunctionPointwise

#from primitives import UnitSphere


class CrispUnion(ImplicitFunctionPointwise):
    def __init__(self, a, b):
        self.a = a
        self.b = b
        assert isinstance(a, ImplicitFunctionPointwise)
        assert isinstance(b, ImplicitFunctionPointwise)

    def implicitFunction(self, p):
        check_vector4(p)
        va = self.a.implicitFunction(p)
        vb = self.b.implicitFunction(p)
        if va > vb:
            v = va
        else:
            v = vb
        return v

    def implicitGradient(self, p):
        check_vector4(p)
        va = self.a.implicitFunction(p)
        vb = self.b.implicitFunction(p)
        if va > vb:
            grad = self.a.implicitGradient(p)
        else:
            grad = self.b.implicitGradient(p)
        check_vector4(grad)
        return grad

    def hessianMatrix(self, p):
        check_vector4(p)
        va = self.a.implicitFunction(p)
        vb = self.b.implicitFunction(p)
        if va > vb:
            h = self.a.hessianMatrix(p)
        else:
            h = self.b.hessianMatrix(p)
        check_matrix3(h)
        return h


class CrispIntersection(ImplicitFunctionPointwise):
    def __init__(self, a, b):
        self.a = a
        self.b = b

    def implicitFunction(self, p):
        check_vector4(p)
        va = self.a.implicitFunction(p)
        vb = self.b.implicitFunction(p)
        if va > vb:
            v = vb
        else:
            v = va
        return v

    def implicitGradient(self, p):
        check_vector4(p)
        va = self.a.implicitFunction(p)
        vb = self.b.implicitFunction(p)
        if va > vb:
            grad = self.b.implicitGradient(p)
        else:
            grad = self.a.implicitGradient(p)
        check_vector4(grad)
        return grad

    def hessianMatrix(self, p):
        check_vector4(p)
        va = self.a.implicitFunction(p)
        vb = self.b.implicitFunction(p)
        if va > vb:
            h = self.b.hessianMatrix(p)
        else:
            h = self.a.hessianMatrix(p)
        check_matrix3(h)
        return h


class CrispSubtract(ImplicitFunctionPointwise):
    def __init__(self, a, b):
        self.a = a
        self.b = b
        assert isinstance(a, ImplicitFunctionPointwise)
        assert isinstance(b, ImplicitFunctionPointwise)

    def implicitFunction(self, p):
        check_vector4(p)
        va = self.a.implicitFunction(p)
        vb = - self.b.implicitFunction(p)
        if va > vb:
            v = vb
        else:
            v = va
        return v

    def implicitGradient(self, p):
        check_vector4(p)
        va = self.a.implicitFunction(p)
        vb = - self.b.implicitFunction(p)
        if va > vb:
            grad = -self.b.implicitGradient(p)
            grad[3] = 1
        else:
            grad = self.a.implicitGradient(p)
        check_vector4(grad)
        return grad

    def hessianMatrix(self, p):
        check_vector4(p)
        va = self.a.implicitFunction(p)
        vb = - self.b.implicitFunction(p)
        if va > vb:
            h = -self.b.hessianMatrix(p)
        else:
            h = self.a.hessianMatrix(p)
        check_matrix3(h)
        return h


__all__ = ['CrispUnion', 'CrispIntersection', 'CrispSubtract']
