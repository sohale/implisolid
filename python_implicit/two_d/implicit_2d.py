def check_vect2(x):
    assert x.ndim == 2
    assert x.shape[1] == 2
    assert x.shape[0] > 0

class Implicit2D(object):
    def implicitFunction(self, x):
        raise NotImplementedError()

    def implicitGradient(self, x):
        raise NotImplementedError()

    def hessianMatrix(self, x):
        raise NotImplementedError()

    def integrity_invariant(self):
        raise NotImplementedError()

def _slice_x(x, z):
    check_vect2(x)
    x3 = np.concatenate((x, z*ones(x.shape[0],1) ), axis=1)
    return x3

class Slice(Implicit2D):
    #not tested: quick draft
    def __init__(self, iobj3d, z):
        self.slice_z = z # todo: any plane U,V
        self.iobj3d = iobj3d

    def integrity_invariant(self):
        integrity = self.iobj3d.integrity_invariant()
        raise NotImplementedError()

    def implicitFunction(self, x):
        x3 = _slice_x(x)
        return iobj3d.implicitFunction(x3)

    def implicitGradient(self, x):
        x3 = _slice_x(x)
        g3d = iobj3d.implicitGradient(x3)
        #projected = 
        return projected

    def hessianMatrix(self, x):
        x3 = _slice_x(x)
        h3d = iobj3d.hessianMatrix(x3)
        return projected

