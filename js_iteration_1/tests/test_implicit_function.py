import numpy as np

class ImplicitFunction(object):
    """ Functions in this type receive numpy vectors of size Nx4 """

    def implicitFunction(self, pv):
        raise Exception()

    def implicitGradient(self, pv):
        """ Returns a vector of size N x 3 where N is the number of points.  """
        raise Exception()

    def hessianMatrix(self, pv):
        """ Returns a vector of size N x 3 x 3 where N is the number of points"""
        raise Exception()

    def integrity_invariant(self):
        return False


class UnitSphere(ImplicitFunction):
    def implicitFunction(self, pv):

        return 0.25 - np.sum(pv * pv, axis=1)

    def implicitGradient(self, pv):
        assert pv.ndim == 2
        assert pv.shape[1:] == (3,)
        grad = -2*pv

        return grad


def make_vector3(x, y, z):
    xyz = np.array([x, y, z])
    assert not np.any(np.isnan(xyz))
    assert not np.any(np.isinf(xyz))

    if issubclass(type(x), np.ndarray):
        return np.array((float(x[0]), float(y[1]), float(z[1])))

    return np.array((float(x), float(y), float(z)))


class UnitCube1(ImplicitFunction):
    def __init__(self, size=0.5):
        self.p0 = []
        self.n0 = []

        def side(x, y, z):
            p0 = make_vector3(x, y, z)
            p0 = p0 / 2.0 * size
            n0 = -make_vector3(x, y, z)
            n0 = n0
            self.p0 += [p0]
            self.n0 += [n0]
            # print(self.p0[-1])

            def norm2(v):
                return v[0]*v[0]+v[1]*v[1]+v[2]*v[2]
            assert norm2(self.n0[-1]) - 1. == 0.0

        side(1, 0, 0)
        side(-1, 0, 0)
        side(0, 1, 0)
        side(0, -1, 0)
        side(0, 0, 1)
        side(0, 0, -1)

    def implicitFunction(self, p):

        sides = len(self.p0)
        n = p.shape[0]
        temp = np.zeros((n, sides))
        for i in range(sides):
            p0 = self.p0[i]
            n0 = self.n0[i]
            sub = p - np.tile(p0[np.newaxis, :], (n, 1))
            vi = np.dot(sub, n0)
            temp[:, i] = vi
        va = np.amin(temp, axis=1)
        return va

    def implicitGradient(self, p):

        sides = 6
        na = np.zeros((sides, 3))
        n = p.shape[0]
        temp = np.zeros((n, sides))
        for i in range(len(self.p0)):
            p0 = self.p0[i]
            n0 = self.n0[i]
            sub = p - np.tile(p0[np.newaxis, :], (n, 1))
            vi = np.dot(sub, n0)

            temp[:, i] = vi

            na[i, :] = n0

        ia = np.argmin(temp, axis=1)

        assert ia.shape == (n,)

        g = na[ia, :]

        return g


def sphere_example(scale=1.):
    iobj = UnitCube1()

    return iobj


def make_example_vectorized(name, scale=1.0):
    res = globals()[name](scale)
    return res


def vtk_mc(gridvals, rrr):
    (RANGE_MIN, RANGE_MAX, STEPSIZE) = rrr

    # from vtk import *

    from vtk import vtkImageData, vtkContourFilter

    data_numpy = gridvals

    spacing = 1

    x_dim, y_dim, z_dim = data_numpy.shape

    data_numpy = data_numpy.ravel(order='C')

    from vtk.util.numpy_support import numpy_to_vtk
    numpy_to_vtk_data = numpy_to_vtk(data_numpy)

    vtk_image_grid = vtkImageData()  # sample , vtk_Data
    vtk_image_grid.SetExtent(0, x_dim-1, 0, y_dim-1, 0, z_dim-1)
    vtk_image_grid.SetSpacing(spacing, spacing, spacing)  # check this
    vtk_image_grid.SetOrigin(0, 0, 0)
    vtk_image_grid.SetNumberOfScalarComponents(1)
    vtk_image_grid.Update()

    vtk_image_grid.GetPointData().SetScalars(numpy_to_vtk_data)

    contour = vtkContourFilter()
    contour.SetInput(vtk_image_grid)
    contour.GenerateValues(1, 0, 1)
    contour.Update()

    mesh_data = contour.GetOutput()

    vcount = mesh_data.GetNumberOfPoints()
    va = np.zeros((vcount, 3))
    for i in range(vcount):
        p = mesh_data.GetPoint(i)
        va[i, :] = p
        va[i, :] = va[i, :] * STEPSIZE + RANGE_MIN

    fcount = mesh_data.GetNumberOfCells()
    fa = np.zeros((fcount, 3), dtype=int)
    for i in range(fcount):
        c = mesh_data.GetCell(i)
        (fa[i, 0], fa[i, 1], fa[i, 2]) = (c.GetPointId(0), c.GetPointId(1), c.GetPointId(2))

    va = va[:, [1, 0, 2]]

    return va, fa


def _prepare_grid(rng):
    """
    rng: like the output of np.arange()
    """
    assert rng.size < 200

    (yy, xx, zz) = np.meshgrid(rng, rng, rng)
    xyz_nparray = np.transpose(np.vstack([xx.ravel(), yy.ravel(), zz.ravel()]))
    assert xyz_nparray.shape[1:] == (3,)

    return xyz_nparray


def make_grid(iobj, rng, old=None, return_xyz=False):
    xyz_nparray = _prepare_grid(rng)

    vgrid_v = iobj.implicitFunction(xyz_nparray)

    f_out = open("/home/solene/Desktop/mp5-private/implisolid/js_iteration_1/tests/value_implicit_function_python.txt", "w")
    f_out.write("Implicit values :" + '\n')
    for i in range(vgrid_v.shape[0]):
        f_out.write(str(vgrid_v[i]) + '\n')
    f_out.close()
    vgrid = np.reshape(vgrid_v, (len(rng), len(rng), len(rng)), order='C')

    if return_xyz:
        return vgrid, xyz_nparray
    else:
        return vgrid


def make_mc_values_grid(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE, old=True):
    rng = np.arange(RANGE_MIN, RANGE_MAX, STEPSIZE)
    vgrid = make_grid(iobj, rng, old=old)
    if old:
        return np.swapaxes(vgrid, 0, 1)
    else:
        print ("*********************************************")
        vgrid = np.swapaxes(vgrid, 1, 2)
        vgrid = np.swapaxes(vgrid, 0, 1)
        return vgrid


def demo_sphere():
    object_name = "sphere_example"
    iobj = make_example_vectorized(object_name)

    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-1, +1, 0.4)

    gridvals = make_mc_values_grid(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE, old=False)
    vertex, faces = vtk_mc(gridvals, (RANGE_MIN, RANGE_MAX, STEPSIZE))


if __name__ == '__main__':

    demo_sphere()
