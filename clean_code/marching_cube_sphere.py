import numpy as np
import sys
from mayavi import mlab


def display_simple_using_mayavi_2(vf_list, pointcloud_list, minmax=(-1, 1), mayavi_wireframe=False, opacity=1.0, separate=True,
 gradients_at=None, gradients_from_iobj=None, pointsizes=None, pointcloud_opacity=1.):
    """Two separate panels"""

    print("Mayavi.")
    sys.stdout.flush()

    if pointsizes is None:
        pointsizes = [0.2]*10

    if type(opacity) is list:
        opacities = opacity  # 1.0
    else:
        opacities = [opacity] + [0.2]*(len(vf_list)-1)  # 1.0, 0.2 #0.1

    for fi in range(len(vf_list)):
        if separate:
            mlab.figure()

        vf = vf_list[fi]
        vertex, faces = vf

        if vertex is None:

            continue
        if vertex.size == 0:
            print("Warning: empty vertex")
            continue
        if faces.size == 0:
            print("Warning: no faces")
            continue

        assert vertex.ndim == 2
        assert faces.ndim == 2
        assert vertex.shape == (vertex.shape[0], 3), str(vertex.shape)
        assert faces.shape == (faces.shape[0], 3), str(faces.shape)
        if type(mayavi_wireframe) is list:
            wire_frame1 = mayavi_wireframe[fi]
            assert len(mayavi_wireframe) == len(vf_list)
        else:
            wire_frame1 = mayavi_wireframe
        mlab.triangular_mesh([vert[0] for vert in vertex],
                         [vert[1] for vert in vertex],
                         [vert[2] for vert in vertex],
                         faces,
                         representation="surface" if not wire_frame1 else "wireframe",
                         opacity=opacities[fi], scale_factor=100.0)

        color_list = [(1, 0, 0), (0, 0, 0), (1, 1, 0), (0, 0, 1), (0, 1, 0)]
        i = 0
        for c in pointcloud_list:
            mlab.points3d(c[:, 0], c[:, 1], c[:, 2], color=color_list[i], scale_factor=pointsizes[i], opacity=pointcloud_opacity)
            i += 1
        del i

        if minmax is not None:
            (RANGE_MIN, RANGE_MAX) = minmax
            x = np.linspace(RANGE_MIN, RANGE_MAX, 2).reshape(2, 1)
            y = np.zeros((2, 1))
            z = np.zeros((2, 1))

            mlab.plot3d(x, y, z, line_width=3, name="x-axis")
            mlab.plot3d(y, x, z, line_width=3, name="y-axis")
            mlab.plot3d(z, y, x, line_width=3, name="z-axis")

            mlab.text3d(RANGE_MAX, 0, 0, "x", scale=0.3)
            mlab.text3d(0, RANGE_MAX, 0, "y", scale=0.3)
            mlab.text3d(0, 0, RANGE_MAX, "z", scale=0.3)
            mlab.text3d(RANGE_MIN, 0, 0, "-x", scale=0.3)
            mlab.text3d(0, RANGE_MIN, 0, "-y", scale=0.3)
            mlab.text3d(0, 0, RANGE_MIN, "-z", scale=0.3)

    mlab.show()
    return


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

        return 1.0 - np.sum(pv * pv, axis=1)

    def implicitGradient(self, pv):
        assert pv.ndim == 2
        assert pv.shape[1:] == (3,)
        grad = -2*pv

        return grad


def sphere_example(scale=1.):
    iobj = UnitSphere()

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

    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-3, +5, 0.2)

    gridvals = make_mc_values_grid(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE, old=False)
    vertex, faces = vtk_mc(gridvals, (RANGE_MIN, RANGE_MAX, STEPSIZE))

    print vertex, faces
    print("MC calculated.")
    sys.stdout.flush()
    display_simple_using_mayavi_2([(vertex, faces), (vertex, faces), ],
       pointcloud_list=[],
       mayavi_wireframe=[False, True, ], opacity=[1, 1, 0.9], gradients_at=None, separate=False, gradients_from_iobj=None,
       minmax=(RANGE_MIN, RANGE_MAX))

    exit()


if __name__ == '__main__':

    demo_sphere()
