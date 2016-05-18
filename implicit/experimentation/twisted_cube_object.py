from example_objects import twisted_cube_example
from ipdb import set_trace
import numpy as np
from optimize_dual_mesh import MeshOptimizer
import numpy as np
import numexpr as ne
from sympy import lambdify, symbols
import sys
from ipdb import set_trace
from mayavi import mlab
from vtk_mc import vtk_mc

VERBOSE = True


def _prepare_grid(rng):
    assert rng.size < 200
    if rng.size > 200:
        raise PolygonizationError(("Grid too large ( >200 ): ", rng.size))

    (yy, xx, zz) = np.meshgrid(rng, rng, rng)
    xyza = np.transpose(np.vstack([xx.ravel(), yy.ravel(), zz.ravel(), (xx * 0 + 1).ravel()]))
    assert xyza.shape[1:] == (4,)

    if VERBOSE:
        print(xyza.shape)
        print("done alloc")
        sys.stdout.flush()

    return xyza


def make_grid(iobj, rng, old=None):
    # assert old is not None
    if old:
        xyza = _prepare_grid_old(rng)
    else:
        xyza = _prepare_grid(rng)
    #slow_grid__dont_use()

    vgrid_v = iobj.implicitFunction(xyza)
    vgrid = np.reshape(vgrid_v, (len(rng), len(rng), len(rng)), order='C')

    if np.sum(vgrid_v > 0) == 0:
        raise PolygonizationError("The shape is empty. No interior points detected")
    if VERBOSE:
        print("interior points:", np.sum(vgrid_v > 0))
    return vgrid


def make_mc_mesh_scikit(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE):
    """ Uses Scikit's MC algorithm,which has minor bugs. """
    rng = np.arange(RANGE_MIN, RANGE_MAX, STEPSIZE)
    import mc_utils
    vgrid = mc_utils.make_grid(iobj, rng, old=True)
    verts, faces = vtk_mc(vgrid, (RANGE_MIN, RANGE_MAX, STEPSIZE))
    verts = ((verts) * STEPSIZE + rng[0])
    print("OLD: swapping x,y")
    verts = np.concatenate((verts[:, 1, np.newaxis], verts[:, 0, np.newaxis], verts[:, 2, np.newaxis]), axis=1)
    return verts, faces

if __name__ == "__main__":
    iobj = twisted_cube_example()
    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-3, 3, 0.05)
    verts, faces = make_mc_mesh_scikit(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE)
    mlab.triangular_mesh([vert[0] for vert in verts],
                         [vert[1] for vert in verts],
                         [vert[2] for vert in verts],
                         faces,
                         representation="surface")
    mlab.show()
