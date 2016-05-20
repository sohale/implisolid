""" Utilities to facilitate using the Marching Cubes algorithm.
The MC algorithm itself is not called in this module. """

import numpy as np
import sys


VERBOSE = True


class PolygonizationError(Exception):
    pass


def _prepare_grid(rng):
    """
    rng: like the output of np.arange()
    """
    assert rng.size < 200
    if rng.size > 200:
        raise PolygonizationError(("Grid too large ( >200 ): ", rng.size))

    (xx, yy, zz) = np.meshgrid(rng, rng, rng)
    xyz_nparray = np.transpose(np.vstack([xx.ravel(), yy.ravel(), zz.ravel()]))
    assert xyz_nparray.shape[1:] == (3,)

    return xyz_nparray


def make_grid(iobj, rng, old=None, return_xyz=False):
    xyz_nparray = _prepare_grid(rng)

    vgrid_v = iobj.implicitFunction(xyz_nparray)
    vgrid = np.reshape(vgrid_v, (len(rng), len(rng), len(rng)), order='C')

    if np.sum(vgrid_v > 0) == 0:
        raise PolygonizationError("The shape is empty. No interior points detected")
    if VERBOSE:
        print("interior points:", np.sum(vgrid_v > 0))
    if return_xyz:
        return vgrid, xyz_nparray
    else:
        return vgrid
