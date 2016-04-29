""" Utilities to facilitate using the Marching Cubes algorithm.
The MC algorithm itself is not called in this module. """

import numpy as np
import sys

VERBOSE = True


class PolygonizationError(Exception):
    pass


def _prepare_grid(rng):
    assert rng.size < 200
    if rng.size > 200:
        raise PolygonizationError(("Grid too large ( >200 ): ", rng.size))

    (yy, xx, zz) = np.meshgrid(rng, rng, rng)
    xyza = np.transpose(np.vstack([xx.ravel(), yy.ravel(), zz.ravel()]))
    assert xyza.shape[1:] == (3,)

    if VERBOSE:
        print(xyza.shape)
        print("done alloc")
        sys.stdout.flush()

    return xyza

@profile
def make_grid(iobj, rng, old=None):
    assert old is not None
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
