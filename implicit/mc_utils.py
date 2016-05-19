""" Utilities to facilitate using the Marching Cubes algorithm.
The MC algorithm itself is not called in this module. """
import profile_support
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
    xyz_nparray = np.transpose(np.vstack([xx.ravel(), yy.ravel(), zz.ravel(), (xx*0+1).ravel()]))
    assert xyz_nparray.shape[1:] == (4,)

    if VERBOSE:
        print(xyz_nparray.shape)
        print("done alloc")
        sys.stdout.flush()

    return xyz_nparray


def _prepare_grid_old(rng):
    assert rng.size < 200
    if rng.size > 200:
        raise PolygonizationError(("Grid too large ( >200 ): ", rng.size))

    (xx, yy, zz) = np.meshgrid(rng, rng, rng)
    #xyz = np.mgrid( rng, rng, rng )
    #xyz_nparray = xyz.reshape((len(rng)**3, 3))
    #xyz_nparray = np.concat( ( np.expand_dims( xyz_nparray, axis=3 ), ones(len(rng)**3,1)  ), axis=3 )
    #assert xyz_nparray.shape == (len(rng), len(rng), len(rng), 4)

    #X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j] #??

    xyz_nparray = np.transpose(np.vstack([xx.ravel(), yy.ravel(), zz.ravel(), (xx*0+1).ravel()]))
    assert xyz_nparray.shape[1:] == (4,)

    if VERBOSE:
        print(xyz_nparray.shape)
        print("done alloc")
        sys.stdout.flush()

    return xyz_nparray


def slow_grid__dont_use():
    """ #Slow way:
    """
    vgrid = xx
    for i in range( len(rng) ):
        for j in range( len(rng) ):
            for k in range( len(rng) ):
                x_ = make_vector4(rng[i], rng[j], rng[k] )
                #vgrid[i,j,k] = i>1 and j>1 and k>1  and i<len(rng)-1 and j<len(rng)-1 and k<len(rng)-1  #iobj.implicitFunction( x )
                vgrid[i,j,k] = iobj.implicitFunction( x_ )

#@profile
def make_grid(iobj, rng, old, return_xyz=False):
    assert old is not None
    if old:
        xyz_nparray = _prepare_grid_old(rng)
    else:
        xyz_nparray = _prepare_grid(rng)
    #slow_grid__dont_use()

    vgrid_v = iobj.implicitFunction(xyz_nparray)
    vgrid = np.reshape(vgrid_v, (len(rng), len(rng), len(rng)), order='C')

    if np.sum(vgrid_v > 0) == 0:
        raise PolygonizationError("The shape is empty. No interior points detected")
    if np.sum(vgrid_v < 0) == 0:
        raise PolygonizationError("The solid volume fills all the space. No exterior points detected")
    if VERBOSE:
        print("interior points:", np.sum(vgrid_v > 0))

    if return_xyz:
        return vgrid, xyz_nparray
    else:
        return vgrid


from basic_types import make_vector4

def make_grid_pointwise(iobj, rng):
    """ An inefficient implementation kept for histotrical reasons"""
    assert rng.size < 200
    (xx, yy, zz) = np.meshgrid(rng, rng, rng)

    vgrid = xx

    aany = False
    for i in range(len(rng)):
        for j in range(len(rng)):
            for k in range(len(rng)):
                x = make_vector4(rng[i], rng[j], rng[k])
                vgrid[i, j, k] = iobj.implicitFunction(x)
                if vgrid[i, j, k] > 0:
                    aany = True

    assert aany, "No point detected"
    return vgrid
