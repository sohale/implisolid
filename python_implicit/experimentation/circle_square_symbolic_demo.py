from mayavi import mlab
import numpy as np
import types

# x,y,z = np.mgrid[-3:3:20j,-3:3:20j,-3:3:20j]


def cube_circle_fnc(x,y,z):
    return np.maximum(1 - np.minimum(np.maximum(abs(x),abs(y)),abs(z)),1 - x**2 - (y-1)**2 - (z + 1)**2)

#
# mlab.figure()
# mlab.contour3d(x,y,z,cube_fnc,contours=[0])
# mlab.show()
#

def cube_fnc(x, y, z):
   return 1 - (x ** 20 + y ** 20 + z ** 20)**(1./20)

def torus_fnc(x, y, z, r=1):
    r = 2
    rx = 0.1
    ry = 0.1
    rz = 0.1
    return (1 - (r - ((x / rx) ** 20 + (y / ry) ** 20)**0.05 ) ** 20 - (z / rz) ** 20 )

def contour_fnc(dim,slices,fnc):
    """
    Creates a fast contour of fnc in 3d space
    ARGUMENTS

    dim and slices: arguments for np.mgrid
    fnc           : is the desired function, this should be a callable

    NOTES:

    This function assumes that mayavi.mlab is imported:
    if not:
        from mayavi import mlab
    if mayavi is not installed:
        install with pip :D
    """

    assert isinstance(fnc, types.FunctionType), "Argument fnc should be a callable"

    x, y, z = np.mgrid[-dim:dim:slices,-dim:dim:slices,-dim:dim:slices]
    mlab.figure()
    mlab.contour3d(x, y, z, fnc, contours=[0])
    mlab.show()


def main():
    # contour_fnc(1.2,150j,cube_fnc)
    a = np.random.random((100000,3))


if __name__ == "__main__":
    main()
