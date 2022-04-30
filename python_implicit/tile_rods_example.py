#tile rods
import sys
sys.path.append("vec4")

from easy_visualise_2 import *


import numpy as np
#global STEPSIZE
from tile1 import Tile1D

from example_objects import make_example_vectorized
ifunc = make_example_vectorized("cube_example")
ifunc = Tile1D(ifunc, start=np.array([1., 1., 1.])*(1.), direction=np.array([-1, -1, -1])*1., tilecount=4)
(RANGE_MIN, RANGE_MAX, STEPSIZE) = (-10., +10., 0.2*2)



#ifunc = SimpleCylinder(self, A, w, u, radius_u, radius_v, c_len)
#ifunc = Tile1D(ifunc, start=np.array([1., 1., 1.])*(1.), direction=np.array([-1, -1, -1])*1., tilecount=4)
#(RANGE_MIN, RANGE_MAX, STEPSIZE) = (-10., +10., 0.2*2)

def rods_row(N=6,rod_len = 10.):
    from basic_types import make_vector4
    from vectorized import SimpleCylinder
    rod_r = 1.
    cage_r = 10.
    delta_d = 0.4

    #return cage_rods(rod_r, rod_len, cage_r, N), (-32 / 2, +32 / 2, 1.92 / 4.0)
    #x, y = cage_r * math.sin(th), cage_r * math.cos(th)
    #A = make_vector4(x, y, -rod_len / 2.)
    A = make_vector4(-rod_r*2*(N/2), 0, -rod_len/2.)
    w = make_vector4(0, 0, 1)
    u = make_vector4(1, 0, 0)
    rod = SimpleCylinder(A, w, u, rod_r, rod_r, rod_len)
    ifunc = Tile1D(rod, start=A[:3]-np.array([1., 0., 0.])*rod_r, direction=np.array([1., 0, 0])*(2*rod_r+delta_d), tilecount=N)
    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-10., +10., 0.2*2)
    return ifunc, (RANGE_MIN, RANGE_MAX, STEPSIZE)

ifunc, (RANGE_MIN, RANGE_MAX, STEPSIZE) = rods_row(N=6)

from twist_z import TwistZ
ifunc, (RANGE_MIN, RANGE_MAX, STEPSIZE) = rods_row(N=6, rod_len=10.)
ifunc = TwistZ(ifunc, 0.25/10.)
(RANGE_MIN, RANGE_MAX, STEPSIZE) = (-10., +10., 0.2*2)

#todo: negative of this shape which goes inside it.

#todo: different lengths

#todo: length-preserving twist.

#coool!
from twist_z import TwistZ
ifunc, (RANGE_MIN, RANGE_MAX, STEPSIZE) = rods_row(N=6, rod_len=30.)
ifunc = TwistZ(ifunc, 0.25/10.)
(RANGE_MIN, RANGE_MAX, STEPSIZE) = (-16., +16., 0.2*2)


"""
#nice tall:
from twist_z import TwistZ
ifunc, (RANGE_MIN, RANGE_MAX, STEPSIZE) = rods_row(N=3, rod_len=30.)
ifunc = TwistZ(ifunc, 0.25/10.)
(RANGE_MIN, RANGE_MAX, STEPSIZE) = (-16., +16., 0.2*2)

tr = Transformed(ifunc);
tr.move(center[0], center[1], center[2])
tr.rotate(-20, along=make_vector3(1, 1, 1), units="deg") .move(0.2 * scale, 0, 0)
ifunc = tr
(RANGE_MIN, RANGE_MAX, STEPSIZE) = (-16., +16., 0.2*2)
"""


from twist_z import TwistZ
ifunc, (RANGE_MIN, RANGE_MAX, STEPSIZE) = rods_row(N=3, rod_len=30.)
ifunc = TwistZ(ifunc, 0.25/10.)
d =np.array([0., 1, 0])*8.
ifunc = Tile1D(ifunc, start=d*(-0.5), direction=d, tilecount=2+1, back_count=1)
(RANGE_MIN, RANGE_MAX, STEPSIZE) = (-16., +16., 0.2*2)


# def cage_rods(rod_r, rod_len, cage_r, N):
#     import math
#     un = None
#     for i in range(N):
#         th = i / float(N) * np.pi * 2
#         x, y = cage_r * math.sin(th), cage_r * math.cos(th)
#         A = make_vector4(x, y, -rod_len / 2.)
#         w = make_vector4(0, 0, 1)
#         u = make_vector4(1, 0, 0)
#         c = SimpleCylinder(A, w, u, rod_r, rod_r, rod_len)

#         if un is None:
#             un = c
#         else:
#             un = vectorized.CrispUnion(un, c)
#     return un

#print ifunc


from stl_tests import make_mc_values_grid
from vtk_mc import vtk_mc

#ifunc, (RANGE_MIN, RANGE_MAX, STEPSIZE) = getifunc()

from stl_tests import make_mc_values_grid
gridvals = make_mc_values_grid(ifunc, RANGE_MIN, RANGE_MAX, STEPSIZE, old=False)
print np.sum(gridvals.ravel()<0)
#exit()
verts, facets = vtk_mc(gridvals, (RANGE_MIN, RANGE_MAX, STEPSIZE))
#return verts, facets

#print verts
#print facets
#exit()

display_simple_using_mayavi_2([(verts, facets),], pointcloud_list=[], minmax=(RANGE_MIN, RANGE_MAX) )
        #mayavi_wireframe=[False], opacity=[1.0],
        #separate_panels=True, gradients_at=None, gradients_from_iobj=None, pointsizes=None, pointcloud_opacity=1.,
        #add_noise=[], noise_added_before_broadcast=False):
