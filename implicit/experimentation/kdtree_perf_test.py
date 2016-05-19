import sys
sys.path.append("..")

import numpy as np
from ipdb import set_trace

#kdtree_perf.py
from vtk_mc import vtk_mc
import mc_utils
from basic_types import check_vector4_vectorized


def augment4(x):
    v4 = np.concatenate((x, np.ones((x.shape[0], 1))), axis=1)
    check_vector4_vectorized(v4)
    return v4



objname = "bowl_15_holes"
    #"rcube_vec"
    #"rdice_vec"  #
    #"cube_example" # problem: zero facet areas.  otherwise, it works.
    #"ell_example1"  #+
    #"bowl_15_holes"  # works too. But too many faces => too slow, too much memory. 32K?
    #"french_fries_vectorized"

global STEPSIZE
from example_objects import make_example_vectorized
iobj = make_example_vectorized(objname)
(RANGE_MIN, RANGE_MAX, STEPSIZE) = (-3, +5, 0.2)
#STEPSIZE = STEPSIZE / 2.


rng = np.arange(RANGE_MIN, RANGE_MAX, STEPSIZE)
#vgrid = mc_utils.make_grid(iobj, rng, old=old)
gridvals, xyz = mc_utils.make_grid(iobj, rng, old=False, return_xyz=True)
v_xyz = gridvals.ravel()

#from stl_tests import make_mc_values_grid
#gridvals = make_mc_values_grid(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE, old="3")
verts, facets = vtk_mc(gridvals, (RANGE_MIN, RANGE_MAX, STEPSIZE))
print("MC calculated.");sys.stdout.flush()


if False:
    import visual5
    visual5. display_simple_using_mayavi_2(
        [(verts, facets), (verts, facets)],
        mayavi_wireframe=[False, True], opacity=[0.4, 0.3],
        #gradients_at=c3,
        separate_panels=False,
        #gradients_from_iobj=iobj,
        minmax=(RANGE_MIN, RANGE_MAX),
        #pointcloud_list=[],
        #add_noise=[0, 0], noise_added_before_broadcast=True,
        #pointcloud_list=[new_centroids[z12, :]], pointsizes=[0.02], #pointcloud_list=[point_collector.get_as_array()], pointsizes=[0.01],
        #pointcloud_list=[new_centroids[nzeros_c, :]], pointsizes=[0.02], #pointcloud_list=[point_collector.get_as_array()], pointsizes=[0.01],
        #labels=(new_centroids, z12), grad_arrow_len=0.2/2.)
        #labels=(new_centroids, nzeros_c), grad_arrow_len=average_edge*1.  #
    )


from utils import Timer

from sklearn.neighbors import KDTree

#v_from = verts[:, :3]
centroids3 = np.mean(verts[facets, :], axis=1)
#v_from = centroids3[:, :]

#negative_indices = np.nonzero(v_xyz < -0.0001)[0]
#v_from = xyz[negative_indices, :3]
negative_indices = np.nonzero(iobj.implicitFunction(augment4(centroids3)) < -0.0001)[0]
v_from = centroids3[negative_indices, :3]


positive_indices = np.nonzero(v_xyz > -0.0001)[0]
struc = xyz[positive_indices, :3]
k = 1
with Timer() as t1:
    kdt = KDTree(struc, leaf_size=30, metric='euclidean')
    dist3, ind3 = kdt.query(v_from, k=k, return_distance=True)
dist3, ind3 = dist3[:, k-1], ind3[:, k-1]
print "time:", t1.interval

v_to = struc[ind3, :]

print v_from.shape, v_to.shape

import visual5
visual5. display_simple_using_mayavi_2(
    [(verts, facets), (verts, facets)],
    mayavi_wireframe=[False, True], opacity=[0.2, 0.2],
    separate_panels=False,
    minmax=(RANGE_MIN, RANGE_MAX),
    fromto=(v_from[:10000, :], v_to[:10000, :]),
    #pointcloud_list=[],
    #add_noise=[0, 0], noise_added_before_broadcast=True,
    #pointcloud_list=[new_centroids[z12, :]], pointsizes=[0.02], #pointcloud_list=[point_collector.get_as_array()], pointsizes=[0.01],
    #pointcloud_list=[new_centroids[nzeros_c, :]], pointsizes=[0.02], #pointcloud_list=[point_collector.get_as_array()], pointsizes=[0.01],
    #labels=(new_centroids, z12), grad_arrow_len=0.2/2.)
    #labels=(new_centroids, nzeros_c), grad_arrow_len=average_edge*1.  #
)
