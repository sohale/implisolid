""" Visualises based on matplotlib and MC from skimage:
Vectorized version.
"""

"""
Speedup 50x ! :
Non-vectorized version:  eval: 17.838456  and MC: 0.057518
    Vectorized version:  eval: 0.351882  and MC: 0.063542
"""


import sys
from timeit import default_timer as dtimer
import numpy as np

import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import Axes3D
from skimage import measure

#from basic_types import normalize_vector4_vectorized
from basic_types import is_python3
import example_objects
import mc_utils


#RANGE_MIN = -4; RANGE_MAX = 4; STEPSIZE = 0.15
#RANGE_MIN = -1; RANGE_MAX = 1; STEPSIZE = 0.15/2.0
#for french fries:
#RANGE_MIN = -1; RANGE_MAX = 2; STEPSIZE = 0.1
#RANGE_MIN = -1; RANGE_MAX = 2+1; STEPSIZE = 0.1*3

(RANGE_MIN,RANGE_MAX, STEPSIZE) = (-1, +3, 0.1*3)
exname = "blend_example2_discs"

(RANGE_MIN,RANGE_MAX, STEPSIZE) = (-1, +3, 0.1)
exname = "screw1"

rng = np.arange(RANGE_MIN, RANGE_MAX, STEPSIZE)

plot_centroids = False

""" Choose the object """
#exname = "bowl_15_holes"  # "blend_example2_discs" "french_fries_vectorized" "cube_example"
#exname = "blend_example2_discs"
iobj = example_objects.make_example_vectorized(exname)

print("Starting evaluation of implicit on the Grid.")
sys.stdout.flush()
t1s = dtimer()
vgrid = mc_utils.make_grid(iobj, rng, old=True)
assert vgrid.shape == (len(rng), len(rng), len(rng))
t1 = dtimer() - t1s
print('done grid')
sys.stdout.flush()


""" Use the marching cubes: Usually very fast """
t2s = dtimer()
verts, faces = measure.marching_cubes(vgrid, 0)
verts = (verts) * STEPSIZE + RANGE_MIN
t2 = dtimer() - t2s
print("Marching cubes done.")
sys.stdout.flush()
print("Timings: evaluation: %f " % (t1,)+" and MC: %f" % (t2,))

import numerical_utils

t3s = dtimer()
radius_of_max_change = STEPSIZE/2.0
print(verts.shape)
print(numerical_utils.error_i_vectorized(verts, iobj, 'abs') )
#verts2 = numerical_utils.optimize_vertices(verts2, iobj, radius_of_max_change)
#print(numerical_utils.error_i_vectorized(verts, iobj, 'abs') )

t3 = dtimer() - t3s
print("Finished correction.")
sys.stdout.flush()


def display_using_matplotlib(verts, faces):
    # Display resulting triangular mesh using Matplotlib. This can also be done
    # with mayavi (see skimage.measure.marching_cubes docstring).
    fig = plt.figure(figsize=(10, 12))
    ax = fig.add_subplot(111, projection='3d')
    if is_python3():
        ax.axis('square')

    # Fancy indexing: `verts[faces]` to generate a collection of triangles
    mesh = Poly3DCollection(verts[faces], alpha=0.2)
    mesh.set_facecolor([1, 0.5, 0.5])
    mesh.set_linewidth(0.2)
    mesh.set_antialiased(True)

    ax.add_collection3d(mesh)

    import mesh_utils  # used for displaying the centroids

    centroids = mesh_utils.centroids(verts, faces) #* 1.04
    #centroids= verts;
    if plot_centroids:
        ax.scatter(centroids[:,0], centroids[:,1], centroids[:,2], c='r', marker='o')

    ax.set_xlabel("x-axis")
    ax.set_ylabel("y-axis")
    ax.set_zlabel("z-axis")

    ax.set_xlim(np.min(rng), np.max(rng))
    ax.set_ylim(np.min(rng), np.max(rng))
    ax.set_zlim(np.min(rng), np.max(rng))

    plt.show()


display_using_matplotlib(verts, faces)


def plot_using_isosurface():
    """ not tested"""
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.gca(projection='3d')


    z = calc_iso_surface( my_array, my_value=0.0, zs=zs, interp_order=6 )
    ax.plot_surface( xx, yy, vgrid, cstride=4, rstride=4, color='b')

    plt.ion()
    plt.show()
