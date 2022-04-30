""" Visualises based on matplotlib and MC from skimage"""
import sys
from timeit import default_timer as dtimer
import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from skimage import measure

from basic_types import check_vector4
from basic_types import is_python3
import mc_utils




def display_uting_matplotlib(verts, faces):
    # Display resulting triangular mesh using Matplotlib. This can also be done
    # with mayavi (see skimage.measure.marching_cubes docstring).
    fig = plt.figure(figsize=(10, 12))
    ax = fig.add_subplot(111, projection='3d')
    if is_python3():
        ax.axis('square')

    # Fancy indexing: `verts[faces]` to generate a collection of triangles
    mesh = Poly3DCollection(verts[faces])
    mesh.set_facecolor([1, 0.5, 0.5])
    mesh.set_linewidth(0.2)
    mesh.set_antialiased(True)

    ax.add_collection3d(mesh)

    ax.set_xlabel("x-axis")
    ax.set_ylabel("y-axis")
    ax.set_zlabel("z-axis")
    ax.set_xlim(RANGE_MIN, RANGE_MAX)
    ax.set_ylim(RANGE_MIN, RANGE_MAX)
    ax.set_zlim(RANGE_MIN, RANGE_MAX)

    plt.show()

    #alternative
    #import matplotlib.pyplot as plt
    #from mpl_toolkits.mplot3d import Axes3D
    #fig = plt.figure()
    #ax = fig.gca(projection='3d')
    ##z = calc_iso_surface( my_array, my_value=0.0, zs=zs, interp_order=6 )
    ##ax.plot_surface( xx, yy, vgrid, cstride=4, rstride=4, color='b')
    #plt.ion()
    #plt.show()

def display_two_meshes_matplotlib(verts1, faces1, verts2, faces2):
    # Display resulting triangular mesh using Matplotlib. This can also be done
    # with mayavi (see skimage.measure.marching_cubes docstring).
    fig = plt.figure(figsize=(10, 12))
    ax = fig.add_subplot(1,2,1, projection='3d')
    if is_python3():
        ax.axis('square')

    # Fancy indexing: `verts[faces]` to generate a collection of triangles
    mesh = Poly3DCollection(verts1[faces1])
    mesh.set_facecolor([1, 0.5, 0.5])
    mesh.set_linewidth(0.2)
    mesh.set_antialiased(True)

    ax.add_collection3d(mesh)

    ax.set_xlabel("x-axis")
    ax.set_ylabel("y-axis")
    ax.set_zlabel("z-axis")
    ax.set_xlim(RANGE_MIN, RANGE_MAX)
    ax.set_ylim(RANGE_MIN, RANGE_MAX)
    ax.set_zlim(RANGE_MIN, RANGE_MAX)

    ax2 = fig.add_subplot(1,2,2, projection='3d')
    if is_python3():
        ax.axis('square')

    # Fancy indexing: `verts[faces]` to generate a collection of triangles
    mesh2 = Poly3DCollection(verts2[faces2])
    mesh2.set_facecolor([1, 0.5, 0.5])
    mesh2.set_linewidth(0.2)
    mesh2.set_antialiased(True)

    ax2.add_collection3d(mesh2)

    ax2.set_xlabel("x-axis")
    ax2.set_ylabel("y-axis")
    ax2.set_zlabel("z-axis")
    ax2.set_xlim(RANGE_MIN, RANGE_MAX)
    ax2.set_ylim(RANGE_MIN, RANGE_MAX)
    ax2.set_zlim(RANGE_MIN, RANGE_MAX)

    plt.show()

import example_objects

#example_nbame = 'csg_example1'
#example_nbame = 'bowl_hole'
#example_nbame = 'rdice'
example_nbame = "french_fries"

iobj = example_objects.make_example_nonvec(example_nbame)

""" Prepare the grid data arrays """
#STEPSIZE = 0.15
#RANGE_MAX = 4
#RANGE_MIN = -4

#for dice only
STEPSIZE = 0.15/2.0*2.0
RANGE_MAX = 3
RANGE_MIN = -1
rng = np.arange(RANGE_MIN, RANGE_MAX, STEPSIZE)

print("Starting evaluation of implicit on the Grid")
sys.stdout.flush()
t1s = dtimer()

vgrid = mc_utils.make_grid_pointwise(iobj, rng)

t1 = dtimer() - t1s
print("done grid")
sys.stdout.flush()
if np.sum(np.ravel(vgrid) > 0) == 0:
    print("No point detected")
    raise Error("No point detected")


# Use marching cubes to obtain the surface mesh of these ellipsoids
t2s = dtimer()

verts, faces = measure.marching_cubes(vgrid, 0)
verts = (verts) * STEPSIZE + RANGE_MIN

t2 = dtimer() - t2s
print("Marching cubes done.")
sys.stdout.flush()

print("Timings: evaluation: %f " % (t1,)+" and MC: %f" % (t2,))
sys.stdout.flush()



#face_area(verts, faces)

#display_uting_matplotlib(verts, faces)

verts2 = verts.copy()
faces2 = faces.copy()

#np.save("verts", verts)
#np.save("faces", faces)

print("starting correction.")
sys.stdout.flush()

import numerical_utils

t3 = 0
if True:
    t3s = dtimer()

    radius_of_max_change = STEPSIZE/2.0
    print(numerical_utils.error_i(verts, iobj, 'abs') )
    verts2 = numerical_utils.optimize_vertices(verts2, iobj, radius_of_max_change)
    print(numerical_utils.error_i(verts, iobj, 'abs') )

    t3 = dtimer() - t3s
    print("Finished correction.")
    sys.stdout.flush()

print("Timings: evaluation: %f " % (t1,)+", MC: %f" % (t2,)+", Z-flow: %f" %(t3,))
sys.stdout.flush()

display_two_meshes_matplotlib(verts, faces, verts2, faces2)
