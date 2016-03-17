""" Visualises based on Mayavi and MC from skimage.
"""

import sys
from timeit import default_timer as dtimer
import numpy as np

# Set toolkit to wx if you get an error with default toolkit (qt4)
# from traits.etsconfig.api import ETSConfig
# ETSConfig.toolkit = 'wx'

from mayavi import mlab
from skimage import measure

#from basic_types import normalize_vector4_vectorized
import mc_utils
import example_objects




#120 2
RANGE_MIN = -4; RANGE_MAX = 4; STEPSIZE = 0.1
rng = np.arange(RANGE_MIN, RANGE_MAX, STEPSIZE)



""" Choose the object """
iobj = example_objects.make_example_vectorized("bowl_15_holes")

print("Starting evaluation of implicit on the Grid.")
sys.stdout.flush()
t1s = dtimer()
vgrid = mc_utils.make_grid(iobj, rng, old=True)
assert vgrid.shape == (len(rng),len(rng),len(rng))
t1 = dtimer() - t1s
print('done grid')
sys.stdout.flush()


""" Use the marching cubes: Usually very fast """
t2s = dtimer()
verts, faces = measure.marching_cubes(vgrid, 0)
verts = (verts) * STEPSIZE + RANGE_MIN
verts = np.concatenate(( verts[:,1,np.newaxis], verts[:,0,np.newaxis],verts[:,2,np.newaxis] ) , axis=1)



x = np.linspace(-10,10,2).reshape(2,1)
y = np.zeros((2,1))
z = np.zeros((2,1))

### Plotting with mayavi
fig = mlab.figure()
mlab.triangular_mesh([vert[0] for vert in verts],
                     [vert[1] for vert in verts],
                     [vert[2] for vert in verts],faces,representation="surface",opacity=1,figure=fig,scale_factor = 100.0)


# set the axes

centroids = np.mean( verts[faces[:],:], axis=1 )
centroids_hom = np.c_[centroids,np.ones((len(centroids),1))]

x_verts = [vert[0] for vert in centroids_hom]
y_verts = [vert[1] for vert in centroids_hom]
z_verts = [vert[2] for vert in centroids_hom]

print verts.shape
UVW = iobj.implicitGradient(centroids_hom)
UVW_normals =  - UVW / np.linalg.norm(UVW,axis = 1).reshape(len(centroids),1)[0]
#todo: use normalize function
assert(UVW_normals.shape[0] == len(centroids))
mlab.quiver3d(x_verts,y_verts,z_verts, UVW_normals[:,0],UVW_normals[:,1],UVW_normals[:,2],color=(0,0,0))
mlab.plot3d(x,y,z,line_width=3,name="x-axis")
mlab.plot3d(y,x,z,line_width=3,name="y-axis")
mlab.plot3d(z,y,x,line_width=3,name="z-axis")
mlab.outline()
# plot the normals


t2 = dtimer() - t2s
print("Marching cubes done.")
sys.stdout.flush()

import mesh_utils
#print(mesh_utils.centroids(verts, faces))

mlab.show()


print("Timings: evaluation: %f " % (t1,)+" and MC: %f" % (t2,))
