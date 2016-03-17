""" Visualise and test an implementation of Ohtake & Belyaev 2002
"""

import sys
from timeit import default_timer as dtimer
import numpy as np

use_mayavi = True
if use_mayavi:
    from mayavi import mlab
else:
    # import matplotlib
    # #import matplotlib
    # #matplotlib.use('Qt4Agg')
    # #matplotlib.use('GTKAgg')  # will require pygtk
    # import matplotlib.pyplot as plt
    # from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    # from mpl_toolkits.mplot3d import Axes3D
    pass


from skimage import measure



from basic_types import normalize_vector4_vectorized
from basic_types import is_python3,check_vector4_vectorized
import example_objects
import mc_utils
import numerical_utils

#RANGE_MIN = -4; RANGE_MAX = 4; STEPSIZE = 0.15
#RANGE_MIN = -1; RANGE_MAX = 1; STEPSIZE = 0.15/2.0
#for french fries:
#RANGE_MIN = -1; RANGE_MAX = 2; STEPSIZE = 0.1
#RANGE_MIN = -1; RANGE_MAX = 2+1; STEPSIZE = 0.1*3*2

#RANGE_MIN = -4; RANGE_MAX = 4; STEPSIZE = 0.5
#RANGE_MIN = -2; RANGE_MAX = 2; STEPSIZE = 0.3

#RANGE_MIN = -2; RANGE_MAX = 4; STEPSIZE = 0.2 #blend
#RANGE_MIN = -2; RANGE_MAX = 4; STEPSIZE = 0.4 #blend low res

#RANGE_MIN = -3; RANGE_MAX = 3; STEPSIZE = 0.2 #bowl+holes
#exname = "bowl_15_holes"

#this will have false sharp edges
(RANGE_MIN,RANGE_MAX, STEPSIZE) = (-2, +4, 0.4)
exname = "blend_example2_discs"
rootfinding_lambda = 0.5 #default


(RANGE_MIN,RANGE_MAX, STEPSIZE) = (-2, +4, 0.4/2)
exname = "blend_example2_discs"
rootfinding_lambda = 0.5/3 #default: 0.5

#(RANGE_MIN,RANGE_MAX, STEPSIZE) = (-1, +2, 0.1)
#exname = "french_fries_vectorized"

#(RANGE_MIN,RANGE_MAX, STEPSIZE) = (-1.3, +1.3, 0.1*2*2)
#exname = "rdice_vec"

#DOES NOT WORK
(RANGE_MIN,RANGE_MAX, STEPSIZE) = (-1.5, +1.5, 0.3/3)
exname = "rdice_vec"
rootfinding_lambda = 0.5/2.0*2.0

raw_mc_plot = False
if raw_mc_plot:
    #bad
    (RANGE_MIN,RANGE_MAX, STEPSIZE) = (-2.5, +2.5, 0.3/3)
    exname = "screw3"
    rootfinding_lambda = 0.5
    do_qem = False
    do_vertex_resampleing = False
    mayavi_plot_centroids = False
    #qem_breakdown irrelevant
    do_project_centroids = False
else:
    #good
    (RANGE_MIN,RANGE_MAX, STEPSIZE) = (-2.5, +2.5, 0.3/3)
    exname = "screw3"
    rootfinding_lambda = 0.1 * 5
    do_qem = True
    do_vertex_resampleing = False
    mayavi_plot_centroids = False
    qem_breakdown = False
    do_project_centroids = True #crucial

#(RANGE_MIN,RANGE_MAX, STEPSIZE) = (-2.5, +2.5, 0.3/3)
#exname = "screw3"
#rootfinding_lambda = 0.5

#(RANGE_MIN,RANGE_MAX, STEPSIZE) = (-2.5, +2.5, 0.3/3)
#exname = "french_fries_vectorized"
#rootfinding_lambda = 0.5/3


#RANGE_MIN = -3.3; RANGE_MAX = 3.3; STEPSIZE = 0.2 #bowl+holes
#exname = "bowl_15_holes"
#rootfinding_lambda = 0.2 #*2*2  #large values -> moves some of them to wrong positions
##infinite

rng = np.arange(RANGE_MIN, RANGE_MAX, STEPSIZE)

plot_normals = False
plot_centroids = False
plot_ranks = False
plot_rank3 = False
plot_random_interior_points = False
plot_raycasts = False

#do_project_centroids = False
#do_qem = True
#qem_breakdown = True

mayavi_wireframe = False


""" Choose the object """
#exname = "bowl_15_holes"  # "blend_example2_discs" "french_fries_vectorized" "cube_example"
#exname = "blend_example2_discs" # 
#exname ="ell_example1" #
#exname = "first_csg"
#exname = "bowl_15_holes"
iobj = example_objects.make_example_vectorized(exname)
#iobj = example_objects.make_example_nonvec(exname)


#import vectorized
#iobj = cube1(vectorized)

print("Starting evaluation of implicit on the Grid.")
sys.stdout.flush()
t1s = dtimer()
vgrid = mc_utils.make_grid(iobj, rng, old=True)
#vgrid = mc_utils.make_grid_pointwise(iobj, rng)
assert vgrid.shape == (len(rng), len(rng), len(rng))
t1 = dtimer() - t1s
print('done grid')
sys.stdout.flush()


""" Use the marching cubes: Usually very fast """
t2s = dtimer()
verts, faces = measure.marching_cubes(vgrid, 0)
#print(verts)
#verts = (verts) * STEPSIZE + RANGE_MIN
#verts = (verts) * (rng[1]-rng[0]) + rng[0]
verts = ( (verts) * STEPSIZE + rng[0] )
verts = np.concatenate(( verts[:,1,np.newaxis], verts[:,0,np.newaxis],verts[:,2,np.newaxis] ) , axis=1)
t2 = dtimer() - t2s
print("Marching cubes done.")
sys.stdout.flush()
print("Timings: evaluation: %f " % (t1,)+" and MC: %f" % (t2,))
print(verts[0:5,:])
#exit()

def test_gradients(iobj, points, grads):
    n = points.shape[0]
    for i in range(n):
        v14 = points[i,:][np.newaxis,:]
        v4 = points[i,:]
        g = iobj.implicitGradient(v14)
        g1 = grads[i,:][np.newaxis,:]
        g2 = numerical_utils.numerical_gradient(iobj, v4, is_vectorized=True)
        THRESHOLD = 0.0001
        assert np.sum(np.abs(g2 - g1)) < THRESHOLD

def display_simple(verts, faces):
    fig = plt.figure(figsize=(10, 12))
    ax = fig.add_subplot(111, projection='3d')
    if is_python3():
        ax.axis('square')
    mesh = Poly3DCollection(verts[faces], alpha=0.2)
    mesh.set_facecolor([1, 0.5, 0.5])
    mesh.set_linewidth(0.2)
    mesh.set_antialiased(True)
    ax.add_collection3d(mesh)
    ax.set_xlim(np.min(rng), np.max(rng))
    ax.set_ylim(np.min(rng), np.max(rng))
    ax.set_zlim(np.min(rng), np.max(rng))
    ax.set_xlabel("x-axis")
    ax.set_ylabel("y-axis")
    ax.set_zlabel("z-axis")
    plt.show()

def display_simple_using_mayavi(verts, faces, centroids):
    # fig = mlab.figure()
    mlab.triangular_mesh(
        [vert[0] for vert in verts],
        [vert[1] for vert in verts],
        [vert[2] for vert in verts],
        faces,
        representation="surface" if not mayavi_wireframe else "wireframe",
        opacity=1,scale_factor = 100.0)
    #fancymesh

    if mayavi_plot_centroids:
        mlab.points3d(centroids[:,0], centroids[:,1], centroids[:,2], color=(1,0,0), scale_factor=0.01)

    x = np.linspace(RANGE_MIN,RANGE_MAX,2).reshape(2,1)
    y = np.zeros((2,1))
    z = np.zeros((2,1))

    mlab.plot3d(x,y,z,line_width=3,name="x-axis")
    mlab.plot3d(y,x,z,line_width=3,name="y-axis")
    mlab.plot3d(z,y,x,line_width=3,name="z-axis")

    mlab.show() #figure=fig,


def display_using_matplotlib(verts, faces, m, cloud1, new_cloud2, iobj_=None):
    """ Displays vertices, centroids (optional) and normal vectors """
    # Display resulting triangular mesh using Matplotlib.
    # Display the centroids
    # Display the normal vectors at centroids
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


    def add_random_interior_points(ax, iobj):
        """ Adding random points """
        n=1000
        import basic_types
        x = basic_types.make_random_vector_vectorized(n, 2, 1, type="rand", normalize=False)
        v = iobj.implicitFunction(x)
        x_sel =  x[ v >= 0 , :]
        print(".shape: ",x_sel.shape)
        ax.scatter(x_sel[:,0], x_sel[:,1], x_sel[:,2], c='b', marker='.')

    if plot_random_interior_points:
        add_random_interior_points(ax, iobj)


    #centroids = mesh_utils.centroids(verts, faces) #* 1.04
    centroids = m.centroids
    #centroids= verts;
    if plot_centroids:
        ax.scatter(centroids[:,0], centroids[:,1], centroids[:,2], c='r', marker='o')



    #if not cloud1 is None:
    #    ax.scatter(cloud1[:,0], cloud1[:,1], cloud1[:,2], c='b', marker='.')
    #if not new_cloud2 is None:
    #    ax.scatter(new_cloud2[:,0], new_cloud2[:,1], new_cloud2[:,2], c='r', marker='o')

    if plot_raycasts:
        import basic_types
        npoints = 1000
        radius = 5
        x = basic_types.make_random_vector_vectorized(npoints, radius, 1, type="rand", normalize=False)
        ray_x = x
        ray_x[:,0] += 4
        ray_x[:,1] += 10
        ray_n = -x # iobj_.implicitGradient(x)
        ray_n[:,3] = 1
        #, rayscast=None
        xx = numerical_utils.numerical_raycast_bisection_vectorized(iobj_, ray_x, ray_n)
        #xx = rayscast
        ax.scatter(xx[:,0], xx[:,1], xx[:,2], c='k', marker='.')

    if plot_ranks:
        print(m.verts_ranks.shape)
        xx = verts[m.verts_ranks==2,:]
        ax.scatter(xx[:,0], xx[:,1], xx[:,2], c='g', marker='o')
        xx = verts[m.verts_ranks==1,:]
        ax.scatter(xx[:,0], xx[:,1], xx[:,2], c='r', marker='o')
    if plot_rank3:
        xx = verts[m.verts_ranks==3,:]
        ax.scatter(xx[:,0], xx[:,1], xx[:,2], c='k', marker='o')

    pos = m.centroids
    pnormals = iobj.implicitGradient(pos) #- 10  #m.centroid_gradients
    if False:
        #slow:
        test_gradients(iobj, pos, pnormals)

    ons = np.ones((verts.shape[0],1))
    pos = np.concatenate( (verts, ons), axis=1) #verts nx3

    #import basic_types
    #pos = basic_types.make_random_vector_vectorized(100, 0.2,1)
    #pos[:,1] += 1
    #pos[:,0] += 2

    ival = iobj.implicitFunction(pos)
    print(np.mean(np.abs(ival)))


    #ax.scatter(pos[:,0], pos[:,1], pos[:,2], c='r', marker='o')

    pnormals = - iobj.implicitGradient(pos) #- 10  #m.centroid_gradients


    pnormals = normalize_vector4_vectorized(pnormals) 
    lm = STEPSIZE  # 0.2 * 2

    xyz = pos
    check_vector4_vectorized(xyz)
    uvw = pnormals

    if plot_normals:
        xx,yy,zz = xyz[:,0], xyz[:,1], xyz[:,2]
        uu,vv,ww = uvw[:,0], uvw[:,1], uvw[:,2]
        ax.quiver( xx,yy,zz,   uu,vv,ww,  length=np.abs(lm), arrow_length_ratio=0.3, alpha=0.3, pivot="tail")
        #arrow_length_ratio=   length=np.abs(lm)
        #pivot: tail | middle | tip




    ax.set_xlabel("x-axis")
    ax.set_ylabel("y-axis")
    ax.set_zlabel("z-axis")

    ax.set_xlim(np.min(rng), np.max(rng))
    ax.set_ylim(np.min(rng), np.max(rng))
    ax.set_zlim(np.min(rng), np.max(rng))

    plt.show()



import mesh_utils

print("faces,verts: ", faces.shape[0], verts.shape[0] )
mesh_utils.make_neighbour_faces_of_vertex(faces)
if False:
    mesh_utils.mesh_invariant(faces)
    mesh_utils.make_dual(verts, faces)

if False:
    #print(mesh_utils.centroids(verts, faces))
    mesh_utils.test_make_edge_lookup()

if False:
    (edges_of_faces, faces_of_edges, verts_of_edges) = mesh_utils.make_edge_lookup(faces)
    #Axes3D.quiver(*args, **kwargs)


import mesh1

if True:
    m = mesh1.Mesh_1(faces, verts)
    m.build_centroids()
    m.build_neighbours()
    m.evaluate_centroid_gradients(iobj)

    if do_vertex_resampleing:
        m.vertex_resampling()

    print("================================", do_qem)
    if do_qem:
        projection_maxdist = np.infty  # STEPSIZE

        if do_project_centroids:
            #no effect anyway
            m.update_centroids_and_gradients(iobj, lambda_=rootfinding_lambda, maxdist=projection_maxdist,
                method="v0.1_inplace")
            #second time does nothing
            #m.update_centroids_and_gradients(iobj,lambda_=rootfinding_lambda, maxdist=projection_maxdist,
            #    method="v0.1_inplace")
        else:
            pass

        if not qem_breakdown:
            m.quadratic_optimise_vertices(1)
            m.verts = m.new_verts
        else:
            m.quadratic_optimise_vertices(0.2)
            m.verts = m.new_verts
            for i in range(15):
                m.quadratic_optimise_vertices(0.4)
                m.verts = m.new_verts
    else:
        m.new_verts = m.verts

    # #display_using_matplotlib(verts, faces, m, None, None )
    # display_using_matplotlib(verts, faces, m, m.verts, m.new_verts )
    # #display_using_matplotlib(m.new_verts, faces, m, m.verts, m.new_verts )
    #display_using_matplotlib(m.verts, faces, m, m.verts, m.new_verts )

    #display_using_matplotlib(m.verts, faces, m, m.verts, m.new_verts, iobj_=iobj )

if use_mayavi:
    display_simple_using_mayavi(m.verts, faces, m.centroids)
else:
    display_simple(m.verts, faces)

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
