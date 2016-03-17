import numpy as np

from ohtake_surface_projection import set_centers_on_surface__ohtake
from ohtake_surface_projection import display_simple_using_mayavi_

def test_proj_ohtak2():

    from stl_tests import make_mc_mesh_scikit

    exname = "french_fries_vectorized"
    import example_objects
    iobj = example_objects.make_example_vectorized(exname, 8.0)

    #(RANGE_MIN, RANGE_MAX, STEPSIZE) = (-2.*8, +4.*8, 0.4*8/5)
    #(RANGE_MIN, RANGE_MAX, STEPSIZE) = (-16, +32, 0.64*3*2/2)
    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-2*8, +4*8, 0.64*3*2/2/4. *8*0.2 )  # / 8.   #0.768
    #STEPSIZE2 = STEPSIZE * 0.2
    STEPSIZE, STEPSIZE2 = (1.5, 0.384 )

    verts, faces = make_mc_mesh_scikit(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE)

    print verts.shape, faces.shape

    average_edge = STEPSIZE

    # verts = optimise_mesh(verts, faces, iobj)

    c3_0 = np.mean(verts[faces[:], :], axis=1)
    # add extra points
    offsets = [0]
    print c3_0.shape
    #offsets = [0,+0.1,-0.1]
    #c3 = np.concatenate((c3, c3+STEPSIZE*0.1, c3+STEPSIZE*(-0.1)), axis=0)
    c3 = np.zeros((0, 3))
    for x in offsets:
        c3 = np.concatenate((c3, c3_0+STEPSIZE*x, ), axis=0)
    #c3 = np.concatenate((c3,), axis=0)
    centroids = np.concatenate((c3, np.ones((c3.shape[0], 1))), axis=1)

    nones_map = centroids[:, 0]*0 > 100
    new_centroids = centroids.copy()

    set_centers_on_surface__ohtake(iobj, new_centroids, average_edge, nones_map)

    print
    #print np.arange(nones_map.shape[0])[np.logical_not(nones_map)].shape
    print np.sum(np.logical_not(nones_map)), "(success)  + ", np.sum(nones_map), " (failed)"
    new_centroids2 = new_centroids[np.logical_not(nones_map), :]

    #todo: Use "None" for non-converging ones. (and visualise them)

    if False:
        #(RANGE_MIN, RANGE_MAX, STEPSIZE) = (-16, +32, 0.64*3*2/2)
        #verts2, faces2 = make_mc_mesh_scikit(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE * 0.2/2.0) #15 million!
        verts2, faces2 = make_mc_mesh_scikit(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE2)
    verts2, faces2 = None,None


    #m = m2stl_mesh(verts, faces)
    #if ACTUALLY_SAVE:
    #    m.save('implicit6-blend.stl')  # wow

    centroids_failed = centroids[nones_map, :]  #green ones

    #display_simple_using_mayavi_([(verts, faces), (verts2, faces2)], [centroids, new_centroids2], minmax=(RANGE_MIN, RANGE_MAX))
    #display_simple_using_mayavi_([(verts, faces), (verts2, faces2)], [centroids, new_centroids2, centroids_failed], minmax=(RANGE_MIN, RANGE_MAX))
    #display_simple_using_mayavi_([(verts, faces), (verts2, faces2)], [centroids_failed, new_centroids2], minmax=(RANGE_MIN, RANGE_MAX))
    #best:
    #display_simple_using_mayavi_([(verts, faces), (verts2, faces2)], [centroids, new_centroids], minmax=(RANGE_MIN, RANGE_MAX))
    display_simple_using_mayavi_([(verts, faces), (verts2, faces2)], [centroids, new_centroids, centroids_failed], minmax=(RANGE_MIN, RANGE_MAX))


if __name__ == '__main__':
    test_proj_ohtak2()  # 837 (success)  +  245  (failed)

