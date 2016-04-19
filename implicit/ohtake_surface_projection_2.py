import __builtin__
try:
    __builtin__.profile
except AttributeError:
    # No line profiler, provide a pass-through version
    def profile(func): return func
    __builtin__.profile = profile



import numpy as np

from implicit_config import config

from basic_types import check_vector4_vectorized

"""
Order (sequence) of operation:
------
project_single_point2_ohtake -> search_near_using_vector_ohtake
search_near_using_gradient_ohtake -> bisection_prop_2 -> bisection_3
? -> search_near_1d_ohtake -> bisection_prop_2 -> bisection_3

search_near_1d_ohtake and search_near_using_gradient_ohtake are combined into
search_near__ohtake

Mapping between this and Ohtake's implementation
------
project_single_point2_ohtake ->
search_near_using_vector_ohtake -> ??
search_near_using_gradient_ohtake -> searchNearPoint+lambda_hops(grad,MAX)

search_near_1d_ohtake -> searchNearPoint1D +lambda_hops(MAX)
bisection_prop_2 -> proportional_bisection (MAX/2)
bisection_3 -> bisection (2MAX)

Ohtake's call sequence:
setCenterOnSurface -> searchNearPoint-> searchNearPoint1D.
searchNearPoint = ( 0, lambda_hops(grad,MAX) -> proportional_bisection (MAX/2) -> bisection (2MAX) -> giveup ).
searchNearPoint1D = ( 0, lambda_hops(MAX) -> proportional_bisection (MAX/2) -> bisection (2MAX) )

warning: max_dist is not taken into account in Ohtake
"""


def bisection_3_standard(iobj, p1, p2, f1, f2, MAX_ITER):
    TH1 = 0.001
    TH3 = 0.001

    assert p1.shape[0] == 1
    assert p2.shape[0] == 1

    assert f1 < 0
    assert f1*f2 < 0, "Opposite signs required"
    for j in range(MAX_ITER):
        if np.linalg.norm(p1-p2) < TH1:
            return p3, j
        p3 = 0.5 * (p1 + p2)
        f3 = iobj.implicitFunction(p3)
        if np.abs(f3) < TH3:
            return p3, j
        elif f1*f3 >= 0:
            p1 = p3
            f1 = f3
        else:
            p2 = p3
            f2 = f3
    #print "Convergence of the bisection did not happen"
    return None, MAX_ITER


@profile
def bisection_prop_2(iobj, p1, p2, f1, f2, MAX_ITER):
    """ The proportional bisection method. See searchNearPoint1D() """
    # The following is used with MAX_ITER==10 in Ohtake.
    TH2 = 0.0001  # distance between p1,p2
    TH3 = 0.001   # for abs(f3)

    assert p1.shape[0] == 1
    assert p2.shape[0] == 1

    dt = p2 - p1

    assert f1 < 0
    assert f1*f2 < 0, "Opposite signs required"

    for j in range(MAX_ITER):

        assert f1 < 0
        assert f1*f2 < 0
        fp = f1/(f1-f2)
        p3 = p1 + dt * fp
        f3 = iobj.implicitFunction(p3)

        if np.abs(f3) < TH3:
            return True, p3, None, j

        #print f1, f2, " -> fp:", fp  # Sometimes (0.96) it searches too close to f2, and fp converges to 1
        if f3 < 0.:
            if VERBOSE:
                print "A",
            (p1, f1) = (p3, f3)
        else:
            if VERBOSE:
                print "B",
            (p2, f2) = (p3, f3)

        dt = p2 - p1
        #should be moved above

        if np.linalg.norm(dt) < TH2:
            return True, p3, None, j

    if VERBOSE:
        print "Convergence based on the proportional method did not happen"
    #return (False, p1, p2)
    p,ni = bisection_3_standard(iobj, p1, p2, f1, f2, MAX_ITER*4)
    if p is not None:
        if VERBOSE:
            print "Bisection converged "
        return True, p, None, 0  # very unlikely
    else:
        return False, p1, p2, -ni


VERBOSE = False

@profile
def func_test_bisection_p(iobj, sp1, sp2, prop, max_iter_count):

    global count_converged
    global count_not_converged


    from basic_types import check_vector4
    check_vector4(sp1)
    check_vector4(sp2)
    p1 = sp1.reshape((1, 4))
    p2 = sp2.reshape((1, 4))
    f1 = iobj.implicitFunction(p1)
    f2 = iobj.implicitFunction(p2)
    print "f1,f2: ", f1, f2

    #p = bisection_prop_2(iobj, p1, p2, f1, f2, 20)
    if prop:
        converged, p, p2, jj = bisection_prop_2(iobj, p1, p2, f1, f2, max_iter_count)
        assert not p is None
    else:
        p, jj = bisection_3_standard(iobj, p1, p2, f1, f2, max_iter_count)
        converged = p is not None

    if not converged:
        count_not_converged += 1
        if prop:
            if VERBOSE:
                print "*** bisection_prop_2 did not converge. Hence the result is an interval"  # this happens
                print (p, p2, jj)
        if not prop:
            if VERBOSE:
                print "*** bisection_standard did not converge. Hence the result is None"
                print p, jj
    else:
        print(count_converged)
        count_converged += 1
        if VERBOSE:
            print "converged: p", p
            print "iteration: ", jj
        from basic_types import make_vector4_vectorized
        p4 = make_vector4_vectorized(p[0, 0], p[0, 1], p[0, 2])
        f = iobj.implicitFunction(p4)
        if VERBOSE:
            print("f(p)=", f)
            print("Total distance travelled: ", np.linalg.norm(sp1 - p), np.linalg.norm(sp2 - p) )


def test_bisection(prop, num_pairs, max_iter_count):
    import example_objects
    iobj = example_objects.make_example_vectorized("first_csg")
    from basic_types import make_vector4
    func_test_bisection_p(iobj, make_vector4(1, 1, 1), make_vector4(0, 0, 0), prop,  max_iter_count =max_iter_count)
    if VERBOSE:
        print("=============")
    func_test_bisection_p(iobj, make_vector4(10, 10, 10), make_vector4(0, 0, 0), prop,     max_iter_count =max_iter_count)
    if VERBOSE:
        print("=============")
    counter = 0
    while counter < num_pairs:
        from basic_types import make_random_vector_vectorized
        x1 = make_random_vector_vectorized(1, 90/10, 3, "randn", normalize=False)
        x2 = make_random_vector_vectorized(1, 90/10, 3, "randn", normalize=False)
        f1 = iobj.implicitFunction(x1)
        f2 = iobj.implicitFunction(x2)
        if f1*f2 < 0 and f1 < 0:
            func_test_bisection_p(iobj, x1[0, :], x2[0, :], prop, max_iter_count=max_iter_count)
            if VERBOSE:
                print("======================")
                import sys
                sys.stdout.flush()
            counter += 1




def plot_maxiter_vs_convergence():
    """ Plots the probability (count of successful convergence) given different values of MAXITER.
    This should also depend on the mean distance between the points. """
    y1 = np.zeros((100,))
    y2 = np.zeros((100,))
    x = np.zeros((100,)) + np.nan
    for i in range(13):
        maxiter = i

        global count_converged
        global count_not_converged
        count_converged = 0
        count_not_converged = 0

        #test_bisection(True)
        test_bisection(True, 20, max_iter_count=maxiter)
        print count_not_converged, count_converged

        y1[i] = count_not_converged
        y2[i] = count_converged
        x[i] = maxiter
        print("===========================", i)

    print(y1)
    print(y2)
    #import numpy as np
    import matplotlib.pyplot as plt

    plt.plot(x, y1, 'b.-', x,y2, 'rs')
    plt.show()




import math
# Lots of repeated code between search_near_using_gradient_ohtake and search_near_using_gradient_ohtake

@profile
#search_near_1d_ohtake
def search_near_ohtake(iobj, start_x, direction, lambda_val, MAX_ITER):  # max_dist
    """Returns either the point, or None, if not found. lambda_val is the expected distance from the surface.
    The resommended value is half of the average edge length, but a better way is half of the MC'step size (because the expected error is half of the MC grid voxel size).
    Remeber: we may be on a totally irrelevant direction here.
    'direction' should be normalised. lambda*direction is used. Note that lambda is negated by default.
    Note: along_1d mode is in fact the same. It's just initialises direction=gradient(start_x).
    Does both searchNearPoint1D and searchNearPoint()
    :param direction: description
    @param np.array direction
    """
    #lambda should be ~ expected distance?  (that's why it should be half of the average edge size, i.e. half of the MC step size)

    TH1 = 0.001
    # MAX_ITER = 20
    #TH2_L = 0.00000001  # Used by Ohtake. But too small
    TH2_L = 0.00001  #only used in the along_1d mode

    if direction is not None:
        along_1d = True
    else:
        along_1d = False

    if not along_1d:
        direction = iobj.implicitGradient(start_x)  ## what?! why start_x ??
        dn = np.linalg.norm(direction[0:3])
        if dn>0.0:  # 00000001:
            direction = direction/dn

    #print direction
    eval_count = 0

    p1 = start_x
    f1 = iobj.implicitFunction(p1)

    eval_count += 1
    f0 = f1

    if math.fabs(f1) < TH1:
        return p1

    p2 = p1  # no need actually

    negative_f1 = -1 if f1 < 0. else +1
    lambda_ = lambda_val * (-negative_f1)  # negation of input is bad practice

    exit_A = False
    while True:
        # (C) jumps back here.
        for j in range(MAX_ITER):
            p2 = p2 + lambda_ * direction
            p2[:, 3] = 1
            f2 = iobj.implicitFunction(p2)
            eval_count += 1

            if f1*f2 < 0.0:
                # (A)
                exit_A = True
                break
            p1 = p2

        else:
            #for loop ended becasue of MAX_ITER

            #(C): next iteration with adaptively decreased lambda_. Revert and start over again using a smaller lambda
            lambda_ = lambda_ / 2.
            #either quit:
            if np.abs(lambda_) < TH2_L:
                print "(B)", eval_count
                return None   #
            #or back to start
            p1 = start_x
            assert f0 * f1 >= 0
            f1 = f0  #This was missing in Ohtake, because the sign of f1 is not expected to change. So I added the assert above.
            p2 = p1

            #restart the loop
            #(C).

        #(A)
        if exit_A:  # f1*f2 < 0.0:
            break
    #(A)
    assert f1*f2 < 0.0
    if f1 > 0:
        (p1, p2) = (p2, p1)
        (f1, f2) = (f2, f1)

    assert f1 < 0
    assert f2 > 0
    assert math.fabs(f1) < 1000, str(p1)
    assert math.fabs(f2) < 1000, str(p2)

    converged, p1, p2, iter_log = bisection_prop_2(iobj, p1, p2, f1, f2, MAX_ITER/2)
    assert f1 < 0
    assert f2 > 0

    if converged:
        assert p2 is None
        return p1
    else:
        return None


if __name__ == "__main__":
    print search_near__ohtake.__doc__


def test_ohtake1():
    global count_not_converged
    global count_converged

    count_converged = 0
    count_not_converged = 0

    num_pairs = 20
    import example_objects
    iobj = example_objects.make_example_vectorized("first_csg")
    from basic_types import make_vector4
    counter = 0
    while counter < num_pairs:
        from basic_types import make_random_vector_vectorized
        x1 = make_random_vector_vectorized(1, 9/4, 3, "randn", normalize=False)
        x2 = make_random_vector_vectorized(1, 9/4, 3, "randn", normalize=False)

        f1 = iobj.implicitFunction(x1)
        f2 = iobj.implicitFunction(x2)
        if f1*f2 < 0 and f2 < 0:
            (f1, x1, f2, x2) = (f2, x2, f1, x1)  # not tested

        if f1*f2 < 0 and f1 < 0:
            #print "*", f1, f2
            start_x = x1
            dr = x2 - x1
            dr = dr / np.linalg.norm(dr)

            lambda_val = 0.1 # 1.0
            #max_dist = 9
            #p = project_single_point2_ohtake(iobj, start_x, lambda_val, max_dist )
            MAX_ITER = 20
            #p, p2 = search_near__ohtake(iobj, start_x, dr, lambda_val, MAX_ITER)
            #success = p2 is None
            #if not success: # p is None:
            #    #print "not"
            #    count_not_converged += 1
            #else:
            #    assert p is not None
            #    #print "yes"
            #    count_converged += 1
            #    print("Total distance travelled: ", np.linalg.norm(start_x - p), np.linalg.norm(sp2 - p) )
            #print("***********************")
            is_1d = False
            p = search_near__ohtake(iobj, start_x, dr if is_1d else None, lambda_val, MAX_ITER)
            #print p
            #print "=============="
            if p is None:
               count_not_converged += 1
               print "Did not converge ",
               print "*", f1, f2,
               print("***********************************=*=*=*=")
               #341 iterations when it doesn't find it. Such a waste of O(.)
               print 100.0*float(count_converged)/float(count_not_converged+count_converged), "%", "converged"
               #Result: only %19.74% will converge
            else:
               count_converged += 1
               print("cccc +")
               #print start_x, p
               #print "Total distance travelled: ", np.linalg.norm(start_x - p[0,:]), " f(x)=", iobj.implicitFunction(p)



@profile
def project_point_bidir_ohtake(iobj, start_x, lambda_val, max_dist):
    """ max_dist is used.
    See # setCenterOnSurface """
    #""" lambda_val: step size"""
    #max_iter = 20  # config["max_iter"]
    check_vector4_vectorized(start_x)
    assert start_x.shape[0] == 1

    max_iter=20
    #p =
    p1 = search_near__ohtake(iobj, start_x, None, lambda_val, max_iter)
    if p1 is None:
        return None  # Should we return nothing if nothing found??
    f1 = iobj.implicitFunction(p1)

    # Mirror image: search the opposite way and accept only if it is closer than the best found.
    p2 = 2*start_x - p1
    f2 = iobj.implicitFunction(p2)
    p = p1  #None #p1 # None #p1  #default
    if f1*f2 < 0:
        direction = (start_x - p1)  # as in Ohtake
        dn = np.linalg.norm(direction)
        if dn > 0:  #dn>0.000000001:
            direction = direction/dn

            #broken
            p3 = search_near__ohtake(iobj, start_x, direction, lambda_val, max_iter)
            #no max_dist

            if p3 is not None:
                if np.linalg.norm(start_x - p3) > np.linalg.norm(start_x - p1):
                    p = p3
                #else:
                #    p = p1
    #else:
    #    p = p1
    if p is not None:
     if np.linalg.norm(start_x - p) > max_dist:
        return None
    return p

# @profile
# def set_centers_on_surface__ohtake(iobj, centroids, average_edge):
#     #nones_map = centroids[:,0]*0 < 100
#     print "Projecting the centroids:"
#     for i in range(centroids.shape[0]):
#         print i,
#         e = average_edge
#         lm = e/2
#         max_dist = e
#         c = project_point_bidir_ohtake(iobj, centroids[i,np.newaxis,:], lm, max_dist)
#         if c is not None:
#             centroids[i] = c
#         else:
#             nones_map[i] = True
#     #print nones_map



def set_centers_on_surface_ohtake(iobj, centroids, average_edge):
    #here we consider that the max_dist is the average_edge and lambda = average_edge/2
    #new function who is a combination of set_centers_on_surface_ohtake and project_point_bidir_ohtake
    lambda_val = average_edge/2
    check_vector4_vectorized(centroids)
    #definition of the matrix that are gonna be used in the rest of the programm
    p1 = np.ndarray(centroids.shape)
    p2 = np.ndarray(centroids.shape)
    p3 = np.ndarray(centroids.shape)
    p = np.ndarray(centroids.shape)
    f1 = np.ndarray(centroids.shape[0])
    f2 = np.ndarray(centroids.shape[0])
    direction = np.ndarray(centroids.shape)
    dn = np.ndarray(centroids.shape[0])

    max_iter = 20 # used in search_near_ohtake

    for i in range(centroids.shape[0]):
        p1[i,:] = search_near_ohtake(iobj, centroids[i,:].reshape((1,4)), None, lambda_val, max_iter)
        if np.allclose(p1[i, 3], 1, 0.00000000000001) == True: #make sure that p are found by the program and they respect the condition enforce by check_vector4_vectorized
            p1[i,:].reshape(1,4)
            f1[i] = iobj.implicitFunction(p1[i,:].reshape(1,4))
            #print f1[i]
                # Mirror image: search the opposite way and accept only if it is closer than the best found.
            p2[i,:] = 2*centroids[i,:] - p1[i,:] #p2 correspond to S in the paper
            f2[i] = iobj.implicitFunction(p2[i,:].reshape(1,4))
            p[i,:] = p1[i,:]

            if f1[i]*f2[i] < 0:
                direction[i,:] = (centroids[i,:] - p1[i,:])  # as in Ohtake
                dn[i] = np.linalg.norm(direction[i,:])
                if dn[i] > 0:  #dn>0.000000001:
                    direction[i,:] = direction[i,:]/dn[i]
                    p3[i,:] = search_near_ohtake(iobj, centroids[i,:].reshape((1,4)), direction[i,:], lambda_val, max_iter)
                    #no max_dist

                    if p3[i,:] is not None and np.allclose(p3[i, 3], 1, 0.00000000000001) == True:
                        if np.linalg.norm(centroids[i,:] - p3[i,:]) > np.linalg.norm(centroids[i,:] - p1[i,:]):
                            p[i,:] = p3[i,:]
                            #else:
                            #    p = p1
                #else:
                # #    p = p1
            if p[i,:] is not None:
                if np.linalg.norm(centroids[i,:] - p[i,:]) <= average_edge:
                    centroids[i,:] = p[i,:]

def display_simple_using_mayavi_(vf_list, pointcloud_list, minmax=(-1,1), mayavi_wireframe=False, opacity=1.0):
    """ opacity can be either a list of a constant """
    from mayavi import mlab

    if type(opacity) is list:
        opacities = opacity  # 1.0
    else:
        opacities = [opacity] + [0.2]*(len(vf_list)-1)  # 1.0, 0.2 #0.1

    i = 0
    for vf in vf_list:
        verts, faces = vf
        if verts is None:
            continue
        if verts.size == 0:
            print("Warning: empty vertices")
            continue
        if faces.size == 0:
            print("Warning: no faces")
            continue
        mlab.triangular_mesh([vert[0] for vert in verts],
                         [vert[1] for vert in verts],
                         [vert[2] for vert in verts],faces,representation="surface" if not mayavi_wireframe else "wireframe",
                         opacity=opacities[i], scale_factor = 100.0)
        i += 1


    color_list = [(1, 0, 0), (0, 0, 0), (1, 1, 0), (0, 0, 1), (0,1,0)]
    i = 0
    for c in pointcloud_list:
        #print c[:,0:3]
        mlab.points3d(c[:, 0], c[:, 1], c[:, 2], color=color_list[i], scale_factor=0.2)
        i+=1
    del i


    #if pointcloud_list[0].shape[0] == pointcloud_list[1].shape[0]:
    #    print "PLOT3d"
    #    c0 = pointcloud_list[0][np.newaxis, :, :]   # 2 x N x 4
    #    c1 = pointcloud_list[1][np.newaxis, :, :]
    #    c01 = np.concatenate((c0, c1), axis=0)
    #    for i in range(c01.shape[1]):
    #        #print c01[:,i,0], c01[:,i,1], c01[:,i,3]
    #        mlab.plot3d(c01[:,i,0],c01[:,i,1],c01[:,i,3])
    #        #exit()
    #        pass
    #    #mlab.plot3d((c01[:,:,0]), (c01[:,:,1]), (c01[:,:,3]))


    (RANGE_MIN, RANGE_MAX) = minmax
    x = np.linspace(RANGE_MIN,RANGE_MAX,2).reshape(2,1)
    y = np.zeros((2,1))
    z = np.zeros((2,1))

    mlab.plot3d(x, y, z,line_width=3,name="x-axis")
    mlab.plot3d(y, x, z,line_width=3,name="y-axis")
    mlab.plot3d(z, y, x,line_width=3,name="z-axis")

    mlab.text3d(RANGE_MAX,0,0, "x", scale=0.3)
    mlab.text3d(0,RANGE_MAX,0, "y", scale=0.3)
    mlab.text3d(0,0,RANGE_MAX, "z", scale=0.3)
    mlab.text3d(RANGE_MIN,0,0, "-x", scale=0.3)
    mlab.text3d(0,RANGE_MIN,0, "-y", scale=0.3)
    mlab.text3d(0,0,RANGE_MIN, "-z", scale=0.3)

    mlab.show()


def test_proj_ohtak1():
    """ Tests projection of centroids using Ohtake's original method. Projects point-wise.
    Repeats points by applying a small perturbation.
    Blend disc object.
    Uses a high resolution MC for checking the new projected points (black dots).
    The original centroids are shon in red.
    The failed projections are shown in Yello. Their number is writtn in output (9 points)"""

    from stl_tests import make_mc_mesh_scikit

    exname = "blend_example2_discs"  # "blend_example2"
    import example_objects
    iobj = example_objects.make_example_vectorized(exname, 8.0)

    #(RANGE_MIN, RANGE_MAX, STEPSIZE) = (-2.*8, +4.*8, 0.4*8/5)
    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-16, +32, 0.64*3*2/2)
    verts, faces = make_mc_mesh_scikit(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE)

    print verts.shape, faces.shape

    average_edge = STEPSIZE

    # verts = optimise_mesh(verts, faces, iobj)

    c3 = np.mean(verts[faces[:], :], axis=1)
    # add extra points
    c3 = np.concatenate((c3, c3+STEPSIZE*0.1, c3+STEPSIZE*(-0.1)), axis=0)
    c3 = np.concatenate((c3,), axis=0)
    centroids = np.concatenate((c3, np.ones((c3.shape[0], 1))), axis=1)


    new_centroids = centroids.copy()
    set_centers_on_surface__ohtake(iobj, new_centroids, average_edge)

    print
    #print np.arange(nones_map.shape[0])[np.logical_not(nones_map)].shape
    print np.sum(np.logical_not(nones_map)), "(success)  + ", np.sum(nones_map), " (failed)"
    new_centroids2 = new_centroids[np.logical_not(nones_map), :]

    #todo: Use "None" for non-converging ones. (and visualise them)

    #(RANGE_MIN, RANGE_MAX, STEPSIZE) = (-16, +32, 0.64*3*2/2)
    #verts2, faces2 = make_mc_mesh_scikit(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE * 0.2/2.0) #15 million!
    verts2, faces2 = make_mc_mesh_scikit(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE * 0.2)


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



def test_profiler():
    #use for profiling only
    exname = "blend_example2_discs"  # "blend_example2"
    import example_objects
    iobj = example_objects.make_example_vectorized(exname, 8.0)
    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-16, +32, 0.64*3*2/2)
    verts, faces = make_mc_mesh_scikit(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE)
    print verts.shape, faces.shape
    average_edge = STEPSIZE
    c3 = np.mean(verts[faces[:], :], axis=1)
    c3 = np.concatenate( (c3,c3+STEPSIZE*0.1,c3+STEPSIZE*(-0.1)), axis=0)
    centroids = np.concatenate((c3, np.ones((c3.shape[0], 1))), axis=1)
    nones_map = centroids[:,0]*0 > 100
    new_centroids = centroids.copy()
    set_centers_on_surface__ohtake(iobj, new_centroids, average_edge, nones_map)
    new_centroids2 = new_centroids[np.logical_not(nones_map),:]
    verts2, faces2 = make_mc_mesh_scikit(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE * 0.2/1.0)





#  ********************************** Vectorized version ****************************

from basic_types import make_vector4_vectorized

def bisection_3_standard_vectorized(iobj, p1, p2, f1, f2, MAX_ITER):
    TH1 = 0.001
    TH3 = 0.001
    assert p1.shape[0] == 1
    assert p2.shape[0] == 1

    raise ImplementationError()

    assert np.all(f1 < 0)
    assert np.all(f1*f2 < 0), "Opposite signs required"
    for j in range(MAX_ITER):
        if np.linalg.norm(p1-p2, axis=1) < TH1:
            return p3, j
        p3 = 0.5 * (p1 + p2)
        f3 = iobj.implicitFunction(p3)
        if np.abs(f3) < TH3:
            return p3, j
        elif f1*f3 >= 0:
            p1 = p3
            f1 = f3
        else:
            p2 = p3
            f2 = f3
    #print "Convergence of the bisection did not happen"
    return None, MAX_ITER


@profile
def func_test_bisection_p_vectorized(iobj, sp1, sp2, prop, max_iter_count):

    global count_converged
    global count_not_converged

    raise ImplementationError()


    from basic_types import check_vector4
    check_vector4_vectorized(sp1)
    check_vector4_vectorized(sp2)
    n = sp1.shape[0]
    p1 = sp1.reshape((n, 4))
    p2 = sp2.reshape((n, 4))
    f1 = iobj.implicitFunction(p1)
    f2 = iobj.implicitFunction(p2)
    print "f1,f2: ", f1, f2

    #p = bisection_prop_2(iobj, p1, p2, f1, f2, 20)
    if prop:
        converged, p, p2, jj = bisection_prop_2(iobj, p1, p2, f1, f2, max_iter_count)
        assert not p is None
    else:
        p, jj = bisection_3_standard_vectorized(iobj, p1, p2, f1, f2, max_iter_count)
        converged = p is not None

    if not converged:
        count_not_converged += 1
        if prop:
            if VERBOSE:
                print "*** bisection_prop_2 did not converge. Hence the result is an interval"  # this happens
                print (p, p2, jj)
        if not prop:
            if VERBOSE:
                print "*** bisection_standard did not converge. Hence the result is None"
                print p, jj
    else:
        print(count_converged)
        count_converged += 1
        if VERBOSE:
            print "converged: p", p
            print "iteration: ", jj
        from basic_types import make_vector4_vectorized
        p4 = make_vector4_vectorized(p[0, 0], p[0, 1], p[0, 2])
        f = iobj.implicitFunction(p4)
        if VERBOSE:
            print("f(p)=", f)
            print("Total distance travelled: ", np.linalg.norm(sp1 - p), np.linalg.norm(sp2 - p) )


def test_bisection_vectorized(prop, num_pairs, max_iter_count):
    import example_objects
    iobj = example_objects.make_example_vectorized("first_csg")
    from basic_types import make_vector4
    func_test_bisection_p_vectorized(iobj, make_vector4_vectorized(1, 1, 1), make_vector4_vectorized(0, 0, 0), prop,  max_iter_count =max_iter_count)
    if VERBOSE:
        print("=============")
    func_test_bisection_p_vectorized(iobj, make_vector4_vectorized(10, 10, 10), make_vector4_vectorized(0, 0, 0), prop,     max_iter_count =max_iter_count)
    if VERBOSE:
        print("=============")
    counter = 0
    while counter < num_pairs:
        from basic_types import make_random_vector_vectorized
        x1 = make_random_vector_vectorized(1, 90/10, 3, "randn", normalize=False)
        x2 = make_random_vector_vectorized(1, 90/10, 3, "randn", normalize=False)
        f1 = iobj.implicitFunction(x1)
        f2 = iobj.implicitFunction(x2)
        if f1*f2 < 0 and f1 < 0:
            func_test_bisection_p_vectorized(iobj, x1[:, :], x2[:, :], prop, max_iter_count=max_iter_count)
            if VERBOSE:
                print("======================")
                import sys
                sys.stdout.flush()
            counter += 1


# *********************************** main **********************************


if __name__ == '__main__':
    if False:
        #not needed:
        #global count_converged
        #global count_not_converged
        count_converged = 0
        count_not_converged = 0

        test_bisection(True,10,10)
        test_bisection(False,10,10)
        print(count_converged, count_not_converged)
        print("fine")



    #plot_maxiter_vs_convergence()

    #test_ohtake1()
    #Result: only %19.74% will converge

    if True:
        test_proj_ohtak1() #tha main test
    #test_profiler()

    if False:  # vectorized
        count_converged = 0
        count_not_converged = 0
        test_bisection_vectorized(True, 10, 10)
        test_bisection_vectorized(False, 10, 10)
