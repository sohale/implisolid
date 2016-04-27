from ipdb import set_trace
import profile_support
from vtk_mc import vtk_mc

import sys
import math

import numpy as np
from basic_types import check_vector4_vectorized, normalize_vector4_vectorized

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


#@profile
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
def search_near_ohtake_old(iobj, start_x, direction, lambda_val, MAX_ITER):  # max_dist
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

            if not along_1d:
                direction = iobj.implicitGradient(start_x)  ## what?! why start_x ??
                dn = np.linalg.norm(direction)
                if dn>0.0:  # 00000001:
                    direction = direction/dn
                else:
                    pass  # Finding is not going to happen. But it's fine.

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
    #            print i, "(B)", eval_count
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


@profile
def search_near_ohtake(iobj, start_x, lambda_val, MAX_ITER):  # max_dist
    """Vectorized function of the search_near_ohtake_old who take a vector (start_x)
     and return the vector p1_found if we found the centroids"""

    #lambda should be ~ expected distance?  (that's why it should be half of the average edge size, i.e. half of the MC step size)
    TH1 = 0.001
    TH2_L = 0.00001

    #start_x_new = start_x.copy()

    #p1_found = np.ndarray(start_x.shape)
    p1_found = np.zeros(start_x.shape)
#    p2 = np.ndarray(start_x.shape)
#    f0 = np.ndarray(start_x.shape[0])
#    f2 = np.ndarray(start_x.shape[0])
    direction = np.ndarray(start_x.shape)

    p1 = start_x.copy() #adding the copy : solve one problem of the code
    #p1 = start_x
    f1 = iobj.implicitFunction(p1)

    lambda_ = lambda_val * (-np.sign(f1))

#    counting = 0

    for i in range(start_x.shape[0]):
        xxv = start_x[i,:].reshape(1,4)
        check_vector4_vectorized(xxv)
        direction[i,:] = (iobj.implicitGradient(xxv)) #no need for copy here!
        dn = np.linalg.norm(direction[i,0:3]) #calculation of the normalization

        if  dn > 0.0:
            direction[i,:] = (direction[i,:]/dn)
            direction[i,3] = 1

#        f0 = f1[i]
        eval_count = 0
        eval_count += 1


        if math.fabs(f1[i]) < TH1:
            #print i, "Out"
            p1_found[i,:] = p1[i,:]

        else:

        #    negative_f1 = -1 if f1[i] < 0. else 1
        #    lambda_ = lambda_val * (-negative_f1)  # negation of input is bad practice

            p2 = p1[i,:].copy()
            exit_A = False
            not_found = False
            number_of_iteration = 0
            s1 = start_x[i,:].copy()

            while True:

                # (C) jumps back here.
                for j in range(MAX_ITER):
                    number_of_iteration += 1

                    xxv2 = start_x[i,:].copy().reshape(1,4)
                    check_vector4_vectorized(xxv2)
                    dir1 = iobj.implicitGradient(xxv2) #we have to give a copy of start_x to avoid problems
                    dir1 = dir1[0,:]
                    dn1 = np.linalg.norm(dir1[0:3].copy())

                    if  dn1 > 0.0:
                        dir1 = dir1/dn1
                    dir1[3] = 1
                    dir2 = direction[i,:]


                    # if not np.allclose(s2, s1):
                    #     print i, "222"
                    #     print start_x_new[i].shape, start_x_new[i].dtype
                    #     print s1.shape, s1.dtype
                    if not np.allclose(dir1,dir2):
                        print i,"patate!"
                        print dir1, dir2
                        print dir1/dir1[0], dir2/dir2[0]
                        print "ALLCLOSE failed"
                        exit()


                    p2 = p1[i,:] + lambda_[i] * dir2
            #        p2[i,:] = p1[i,:] + lambda_[i] * direction[i]

                    p2[3] = 1
                    f2 =iobj.implicitFunction(p2.reshape(1,4))
                    eval_count += 1

                    if f1[i]*f2 < 0.0:
                        exit_A = True #(A)
                        break

                    p1[i,:] = p2

                    if number_of_iteration == MAX_ITER:
                        number_of_iteration = 0
                        lambda_[i] = lambda_[i] / 2.

                    #    print lambda_[i], TH2_L

                        if np.abs(lambda_[i]) < TH2_L:
                            print i, "(B)", eval_count
                            p1_found[i,3] = 0
                            not_found = True
                            break
                        else:
                            p1[i,:] = start_x[i,:].copy()
                            # assert f0 * f1[i] >= 0
                            # f1[i] = f0  #This was missing in Ohtake, because the sign of f1 is not expected to change. So I added the assert above.
                            p2 = p1[i,:].copy()

                if exit_A:
                    break
                if not_found: #np.abs(lambda_)< TH2_l
                    break
    #
            if not not_found:

                #import ipdb; ipdb.set_trace()
                assert f1[i]*f2 < 0.0
                if f1[i] > 0:
                    (p1[i,:], p2) = (p2, p1[i,:])
                    (f1[i], f2) = (f2, f1[i])

                assert f1[i] < 0
                assert f2 > 0
                assert math.fabs(f1[i]) < 1000, str(p1[i,:])
                assert math.fabs(f2) < 1000, str(p2)

                converged, p1[i,:], p2, iter_log = bisection_prop_2(iobj, p1[i,:].reshape(1,4), p2.reshape(1,4), f1[i], f2, MAX_ITER/2)
                assert f1[i] < 0
                assert f2 > 0

                if converged:
                    p1_found[i,:] = p1[i,:].copy()

                else:
                    p1_found[:,3] = 0
                    p1_found[i,3] = 0

                #    print p1_found[i,:]
        #    print direction[i,:], np.linalg.norm(direction[i,:])
        # if f1[i] == f0:
        #     counting += 1
    #    print p1_found[i] - start_x[i]
#    print counting
    return p1_found

@profile
def set_centers_on_surface_ohtake(iobj, centroids, average_edge):
#    here we consider that the max_dist is the average_edge and lambda = average_edge/2
#    new function who is a combination of sers_on_surface_ohtake and project_point_bidir_ohtake
    lambda_val = average_edge/2
    check_vector4_vectorized(centroids)
    #definition of the matrix that are gonna be used in the rest of the programm
    p1 = np.ndarray(centroids.shape)
    p2 = np.ndarray(centroids.shape)
    p3 = np.ndarray(centroids.shape)
    p = np.ndarray(centroids.shape)
    f1 = np.ndarray(centroids.shape[0])
    f2 = np.ndarray(centroids.shape[0])
    direction_3 = np.ndarray(centroids.shape)
    dn_3 = np.ndarray(centroids.shape[0])


    centroids_new = centroids.copy()
    max_iter = 20

    for i in range(centroids.shape[0]):
        #import ipdb; ipdb.set_trace()
        p1[i,:] = search_near_ohtake_old(iobj, centroids[i,:].reshape((1,4)), None, lambda_val, max_iter)
        if np.allclose(p1[i, 3], 1, 0.00000000000001) == True: #make sure that p are found by the program and they respect the condition enforce by check_vector4_vectorized
            p1[i,:].reshape(1,4)
            f1[i] = iobj.implicitFunction(p1[i,:].reshape(1,4))

                # Mirror image: search the opposite way and accept only if it is closer than the best found.
            p2[i,:] = 2*centroids_new[i,:] - p1[i,:] #p2 correspond to S in the paper
            f2[i] = iobj.implicitFunction(p2[i,:].reshape(1,4))
            p[i,:] = p1[i,:]

#            print f1[i], f2[i]

            if f1[i]*f2[i] < 0:
                direction_3[i,:] = (centroids_new[i,:] - p1[i,:].copy())  # as in Ohtake
                dn_3[i] = np.linalg.norm(direction_3[i,:])
                if dn_3[i] > 0:  #dn>0.000000001:
                    direction_3[i,:] = direction_3[i,:]/dn_3[i]
                    p3 = search_near_ohtake_old(iobj, centroids[i,:].reshape(1,4), direction_3[i,:].reshape(1,4), lambda_val, max_iter)


                    #no max_dist
                if p3 is not None:
                    if np.allclose(p3[:,3], 1, 0.00000000000001) == True:
                        if np.linalg.norm(centroids[i,:] - p3) > np.linalg.norm(centroids[i,:] - p1[i,:]):
                            p[i,:] = p3
                            #else:
                            #    p = p1
                #else:
                # #    p = p1
            #if p[i,:] is not None:
            if np.linalg.norm(centroids[i,:] - p[i,:]) <= average_edge:
                centroids[i,:] = p[i,:]






#evaluate_centroid_gradients
def compute_centroid_gradients(centroids, iobj, normalise=True):
    assert centroids is not None
    check_vector4_vectorized(centroids)
    centroid_gradients = iobj.implicitGradient(centroids)
    assert not np.any(np.isnan(centroid_gradients))
    assert not np.any(np.isinf(centroid_gradients))
    if normalise:
        centroid_normals = normalize_vector4_vectorized(centroid_gradients)
        return centroid_normals
    else:
        return centroid_gradients

def visualise_gradients(mlab, pos, iobj, arrow_size):
    lm = arrow_size  # 1.  # STEPSIZE
    pos4 = np.concatenate((pos, np.ones((pos.shape[0],1))),axis=1)
    pnormals = - iobj.implicitGradient(pos4)
    pnormals = normalize_vector4_vectorized(pnormals)
    check_vector4_vectorized(pos4)
    xyz = pos
    uvw = pnormals [:,0:3] / 2.
    xx, yy, zz = xyz[:, 0], xyz[:, 1], xyz[:, 2]
    uu, vv, ww = uvw[:, 0], uvw[:, 1], uvw[:, 2]
    mlab.quiver3d(xx, yy, zz, uu, vv, ww, color=(0, 0, 0), scale_factor=np.abs(lm), line_width=0.5)


def display_simple_using_mayavi_2(vf_list, pointcloud_list, minmax=(-1,1), mayavi_wireframe=False, opacity=1.0,
        separate=True, gradients_at=None, gradients_from_iobj=None, pointsizes=None, pointcloud_opacity=1.):
    """Two separate panels"""

    print("Mayavi."); sys.stdout.flush()

    from mayavi import mlab

    if pointsizes is None:
        pointsizes = [0.2]*10

    if type(opacity) is list:
        opacities = opacity  # 1.0
    else:
        opacities = [opacity] + [0.2]*(len(vf_list)-1)  # 1.0, 0.2 #0.1


    for fi in range(len(vf_list)):
        if separate:
            mlab.figure()

        vf = vf_list[fi]
        verts, faces = vf

        if verts is None:
            continue
        if verts.size == 0:
            print("Warning: empty vertices")
            continue
        if faces.size == 0:
            print("Warning: no faces")
            continue

        assert verts.ndim == 2
        assert faces.ndim == 2
        assert verts.shape == (verts.shape[0], 3), str(verts.shape)
        assert faces.shape == (faces.shape[0], 3), str(faces.shape)
        if type(mayavi_wireframe) is list:
            wire_frame1 = mayavi_wireframe[fi]
            assert len(mayavi_wireframe) == len(vf_list)
        else:
            wire_frame1 = mayavi_wireframe
        mlab.triangular_mesh([vert[0] for vert in verts],
                         [vert[1] for vert in verts],
                         [vert[2] for vert in verts],
                         faces,
                         representation="surface" if not wire_frame1 else "wireframe",
                         opacity=opacities[fi], scale_factor = 100.0)
        #opacity = 0.2 #0.1


        #allpoints are plottedon all panels?
        color_list = [(1, 0, 0), (0, 0, 0), (1, 1, 0), (0, 0, 1), (0,1,0)]
        i = 0
        for c in pointcloud_list:
            #if separate:
            #    if i != fi:
            #        continue
            #print c[:,0:3]
            mlab.points3d(c[:, 0], c[:, 1], c[:, 2], color=color_list[i], scale_factor=pointsizes[i], opacity=pointcloud_opacity )
            i+=1
        del i

        if minmax is not None:
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


    def add_random_interior_points(ax, iobj):
        """ Adding random points """
        n=10000
        import basic_types
        print avg_edge_len, "WHY USED BEFORE DEFINED?"
        ampl = avg_edge_len
        #ampl = 2
        x = basic_types.make_random_vector_vectorized(n, ampl, 1, type="rand", normalize=False)
        v = iobj.implicitFunction(x)
        x_sel =  x[ v >= 0 , :]
        if x_sel.size ==0:
            print("No points")
            return
        ax.points3d(x_sel[:,0], x_sel[:,1], x_sel[:,2], color=(0,0,0), scale_factor=0.2)

    if gradients_at is not None:
        verts1, faces1 = vf_list[0]
        avg_edge_len = compute_average_edge_length(verts, faces)
        visualise_gradients(mlab, gradients_at, gradients_from_iobj, avg_edge_len / 20.)
    if gradients_from_iobj is not None:
        add_random_interior_points(mlab, gradients_from_iobj)

    mlab.show()
    return

def compute_average_edge_length(verts, faces):
    nfaces = faces.shape[0]
    expand = verts[faces, :]
    assert expand.shape == (nfaces, 3, 3)
    assert expand[:, 2, :].shape == (nfaces, 3)
    ea_sum = 0.
    for i in range(3):
        i1 = i
        i2 = (i+1) % 3
        e1 = np.linalg.norm(expand[:, i1, :] - expand[:, i2, :])
        ea_sum += np.mean(e1)
    return ea_sum / 3.



#@profile
def get_A_b(vertex_id, nlist_numpy, centroids, centroid_gradients):

    nai = nlist_numpy
    center_array = centroids[nai, :]

    #note some centers may not be projected successfully in the previous step
    not_projected_successfully = np.isnan(center_array[:].ravel())
    if np.any(not_projected_successfully):
        pass

    normals = centroid_gradients[nai, :]  #why do we have repeats??
    #note : not normalised. But it works.

    norms = np.linalg.norm(normals, ord=2, axis=1)
    #can be either 0, 1 or Nan
    if np.any(norms < 0.000001):  #can be exactly 0.0
        print("Error: bad normal", normals)

    TH_N = 0.0000001  # 0.000001 = I millions of millimeter = 1 nanometer
    #can be 0,0,0, inf, nonsharp, degenerate, ...
    degenerate_normals = np.logical_or(np.isnan( np.sum(normals, axis=1)), norms < TH_N )


    assert not np.any(np.isnan(normals) )
    assert not np.any(np.isinf(normals) )


    x0 = np.zeros((3, 1))

    assert normals.shape[1] == 4
    #normals = normals   # (4)x4
    #grad = Ax+b
    A = np.zeros((3, 3))
    b = np.zeros((3, 1))
    #assert len(center_array) == len(normals)
    assert normals.shape == center_array.shape
    n_i = normals[:, 0:3, np.newaxis]
    p_i = center_array[:, 0:3, np.newaxis]

    A = np.dot(np.reshape(n_i,(normals.shape[0],3)).T, np.reshape(n_i,(normals.shape[0],3)))

    for i in range(normals.shape[0]):

        assert n_i[i].shape == (3, 1)
        nnt = np.dot(n_i[i], np.transpose(n_i[i]))

        assert nnt.shape == (3, 3)

        assert p_i[i].shape == (3, 1)
        b += -np.dot(nnt, p_i[i] - x0)

    return A, b



#@profile
def vertices_apply_qem3(verts, facets, centroids, vertex_neighbours_list, centroid_gradients):
    #based on quadratic_optimise_vertices(self, alpha=1.0)
    assert not centroids is None
    assert not vertex_neighbours_list is None
    assert not centroid_gradients is None

    #alpha = 1.0
    nvert = verts.shape[0]
    assert nvert == len(vertex_neighbours_list)

    result_verts_ranks = np.zeros((nvert,), dtype=int)
    assert verts.shape == (nvert, 3)
    new_verts = np.zeros((nvert, 3))

    for vertex_id in range(nvert):

        vi = vertex_id
        nlist = vertex_neighbours_list[vertex_id]
        nai = np.array(nlist)
        A, b = get_A_b(vi, nai, centroids, centroid_gradients)
        #print A, b


        u, s, v = np.linalg.svd(A)
        assert np.allclose(A, np.dot(u, np.dot(np.diag(s), v)))
        assert s[0] == np.max(s)


        tau = 10. ** 3.
        s[s / s[0] < 1.0/tau] = 0
        #print(s , s[0] , tau)
        rank = np.sum(s / s[0] > 1.0/tau)

        if not  s[0] > 0.000001:
            print("Warning! sigma_1 == 0" )
            print(s)
            print("A", A)

            #not tested
            result_verts_ranks[vi] = 0
            new_verts[vi, 0:3] = new_x[:, 0]

        assert np.all(s[:rank]/s[0] >= 1.0/tau)

        x = verts[vi, 0:3, np.newaxis]
        assert x.shape == (3, 1)

        y = np.dot(v, x).copy()
        utb = np.dot(-np.transpose(u), b)

        for i in range(rank):
            assert np.dot(-np.transpose(u), b).shape == (3,1)
            #print s[i] , 1.0/tau
            #assert s[i] >= 1.0/tau #fails when s[0] is small
            assert s[i]/s[0] >= 1.0/tau
            y[i] = utb[i] / s[i]
        new_x = np.dot(np.transpose(v), y)
        #print(x.ravel(), " -> ", new_x.ravel())
        #print("    delta=", (new_x - x).ravel())

        new_verts[vi, 0:3] = new_x[:, 0]
        #self.new_verts[vi,3] = 1

        assert x.shape == (3, 1)

        new_verts[vi, 0:3] = new_x[:, 0]

        if not np.all(np.abs(utb.ravel()[rank:] ) < 0.0001):
            #print("s", s.ravel()/s[0], "   utb", utb.ravel()/s[0])
            pass
        result_verts_ranks[vi] = rank

        #exit()
    print("max rank = ", np.max(result_verts_ranks))
    print("min rank = ", np.min(result_verts_ranks))
    if not np.min(result_verts_ranks) >= 1:
        print("Warning: assertion: np.min(result_verts_ranks) >= 1 failed." )

    if False:
        assert np.min(result_verts_ranks) >= 1
    return new_verts


import mesh_utils

@profile
def demo_combination_plus_qem():
    """ Now with QEM """
    curvature_epsilon = 1. / 1000.  # a>eps  1/a > 1/eps = 2000
    VERTEX_RELAXATION_ITERATIONS_COUNT = 0
    SUBDIVISION_ITERATIONS_COUNT = 0  # 2  # 5+4

    from example_objects import make_example_vectorized
    iobj = make_example_vectorized(
        #"rcube_vec")  #
        #"sphere_example")
        #"rdice_vec")  #
        #"cube_example");
        #"screw2")
        "cube_with_cylinders")

    #    "ell_example1")  #
        # "bowl_15_holes")  # works too. But too many faces => too slow, too much memory. 32K?
    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-3, +5, 0.2)


    import vectorized, example_objects
    c2 = vectorized.UnitCube1(1.)
    def rotate_scale_(iobj, scale, center, angle=0.):
        ns = vectorized
        import numpy
        m = numpy.eye(4)
        m[0,0] = 0.1
        iobj = ns.Transformed(iobj, m=m)
        iobj  \
            .resize(scale) \
            .move(center[0], center[1], center[2])
        if angle != 0.:
            iobj.rotate(angle, along=make_vector4(1, 1, 1), units="deg")
        return iobj

    c2 = rotate_scale_(c2, 2., [1,1,1])
    iobj = vectorized.CrispUnion( example_objects.rcube_vec(1.), c2 )


    from stl_tests import make_mc_values_grid
    gridvals = make_mc_values_grid(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE, old=False)
    verts, facets = vtk_mc(gridvals, (RANGE_MIN, RANGE_MAX, STEPSIZE))
    print("MC calculated.");sys.stdout.flush()

    old_verts, old_facets = verts, facets
    #
    # display_simple_using_mayavi_2( [(verts, facets),(verts, facets), ],
    #    pointcloud_list=[],
    #    mayavi_wireframe=[False, True,], opacity=[0.8, 0.2], gradients_at=None, separate=False, gradients_from_iobj=None,
    #    minmax=(RANGE_MIN,RANGE_MAX)  )
    # exit()

    for i in range(VERTEX_RELAXATION_ITERATIONS_COUNT):
        verts, facets_not_used, centroids = process2_vertex_resampling_relaxation(verts, facets, iobj)
        print("Vertex relaxation applied.");sys.stdout.flush()


    total_subdivided_facets = []
    for i in range(SUBDIVISION_ITERATIONS_COUNT):
        e_array, bad_facets_count = compute_facets_subdivision_curvatures(verts, facets, iobj)

        e_array[np.isnan(e_array)] = 0  # treat NaN curvatures as zero curvature => no subdivision

        which_facets = np.arange(facets.shape[0])[ e_array > curvature_epsilon ]

        verts4_subdivided, facets3_subdivided = subdivide_multiple_facets(verts, facets, which_facets)
        global trace_subdivided_facets  # third implicit output
        verts, facets = verts4_subdivided, facets3_subdivided
        print("Subdivision applied.");sys.stdout.flush()

        total_subdivided_facets += trace_subdivided_facets  # old face indices remain valid

        for i in range(VERTEX_RELAXATION_ITERATIONS_COUNT):
            verts, facets_not_used, centroids = process2_vertex_resampling_relaxation(verts, facets, iobj)
            print("Vertex relaxation applied.");sys.stdout.flush()

    average_edge = compute_average_edge_length(verts, facets)

    c3 = np.mean(verts[facets[:], :], axis=1)
    old_centroids = np.concatenate((c3, np.ones((c3.shape[0], 1))), axis=1)

    new_centroids = old_centroids.copy()
    set_centers_on_surface_ohtake(iobj, new_centroids, average_edge)
    #new_centroids is the output

    # display_simple_using_mayavi_2( [(verts, facets),(verts, facets), ],
    #    pointcloud_list=[ new_centroids ], pointcloud_opacity=0.2,
    #    mayavi_wireframe=[False, True,], opacity=[1, 1, 0.9], gradients_at=None, separate=False, gradients_from_iobj=None,
    #    minmax=(RANGE_MIN,RANGE_MAX)  )
    # exit()



    # The two CHOICEs are equaivalent. Two rewrite of the same method.
    CHOICE = 1
    if CHOICE == 1:
        #neighbour_faces_of_vertex
        vertex_neighbours_list = mesh_utils.make_neighbour_faces_of_vertex(facets)
    #    import ipdb; ipdb.set_trace()
        centroid_gradients = compute_centroid_gradients(new_centroids, iobj)
        #nv1  =
        new_verts_qem = \
            vertices_apply_qem3(verts, facets, new_centroids, vertex_neighbours_list, centroid_gradients)
        #verts = nv1
        #new_verts_qem = verts


    alpha = 0.
    new_verts_qem_alpha = (new_verts_qem * alpha + verts * (1-alpha))

    chosen_facet_indices = np.array(total_subdivided_facets, dtype=int)

    centroids2, new_centroids2 = old_centroids[chosen_facet_indices], new_centroids[chosen_facet_indices]

    # move the following code into subdivide_multiple_facets() (?)
    if chosen_facet_indices.size == 0:
        chosen_subset_of_facets = np.zeros((0,), dtype=int)
    else:
        chosen_subset_of_facets = facets[chosen_facet_indices, :]

    highlighted_vertices = np.array([131,  71, 132])  # np.arange(100, 200)
    hv = new_verts_qem[highlighted_vertices, :]

    display_simple_using_mayavi_2( [(new_verts_qem_alpha, facets),(new_verts_qem, facets), ],
       pointcloud_list=[ hv ], pointcloud_opacity=0.2,
       mayavi_wireframe=[False,False], opacity=[0.4*0, 1, 0.9], gradients_at=None, separate=False, gradients_from_iobj=None,
       minmax=(RANGE_MIN,RANGE_MAX)  )

#from timeit import default_timer as dtimer


if __name__ == '__main__':

#    timer_t1s = dtimer()

    demo_combination_plus_qem()  # subdivision + projection + qem
