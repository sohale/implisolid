#based on set_centers_on_surface__ohtake_v3s_002

from utils import optimised_used
TEST = not optimised_used()


def make_grad_directions(iobj, x):

    THRESHOLD_minimum_gradient_len = 0.000001  # kill gradients smaller than this

    # todo: Actually: better try direction toward ... : If positive, negate the gradient.

    g_a = iobj.implicitGradient(x)[:, :3]
    glen_a = np.linalg.norm(g_a, axis=1)
    glen_a[np.abs(glen_a) < THRESHOLD_minimum_gradient_len] = 1.

    g_normalization_factors = 1. / glen_a[:, np.newaxis]
    g_direction_a = g_a * g_normalization_factors
    #The directions toward the center
    return g_direction_a


def set_centers_on_surface__kdtree_v1(iobj, centroids, average_edge, grid, grid_values, debug_vf=None, more_normals=None, more_cv_pairs=[]):
    """ see set_centers_on_surface__ohtake_v3s_002().
    Note: zero means |v| <= th  (note: equality) """
    print "Projecting the centroids: newer age"

    THRESHOLD_zero_interval = 0.0001  # f == TH is NOT zero.
    MAX_ITER = 20
    USE_MESH_NORMALS = True
    EXTREME_ALPHA = False  # To use alpha > 1. that exceeds (violates) max_dim

    max_dist = average_edge

    x = centroids

    fc_a = iobj.implicitFunction(x)

    #visualise_scalar_distribution([f_a, taubin])
    #fc_a = f_a

    g_direction_a = make_grad_directions(iobj, x)

    signs_c = (fc_a > THRESHOLD_zero_interval)*1. - (fc_a < -THRESHOLD_zero_interval)*1.
    x0_v3 = x[:, :3]

    dx0_c_grad = - g_direction_a * signs_c[:, np.newaxis]

    step_size = max_dist * 1.

    alpha_list = []
    assert step_size > 0.001
    while step_size > 0.001:
        step_size = step_size * 0.5
        max_step = min(MAX_ITER, int(math.floor(max_dist/math.fabs(step_size) + 0.001)) )
        assert max_step >= 2  # at least one step
        #if max_step
        #violated only at first time but the first point is already done.
        for i in range(1, max_step+1, 2): #Step size is two, to avoid aready visited points
            alpha = float(i)*step_size
            # print i, alpha/average_edge
            alpha_list += [alpha/average_edge]
            alpha_list += [-alpha/average_edge]
            # alpha is prepared

    if EXTREME_ALPHA:
        # This violates the max_dim condition. This is a bad idea. Set to False.
        alpha_list += [+1., -1., +1.5, -1.5, +2., -2.]

    print "Alphas:"
    for i in range(len(alpha_list)):
        print "[%d]%g"%(i, alpha_list[i]),
    print
    print

    # The algorithm
    n = x0_v3.shape[0]
    best_result_x = np.ones((n, 3))
    active_indices = np.arange(0, n, dtype=int)
    active_count = n
    del n

    still_nonsuccess_indices = active_indices

    print "points: ", active_count, ".",
    already_success = fc_a*0 > 1.  # all False
    success = already_success.copy()  # falses  #.copy() is necessary
    assert not np.any(already_success)

    if USE_MESH_NORMALS:
        assert mesh_normals is not None
        dx0c_mesh_normals = mesh_normals
        assert np.allclose(np.linalg.norm(dx0c_mesh_normals, axis=1), 1.)

    TEST = True  # and not optimised_used():

    #print "left(found)",
    #for it in [1, 0, 2, 3, 4,5,6] if USE_MESH_NORMALS else [0]:
    for it in [0, 1, 2, 3, 4,5,6] if USE_MESH_NORMALS else [0]:

        if it == 0:
            dxc = dx0_c_grad
            alpha_list1 = alpha_list
        elif it == 1:
            dxc = dx0c_mesh_normals.copy() #* 0.5
            print
            print "now mesh normals"
            alpha_list1 = alpha_list[:10]
        elif it == 2:
            n = dx0_c_grad.shape[0]
            r = 0.000001
            perturb = (np.random.rand(n, 3)*2.-1.) * r
            #set_trace()
            z = np.cross(dx0c_mesh_normals, dx0_c_grad+perturb, axis=1)
            z = z / np.linalg.norm(z, axis=1, keepdims=True)
            assert np.allclose(np.linalg.norm(z, axis=1), 1.)
            global z3
            z3 = z
            #set_trace()
            #dxc[:, :] = z
            dxc = z.copy()

        elif it == 3:
            n = dx0_c_grad.shape[0]
            m = dx0c_mesh_normals
            missing = np.nonzero(np.linalg.norm(m, axis=1) < 1.-0.00001)[0]
            m[missing] = np.random.randn(missing.shape[0], 3)
            z2 = np.cross(m, z, axis=1)
            z2 = z2 / np.linalg.norm(z2, axis=1, keepdims=True)
            #set_trace()
            assert np.allclose(np.linalg.norm(z2, axis=1), 1.), "zero or NaN vectors"
            dxc = z2.copy()

        elif it in [4, 5, 6]:
            z3 = z.copy()
            dim_i = it-4
            print
            print "XYZ"[dim_i], "direction:"
            z3[:, :] = 0.
            z3[:, dim_i] = 1.
            #oops!!!
            dxc = z3
        else:
            print "Error"

        seq.append()

    for ii in seq:
        counter = -1
        for alpha in alpha_list1:
            counter += 1
            x1_half = x0_v3 + (max_dist*alpha)*dxc
            FAST = True
            if FAST:
                active_indices = still_nonsuccess_indices
                #set_trace()

                # Todo: For those that have changed sign, check if they are closer actually.
                xa4 = augment4(x1_half[active_indices, :])
                f_a = iobj.implicitFunction(xa4)
                signs_a = (f_a > THRESHOLD_zero_interval)*1. + (f_a < -THRESHOLD_zero_interval)*(-1.)
                #success = signs_a * signs_c <= 0.
                success0 = signs_a * signs_c[active_indices] <= 0.
                success[:] = False
                #success[success0] = True
                assert np.all(success == False)
                assert np.all(success[active_indices] == False)
                success[active_indices] = success0
                #print "success:", np.sum(success),
            else:
                # Todo: For those that have changed sign, check if they are closer actually.
                xa4 = augment4(x1_half)
                f_a = iobj.implicitFunction(xa4)
                signs_a = (f_a > THRESHOLD_zero_interval)*1. + (f_a < -THRESHOLD_zero_interval)*(-1.)
                success = signs_a * signs_c <= 0.
                #print "success:", np.sum(success),
            assert success.ndim == 1
            #print "success", np.sum(success)
            #print "already_success", np.sum(already_success)

            new_success_indices = np.nonzero(np.logical_and(success, np.logical_not(already_success)))[0]
            #already_success === old success
            #print "found ", new_success_indices.size,
            #print "(%d)"%new_success_indices.size,

            #collect zeros
            #efficient verison: only from new_success_indices
            #zeros_boolean = np.abs(f_a) <= THRESHOLD_zero_interval

            #result_x[active_indices[new_success_indices]]
            #active_indices = np.nonzero(np.logical_not(success))[0]
            #nonsuccess_indices = np.nonzero(success)[0]
            still_nonsuccess_indices = np.nonzero(np.logical_and(np.logical_not(success), np.logical_not(already_success)))[0]
            best_result_x[new_success_indices, :] = x1_half[new_success_indices, :]
            #print "new success>>", new_success_indices.size, "<<  ",
            #print "already>>", np.sum(already_success.size), "<<  ",
            #todo: also try som ein already_success and improve by replacing those that are CLOSER.
            #already_success_but_open_to_improvement = ...
            #best_so_far = ...

            #for next round
            already_success = np.logical_or(success, already_success)  # Union
            #print "left:", still_nonsuccess_indices.shape, ".",
            print ("[%d](+%d)%d "%(counter, new_success_indices.size, still_nonsuccess_indices.size,)),
            #print "already_success", np.sum(already_success)
            #active_indices = still_nonsuccess_indices

            if TEST:
                x2 = augment4(best_result_x[already_success, :])
                f2 = iobj.implicitFunction(x2)
                f2[np.abs(f2) < THRESHOLD_zero_interval] = 0.
                x1 = centroids[already_success, :]
                assert x1.shape[1] == 4
                f1 = iobj.implicitFunction(x1)
                f1[np.abs(f1) < THRESHOLD_zero_interval] = 0.
                #print f1
                #print f2
                bad_accepted =  (f1*f2)[f1 * f2 > 0]
                if bad_accepted.size > 0:
                        print "bad accepted f1*f2:", (f1*f2)[f1 * f2 > 0]
                #assert np.all((f1*f2)[f1 * f2 > 0])
                assert np.all(f1 * f2 <= 0)
                del f1, f2, x1, x2

            if still_nonsuccess_indices.shape[0] == 0:
                break

        if still_nonsuccess_indices.shape[0] == 0:
            break  # break the 'for' loop too

    # if still_nonsuccess_indices.shape[0] > 0:
    best_result_x[still_nonsuccess_indices, :] = x0_v3[still_nonsuccess_indices, :]  # failed to converge

    if TEST:
        xa1 = augment4(x0_v3)
        f1_test = iobj.implicitFunction(xa1)
        xa2 = augment4(best_result_x)
        f2_test = iobj.implicitFunction(xa2)
        s = f1_test*f2_test
        #print s[still_nonsuccess_indices]
        s[still_nonsuccess_indices] = -1.
        #print s[s > 0]
        assert np.all(s <= +THRESHOLD_zero_interval)  # can contain almost-zeros. Include the ==equality in zero-ness
        del s
        del f1_test
        del f2_test
        print "OK"

    #centroids[:, :3] = best_result_x[:, :]
    assert np.all(centroids[:, 3] == 1.)
    #print "."; flush()

    # ------------
    # Prepare for bisection: By removing zeros and moving negatives to x1 by swapping.

    #collect the zero ones
    xa1 = augment4(x0_v3)
    f1 = fc_a  # iobj.implicitFunction(xa1)
    assert np.all(f1 == iobj.implicitFunction(xa1), axis=None)
    xa2 = augment4(best_result_x)
    f2 = iobj.implicitFunction(xa2)
    zeros2_bool = np.abs(f2) <= THRESHOLD_zero_interval
    # Copy the zeros onto results
    zeros1_bool = np.abs(f1) <= THRESHOLD_zero_interval
    best_result_x[zeros1_bool, :3] = x0_v3[zeros1_bool, :3]
    zeros1or2 = np.logical_or(zeros1_bool, zeros2_bool)  # output
    del zeros1_bool, zeros2_bool
    assert np.all(np.abs(iobj.implicitFunction(augment4(best_result_x[zeros1or2, :]))) <= THRESHOLD_zero_interval)

    #ROOT_TOLERANCE = 0.000001
    ROOT_TOLERANCE = THRESHOLD_zero_interval  # because of the assert
    #relevants_bool = already_success[not.logical_not(zeros1or2[already_success])]
    #relevants_bool = np.delete(already_success, np.nonzero(zeros1or2)[0])
    relevants_bool = np.logical_and(already_success, np.logical_not(zeros1or2))
    #print np.nonzero(np.logical_not(already_success))[0]
    #print np.nonzero(relevants_bool)[0]
    #print np.nonzero(zeros1or2)[0]
    assert np.all(np.abs(f2[relevants_bool]) > +THRESHOLD_zero_interval)
    assert np.all(np.abs(f2[zeros1or2]) <= +THRESHOLD_zero_interval)

    #print "--"*10
    #print "relevants_bool:", np.sum(relevants_bool)
    #relevants_bool = np.logical_and(already_success, np.logical_not(zeros1or2))
    #set_trace()

    assert np.all( mysign_np(f2[relevants_bool], THRESHOLD_zero_interval) * mysign_np(f1[relevants_bool], THRESHOLD_zero_interval) < 0)
    del f1

    x2_relevant_v4 = augment4(best_result_x[relevants_bool, :])
    x1_relevant_v4 = centroids[relevants_bool, :]
    f1_relevants = iobj.implicitFunction(x1_relevant_v4)  # for assert only
    f2_relevants = iobj.implicitFunction(x2_relevant_v4)

    assert np.all(f1_relevants*f2_relevants <= +THRESHOLD_zero_interval)
    del f1_relevants

    #Swap negatives and positives
    swap_bool = f2_relevants < -THRESHOLD_zero_interval
    # performance: boolean or indices?
    temp = x2_relevant_v4[swap_bool, :]
    x2_relevant_v4[swap_bool, :] = x1_relevant_v4[swap_bool, :]
    x1_relevant_v4[swap_bool, :] = temp
    del temp
    #temp_f = f2_relevants[swap_bool]
    #f2_relevants[swap_bool] = f1_relevants[swap_bool]
    #f1_relevants[swap_bool] = temp_f
    #del temp_f
    del f2_relevants

    x_bisect = bisection_vectorized5_(iobj, x1_relevant_v4, x2_relevant_v4, ROOT_TOLERANCE)
    assert np.all(np.abs(iobj.implicitFunction(x_bisect) ) < THRESHOLD_zero_interval)
    assert x_bisect.shape[0] == np.sum(relevants_bool)
    centroids[relevants_bool, :] = x_bisect[:, :]  # x4
    #print np.nonzero(np.logical_not(relevants_bool))[0]
    #set_trace()
    print "centroids:", centroids.shape, "best_result_x:", best_result_x.shape, "relevants_bool:", np.sum(relevants_bool), "+", np.sum(np.logical_not(relevants_bool)), "zeros1or2=", np.sum(zeros1or2)
    centroids[zeros1or2, :3] = best_result_x[zeros1or2, :]  # x3

    #relevants_bool = np.logical_and(already_success, np.logical_not(zeros1or2))

    """
    if debug_vf is not None:
        _vs, _fs = debug_vf
        from visual5 import *
        c3 = centroids[zeros1or2, :3]
        #(np.zeros((0, 4)), np.zeros((0, 3), dtype=int))
        display_simple_using_mayavi_2([(_vs, _fs)],
                   pointcloud_list=[centroids[zeros1or2, :]], pointsizes=[0.02], #pointcloud_list=[point_collector.get_as_array()], pointsizes=[0.01],
                   mayavi_wireframe=[False,], opacity=[0.4,],
                   gradients_at=c3,
                   #separate_panels=False,
                   gradients_from_iobj=iobj,
                   #minmax=(RANGE_MIN, RANGE_MAX),
                   #add_noise=[0, 0], noise_added_before_broadcast=True,
                   labels=(centroids, zeros1or2), grad_arrow_len=0.2/2.)

        set_trace()
    """
    return
    #return zeros1or2

    #visualise_scalar_distribution([f_plot1, f_a])
    #visualise_scatter(f_plot1, f_a)
    #visualise_scatter(f_plot1, f_plot1*0.2)
    #exit()



