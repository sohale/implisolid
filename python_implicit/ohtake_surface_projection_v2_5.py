import numpy as np
VERBOSE = False
from basic_types import check_vector4_vectorized
FAST_VERSION = True
import math
from ipdb import set_trace
from utils import optimised_used
from utils import flush

# from surface_projection_v2_6 import set_centers_on_surface__ohtake_v3s_003


objname = "ell_example1" #"make_bricks"  # "cube_with_cylinders"


def visualise_scalar_distribution(scalar_array_list):
    STEPSIZE = 0.2
    import math
    import matplotlib.pyplot as plt
    for scalar_array in scalar_array_list:
        print scalar_array.shape
        assert scalar_array.ndim == 1
        #plt.hist( scalar_array , 150)
        weights = np.ones_like(scalar_array)/float(scalar_array.size)
        n, bins, patches = plt.hist( scalar_array , 50, weights=weights, alpha=0.5)  # normed=1

    special_lengths = np.array([STEPSIZE,], dtype=float)
    plt.plot(special_lengths, special_lengths*0+0.1, "*")
    plt.title(objname)
    plt.show()


def visualise_scatter(xa, ya):
    STEPSIZE = 0.2
    NOISE = 0.003*0
    xa = xa + np.random.randn(xa.size)*NOISE
    ya = ya + np.random.randn(ya.size)*NOISE
    ss= STEPSIZE
    import math
    import matplotlib.pyplot as plt
    #matplotlib.style.use('ggplot')
    #matplotlib.use('Agg')

    #plt.subplot(1, 2, 1)
    plt.plot(xa, ya, ".")
    plt.plot([-ss, ss], [-ss, ss], "k--")
    plt.plot([-ss, ss], [ss, -ss], "k--")
    special_lengths = np.array([STEPSIZE,], dtype=float)
    #plt.plot(special_lengths, special_lengths*0+0.1, "*")
    plt.title(objname)
    plt.xlabel('f1')
    plt.ylabel('f2')
    plt.grid(True)

    plt.axes().set_aspect('equal', 'datalim')


    if False:
        plt.subplot(1, 2, 2)

        sx = np.sign(xa)
        xa = np.log(np.abs(xa)) * sx

        sy = np.sign(ya)
        ya = np.log(np.abs(ya)) * sy
        ss = np.log(STEPSIZE)*5 /5

        plt.plot(xa, ya, ".")
        plt.plot([-ss, ss], [-ss, ss], "k--")
        plt.plot([-ss, ss], [ss, -ss], "k--")
        special_lengths = np.array([STEPSIZE,], dtype=float)
        #plt.plot(special_lengths, special_lengths*0+0.1, "*")
        plt.title(objname)
        plt.xlabel('f1')
        plt.ylabel('f2')
        plt.grid(True)
        #plt.axes().set_aspect('equal', 'datalim')
        #plt.axes().set_aspect('equal')


    plt.savefig('foo3.png', bbox_inches='tight', dpi=600.*2.)  #

    #plt.ioff()
    #plt.ion()
    #plt.close()
    plt.show()


def augment4(x):
    return np.concatenate((x, np.ones((x.shape[0], 1))), axis=1)


def set_centers_on_surface__ohtake_v3s_001(iobj, centroids, average_edge, nones_map):
    """ Visualises f(x) of points during evolution. see set_centers_on_surface__ohtake() """
    print "Projecting the centroids: new age"

    print "s",
    flush()

    THRESHOLD_minimum_gradient_len = 0.000001  # kill gradients smaller than this
    THRESHOLD_zero_interval = 0.0001

    x = centroids

    f_a = iobj.implicitFunction(x)
    g_a = iobj.implicitGradient(x)[:, :3]
    glen_a = np.linalg.norm(g_a, axis=1)
    glen_a[np.abs(glen_a) < THRESHOLD_minimum_gradient_len] = 1.

    #taubin = f_a/glen_a
    #visualise_scalar_distribution([f_a, taubin])
    f_plot1 = f_a

    g_normalization_factors = 1. / glen_a[:, np.newaxis]
    g_normalized_a = g_a * g_normalization_factors
    g_normalized_a_ = g_normalized_a
    #The directions toward the center

    step_size = average_edge / 2. * 0.2

    #signs = 1.*(f_a > THRESHOLD_zero_interval) - 1.*(f_a < -THRESHOLD_zero_interval)
    signs = (f_a > THRESHOLD_zero_interval)*step_size - (f_a < -THRESHOLD_zero_interval)*step_size
    #todo: stop moving when hit zero
    x0_v3 = x[:, :3]
    #Move opposite the direction toward center if the value is positive.
    x1 = x0_v3 - g_normalized_a * signs[:, np.newaxis] #* step_size

    x1_ = augment4(x1)
    f_a = iobj.implicitFunction(x1_)

    g_a = iobj.implicitGradient(x1_)[:, :3]
    glen_a = np.linalg.norm(g_a, axis=1)
    glen_a[np.abs(glen_a) < THRESHOLD_minimum_gradient_len] = 1.
    g_normalization_factors = 1. / glen_a[:, np.newaxis]
    g_normalized_a = g_a * g_normalization_factors

    print "."
    flush()

    #visualise_scalar_distribution([f_plot1, f_a])
    visualise_scatter(f_plot1, f_a)
    #visualise_scatter(f_plot1, f_plot1*0.2)

    exit()


from utils import Timer

@profile
def set_centers_on_surface__ohtake_v3s_002(iobj, centroids, average_edge, nones_map, debug_vf=None, mesh_normals=None):
    """ see set_centers_on_surface__ohtake() """
    print "Projecting the centroids: new age"

    #print "s", ; flush()

    THRESHOLD_minimum_gradient_len = 0.000001  # kill gradients smaller than this
    THRESHOLD_zero_interval = 0.0001  # f == TH is NOT zero.
    MAX_ITER = 20
    USE_MESH_NORMALS = True
    EXTREME_ALPHA = False  # To use alpha > 1. that exceeds (violates) max_dim

    max_dist = average_edge

    x = centroids

    fc_a = iobj.implicitFunction(x)
    g_a = iobj.implicitGradient(x)[:, :3]
    glen_a = np.linalg.norm(g_a, axis=1)
    glen_a[np.abs(glen_a) < THRESHOLD_minimum_gradient_len] = 1.

    #visualise_scalar_distribution([f_a, taubin])
    #fc_a = f_a

    g_normalization_factors = 1. / glen_a[:, np.newaxis]
    g_direction_a = g_a * g_normalization_factors
    #The directions toward the center
    del g_a

    #step_size = max_dist / 2.  # average_edge / 2.

    negative_f_c = fc_a < -THRESHOLD_zero_interval
    g_direction_a[negative_f_c, :] = -g_direction_a[negative_f_c, :]
    # The only difference here is the correct direction is tested "first".
    #Issue remains: The "distance-after-bisection" should be compared. But here, we just compare before the subdivision.

    #signs = 1.*(f_a >= THRESHOLD_zero_interval) - 1.*(f_a <= -THRESHOLD_zero_interval)
    #signs_c = (fc_a >= THRESHOLD_zero_interval)*step_size - (fc_a <= -THRESHOLD_zero_interval)*step_size
    signs_c = (fc_a > THRESHOLD_zero_interval)*1. - (fc_a < -THRESHOLD_zero_interval)*1.
    #todo: stop moving when hit zero
    x0_v3 = x[:, :3]   # it's a const
    # x is not used later anymore.
    del x
    #Move opposite the direction toward center if the value is positive.

    #del step_size

    dx0_c_grad = - g_direction_a * signs_c[:, np.newaxis]

    directions_basedon_gradient = dx0_c_grad
    del dx0_c_grad


    """
    f_a = fc_a # ???
    taubin = f_a/glen_a
    x1_taubin = x0_v3 - g_direction_a * taubin[:, np.newaxis]
    x1_half = x0_v3 + 0.5*directions_basedon_gradient * step_size
    x1_half_opposite = x0_v3 - 0.5*directions_basedon_gradient * step_size
    #Opposite search: Ohtake does not loop the opposite direction if it did not find te point in the forward direction.
    # boundary:
    ...
    candidates = [x1]  # [x1, x1_taubin, x1_half, x1_half_opposite]

    for xa in candidates:
        xa4 = augment4(xa)
        f_a = iobj.implicitFunction(xa4)
        signs = (f_a > THRESHOLD_zero_interval)*1. - (f_a < -THRESHOLD_zero_interval)*1.
        signs = signs*step_size
        #zeros, negatives.
        success = f_a * fc_a <= 0.  # May miss forget about accidental zeros. In that case, use the slightely further points.
        assert success.ndim == 1
    """
    step_size = max_dist / 2. * 2.

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

    """ Interesting observation: alpha_list[4] always finds many points. We can bring it forward, after [1] or [2] """
    """ If we use mesh normals first, the results seem better! """
    #lanp = np.array(alpha_list)
    #lanp.sort()
    #print lanp

    #print np.diff(lanp) * (2**9)
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

    # USE_MESH_NORMALS seems not necessary
    if USE_MESH_NORMALS:
        assert mesh_normals is not None
        dx0c_mesh_normals = mesh_normals
        assert np.allclose(np.linalg.norm(dx0c_mesh_normals, axis=1), 1.)

    # !!!!!!!!!!!!!!!
    # Was left unnecessary testing !!
    TEST = not optimised_used()

    # todo: move definition of alpha_list here

    #print "left(found)",
    #for it in [1, 0, 2, 3, 4,5,6] if USE_MESH_NORMALS else [0]:
    for it in [0, 1, 2, 3, 4,5,6] if USE_MESH_NORMALS else [0]:

        # alpha_list1 = alpha_list[:10]  // Is this the default value?
        if it == 0:
            dxc = directions_basedon_gradient
            alpha_list1 = alpha_list
        elif it == 1:
            dxc = dx0c_mesh_normals.copy() #* 0.5
            print
            print "now mesh normals"
            alpha_list1 = alpha_list[:10]
        elif it == 2:
            n = directions_basedon_gradient.shape[0]
            r = 0.000001
            perturb = (np.random.rand(n, 3)*2.-1.) * r
            #set_trace()
            z = np.cross(dx0c_mesh_normals, directions_basedon_gradient+perturb, axis=1)
            z = z / np.linalg.norm(z, axis=1, keepdims=True)
            assert np.allclose(np.linalg.norm(z, axis=1), 1.)
            global z3
            z3 = z
            #set_trace()
            #dxc[:, :] = z
            dxc = z.copy()

        elif it == 3:
            n = directions_basedon_gradient.shape[0]
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
            dxc = z3  # !!missed!
        else:
            print "Error"

        counter = -1
        for alpha in alpha_list1:
            counter += 1
            x1_half = x0_v3 + (max_dist*alpha)*dxc
            FAST = True
            if FAST:
                active_indices = still_nonsuccess_indices
                #set_trace()

                # Todo: For those that have changed sign, check if they are closer actually.
                with Timer() as t1:
                    xa4 = augment4(x1_half[active_indices, :])
                    f_a = iobj.implicitFunction(xa4)
                print "Evaluating on %d points"%xa4.shape[0], "took %g sec"%(t1.interval)
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

    """
    Big todo:
        This part is a mess. refactor into separate functions so that the data flow is known and variable names are un-entangled.
    """
    if TEST:
        # TEST IF all elements in x0_v3 & best_result_x are DUAL (CONJUGATE)
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

    # fix me: Is this a bug? The multiplication needs to be compared with: <= ROOT_TOLERANCE * ROOT_TOLERANCE, not just ROOT_TOLERANCE
    assert np.all(f1_relevants*f2_relevants <= +THRESHOLD_zero_interval)  # todo: Use THRESHOLD_zero_interval**2
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

# set_centers_on_surface__ohtake_v3s = set_centers_on_surface__ohtake_v3s_002
#set_centers_on_surface__ohtake_v3s = set_centers_on_surface__ohtake_v3s_003

def mysign_np(v, ROOT_TOLERANCE):
    return np.sign(v) * (np.abs(v) > ROOT_TOLERANCE)

#Put fastest version of bisection here:

def bisection_vectorized5_(iobj, x1_arr, x2_arr, ROOT_TOLERANCE):
    """ based on bisection_vectorized5. Note that this functin assumes there is no root in x1 and x2."""
    check_vector4_vectorized(x1_arr)
    check_vector4_vectorized(x2_arr)
    assert x1_arr.shape[0] == x2_arr.shape[0]
    v1_arr = iobj.implicitFunction(x1_arr)
    x2_arr[:, 3] = 1
    v2_arr = iobj.implicitFunction(x2_arr)

    result_x_arr = np.ones(x1_arr.shape)

    EPS = 0.000001  # sign

    n = x1_arr.shape[0]
    active_indices = np.arange(0, n)  # mid
    active_count = n
    solved_count = 0

    x_mid_arr = np.ones((active_count, 4))
    v_mid_arr = np.zeros((active_count,))

    iteration = 1
    while True:
        #print "iteration", iteration
        assert np.all(mysign_np(v2_arr[:active_count], ROOT_TOLERANCE) * mysign_np(v1_arr[:active_count], ROOT_TOLERANCE) < 0 - EPS)  # greater or equal
        assert np.all(v1_arr[:active_count] < 0-ROOT_TOLERANCE)
        assert active_indices.shape[0] == x1_arr[:active_count].shape[0]
        assert active_indices.shape[0] == x2_arr[:active_count].shape[0]
        x_mid_arr[:active_count] = ( x1_arr[:active_count] + x2_arr[:active_count] ) / 2.0
        v_mid_arr[:active_count] = iobj.implicitFunction(x_mid_arr[:active_count, :])
        assert active_indices.shape == (active_count,)
        assert active_indices.ndim == 1
        abs_ = np.abs(v_mid_arr[:active_count])
        indices_boundary = np.nonzero(abs_ <= ROOT_TOLERANCE)[0]  #eq
        indices_outside = np.nonzero(v_mid_arr[:active_count] < -ROOT_TOLERANCE)[0]  # gt
        indices_inside  = np.nonzero(v_mid_arr[:active_count] > +ROOT_TOLERANCE)[0]  # -v_mid_arr <  ROOT_TOLERANCE
        indices_eitherside = np.nonzero(abs_ > ROOT_TOLERANCE)[0]
        assert indices_boundary.size + indices_inside.size + indices_outside.size == active_count
        assert indices_eitherside.size + indices_boundary.size == active_count
        which_zeroed = active_indices[ indices_boundary ] # new start = mid
        found_count = indices_boundary.shape[0]
        solved_count += found_count
        assert active_count-found_count+solved_count == n
        result_x_arr[which_zeroed] = x_mid_arr[indices_boundary]
        assert np.all(indices_boundary < active_count)
        v2_arr[indices_inside] = v_mid_arr[indices_inside]
        x2_arr[indices_inside] = x_mid_arr[indices_inside]
        v1_arr[indices_outside] = v_mid_arr[indices_outside]
        x1_arr[indices_outside] = x_mid_arr[indices_outside]
        assert np.all(indices_outside < active_count)
        assert np.all(indices_inside < active_count)

        # ------ next round: --------
        assert active_count == active_indices.size
        active_indices = active_indices[indices_eitherside]
        assert active_count - found_count == active_indices.size
        old_active_count = active_count
        active_count = active_count - found_count
        assert active_count == indices_eitherside.size
        #again: does this hold again? assert active_count == active_indices.size
        iteration += 1
        assert np.all(indices_eitherside < old_active_count)
        v1_arr[:active_count] = v1_arr[indices_eitherside]
        v2_arr[:active_count] = v2_arr[indices_eitherside]
        x1_arr[:active_count] = x1_arr[indices_eitherside]
        x2_arr[:active_count] = x2_arr[indices_eitherside]

        assert active_indices.shape == v1_arr[:active_count].shape
        assert active_indices.shape[0] == active_count

        del old_active_count

        assert len(active_indices) == active_count
        if len(active_indices) == 0:
            break

    assert active_indices.size == 0
    optimisation_used = optimised_used()
    if not optimisation_used:
        v_arr = iobj.implicitFunction(result_x_arr)
        assert np.all(np.abs(v_arr) < ROOT_TOLERANCE)
    return result_x_arr


#global bisection_vectorized_thingy
#bisection_vectorized_thingy = bisection_vectorized5_

def bisection_vectorized5_compact(iobj, x1_arr, x2_arr, ROOT_TOLERANCE):
    """ based on bisection_vectorized5"""
    check_vector4_vectorized(x1_arr)
    check_vector4_vectorized(x2_arr)
    assert x1_arr.shape[0] == x2_arr.shape[0]
    v1_arr = iobj.implicitFunction(x1_arr)
    x2_arr[:, 3] = 1
    v2_arr = iobj.implicitFunction(x2_arr)
    result_x_arr = np.ones(x1_arr.shape)
    n = x1_arr.shape[0]
    #not tested
    active_indices = np.arange(0, n, dtype=int)  # mid
    active_count = n
    x_mid_arr = np.ones((active_count, 4))
    v_mid_arr = np.zeros((active_count,))

    while active_count > 0:
        x_mid_arr[:active_count] = ( x1_arr[:active_count] + x2_arr[:active_count] ) / 2.0
        v_mid_arr[:active_count] = iobj.implicitFunction(x_mid_arr[:active_count, :])
        #assert active_indices.shape == (active_count,)
        abs_ = np.abs(v_mid_arr[:active_count])
        indices_boundary = np.nonzero(abs_ <= ROOT_TOLERANCE)[0]  #eq
        indices_outside = np.nonzero(v_mid_arr[:active_count] < -ROOT_TOLERANCE)[0]  # gt
        indices_inside  = np.nonzero(v_mid_arr[:active_count] > +ROOT_TOLERANCE)[0]  # -v_mid_arr <  ROOT_TOLERANCE
        indices_eitherside = np.nonzero(abs_ > ROOT_TOLERANCE)[0]

        result_x_arr[active_indices[indices_boundary]] = x_mid_arr[indices_boundary]

        v2_arr[indices_inside] = v_mid_arr[indices_inside]
        x2_arr[indices_inside] = x_mid_arr[indices_inside]
        v1_arr[indices_outside] = v_mid_arr[indices_outside]
        x1_arr[indices_outside] = x_mid_arr[indices_outside]
        # ------ next round: --------
        # assert active_count == active_indices.size
        # Room for performance improvement: active_indices is replaced. active_indices[:active_count] = active_indices[indices_eitherside]
        active_indices = active_indices[indices_eitherside]  # only used for storing in result_x_arr
        active_count -= indices_boundary.shape[0]  # found_count = indices_boundary.shape[0]
        assert active_count == indices_eitherside.size
        v1_arr[:active_count] = v1_arr[indices_eitherside]
        v2_arr[:active_count] = v2_arr[indices_eitherside]
        x1_arr[:active_count] = x1_arr[indices_eitherside]
        x2_arr[:active_count] = x2_arr[indices_eitherside]
    assert active_indices.size == 0  # == active_count?
    return result_x_arr



def set_centers_on_surface__ohtake(iobj, centroids, average_edge, nones_map):
    #nones_map = centroids[:,0]*0 < 100
    print "Projecting the centroids:"
    for i in range(centroids.shape[0]):
        print i,
        if i==873:
            #set_trace()
            pass
        e = average_edge
        lm = e/2.
        max_dist = e

        c = project_point_bidir_ohtake_v2(iobj, centroids[i,np.newaxis,:], lm, max_dist)
        if c is not None:
            centroids[i] = c
        else:
            nones_map[i] = True
    #print nones_map
    if np.any(nones_map):
        print "failed projections: = ", np.sum(nones_map)


def check_direction(search_direction):
    if search_direction is not None:
        check_vector4_vectorized(search_direction)
        if not math.fabs(np.linalg.norm(search_direction[:, :3])-1.) < 0.000001:
            print search_direction
            return False
        if not math.fabs(search_direction[:, 3]-1.) < 0.000001:
            print search_direction
            return False
    return True


def project_point_bidir_ohtake(iobj, start_x, lambda_val, max_dist ):
    #set_trace()
    """ max_dist is used.
    See # setCenterOnSurface """
    #""" lambda_val: step size"""
    #max_iter = 20  # config["max_iter"]
    check_vector4_vectorized(start_x)
    assert start_x.shape[0] == 1

    max_iter=20

    if FAST_VERSION:
        p1 = search_near__ohtake_max_dist(iobj, start_x, None, lambda_val, max_iter, max_dist )
    else:
        p1 = search_near__ohtake(iobj, start_x, None, lambda_val, max_iter)
    if p1 is None:
        return None  # Should we return nothing if nothing found??
    f1 = iobj.implicitFunction(p1)

    # Mirror image: search the opposite way and accept only if it is closer than the best found.
    p2 = 2*start_x - p1
    f2 = iobj.implicitFunction(p2)
    p = p1  #None #p1 # None #p1  #default
    #see : setCenterOnSurface
    if f1*f2 < 0:  #WHY 'if' ???????????
        direction = (start_x - p1)  # as in Ohtake
        dn = np.linalg.norm(direction[:, :3])
        if dn > 0:  #dn>0.000000001:
            direction[:, :3] = direction[:, :3]/dn
            #check_vector4_vectorized(direction)
            #assert math.fabs(np.linalg.norm(direction[:, :3])-1.) < 0.000001
            #assert math.fabs(direction[0, 3]-1.) < 0.000001
            assert check_direction(direction)

            #broken
            if FAST_VERSION:
                p3 = search_near__ohtake_max_dist(iobj, start_x, direction, lambda_val, max_iter, max_dist)
            else:
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


def project_point_bidir_ohtake_v2(iobj, start_x, lambda_val, max_dist ):
    #set_trace()
    """ max_dist is used. lambda_val is the step size.
    See # setCenterOnSurface """
    check_vector4_vectorized(start_x)
    assert start_x.shape[0] == 1
    max_iter = 20  # todo: config["proj:max_iter"]

    if FAST_VERSION:
        p1 = search_near__ohtake_max_dist(iobj, start_x, None, lambda_val, max_iter, max_dist )
    else:
        p1 = search_near__ohtake(iobj, start_x, None, lambda_val, max_iter)
    #if p1 is None:
    #    return None  # Should we return nothing if nothing found??
    found1 = p1 is not None
    if not found1:
        p1 = start_x
    f1 = iobj.implicitFunction(p1)

    # Mirror image: search the opposite way and accept only if it is closer than the best found.
    p2 = 2*start_x - p1
    f2 = iobj.implicitFunction(p2)
    p = p1  #None #p1 # None #p1  #default
    #see : setCenterOnSurface

    f0 = iobj.implicitFunction(start_x)

    if f0*f2 < 0:  # if f1*f2 < 0:   # WHY 'if' ???????????  NO
        search_direction = (start_x - p1)  # as in Ohtake
        search_direction[:, 3] = 1.
        dn = np.linalg.norm(search_direction[:, :3])
        if dn > 0:  #dn>0.000000001:
            search_direction[:, :3] = search_direction[:, :3]/dn  #fixed
            assert check_direction(search_direction)

            #broken#
            if FAST_VERSION: #
                p3 = search_near__ohtake_max_dist(iobj, start_x, search_direction, lambda_val, max_iter, max_dist)
            else:
                p3 = search_near__ohtake(iobj, start_x, search_direction, lambda_val, max_iter)
            #no max_dist #
#
            if p3 is not None: #
                if np.linalg.norm(start_x[:, :3] - p3[:, :3]) > np.linalg.norm(start_x[:, :3] - p1[:, :3]): #
                    p = p3 #
                ##else:
                ##    p = p1
    ##else:
    ##    p = p1
    if p is not None:
        if np.linalg.norm(start_x[:, :3] - p[:, :3]) > max_dist:
            return None
    return p


def search_near__ohtake_max_dist(iobj, start_x, search_direction, lambda_val, MAX_ITER, max_dist):
    """Returns either the point, or None, if not found. lambda_val is the expected distance from the surface.
    The resommended value is half of the average edge length, but a better way is half of the MC'step size (because the expected error is half of the MC grid voxel size).
    Remeber: we may be on a totally irrelevant search_direction here.
    'search_direction' should be normalised. lambda*search_direction is used. Note that lambda is negated by default.
    Note: along_1d mode is in fact the same. It's just initialises search_direction=gradient(start_x).
    Does both searchNearPoint1D and searchNearPoint()
    :param search_direction: description
    @param np.array search_direction
    """
    #lambda should be ~ expected distance?  (that's why it should be half of the average edge size, i.e. half of the MC step size)

    TH1 = 0.001
    # MAX_ITER = 20
    #TH2_L = 0.00000001  # Used by Ohtake. But too small
    #TH2_L = 0.00001  #only used in the along_1d mode
    TH2_L = 0.1/2. #fast version
    #TH2_L = 0.00001

    if search_direction is not None:
        along_1d = True
    else:
        along_1d = False

    assert check_direction(search_direction)

    if not along_1d:
        direction_fixed = iobj.implicitGradient(start_x)  ## what?! why start_x ??
        dn = np.linalg.norm(direction_fixed[:, 0:3])
        assert direction_fixed.shape[0] == 1
        if dn>0.0:  # 00000001:
            direction_fixed[:, 0:3] = direction_fixed[:, 0:3]/dn #again?
        else:
                pass  # Finding is not going to happen. But it's fine.

    if not along_1d:
        assert check_direction(direction_fixed)

    eval_count = 0

    p1 = start_x
    f1 = iobj.implicitFunction(p1)
    eval_count += 1
    f0 = f1

    if math.fabs(f1) < TH1:
        return p1

    p2 = p1  # no need actually

    negative_f1 = -1 if f1 < 0. else +1
    # -1 => I am negative, searching for positive.
    # +1 => I am positive, searching for negative.
    lambda_ = lambda_val * (-negative_f1)  # negation of input is bad practice
    del lambda_val
    del negative_f1

    prev_start_x = None

    exit_A = False
    while True:
        max_iter_ = max(min(MAX_ITER, int(math.ceil(max_dist/math.fabs(lambda_) + 0.001)) ), 2)
        #print max_iter_
        assert max_iter_ >= 2
        #max_iter_ =MAX_ITER
        #print MAX_ITER, max_iter_, " > "*100

        # (C) jumps back here.
        for j in range(max_iter_):
            #search_direction = iobj.implicitGradient(start_x)  # what?!
            #dn = np.linalg.norm(search_direction)
            #if dn>0.0:  # 00000001:
            #    search_direction = search_direction/dn
            if not along_1d:
                #if prev_start_x is not None:
                #    assert np.linalg.norm(prev_start_x - start_x) == 0
                #prev_start_x = start_x.copy()
                #search_direction = iobj.implicitGradient(start_x)  ## what?! why start_x ??
                #dn = np.linalg.norm(search_direction[0:3])
                #if dn>0.0:  # 00000001:
                #    search_direction = search_direction/dn
                #else:
                #    pass  # Finding is not going to happen. But it's fine.
                #assert np.linalg.norm(search_direction - direction_fixed) == 0
                search_direction = direction_fixed
                #print "reused"
                #print "="*100
            else:
                pass

            assert check_direction(search_direction)

            p2 = p2 + lambda_ * search_direction
            p2[:, 3] = 1
            f2 = iobj.implicitFunction(p2)
            eval_count += 1


            # If sign changed.
            if f1*f2 < 0.0:  # If I am negative (f1<0), but f2 changes sign.
                # (A)
                exit_A = True  # success
                break

            #todo: if f got worse, exit_C
            #Sign has not changed:
            if math.fabs(f2) > math.fabs(f1):
                print "got worse!"
                assert f1*f2 > 0.0

            p1 = p2

        else:
            #for loop ended becasue of MAX_ITER

            #(C): next iteration with adaptively decreased lambda_. Revert and start over again using a smaller lambda
            lambda_ = lambda_ / 2.
            #either quit:
            if np.abs(lambda_) < TH2_L:
                print "(B)", eval_count
                return None   # ultimate failure.
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
        else:
            pass  # failed. Next round.
    #(A)
    assert f1*f2 < 0.0
    if f1 > 0:
        (p1, p2) = (p2, p1)
        (f1, f2) = (f2, f1)

    # f_outisde, f_inside

    assert f1 < 0
    assert f2 > 0
    assert math.fabs(f1) < 1000, str(p1)
    assert math.fabs(f2) < 1000, str(p2)

    converged, p1, p2, iter_log = bisection_prop_2(iobj, p1, p2, f1, f2, MAX_ITER/2)
    assert f1 < 0
    assert f2 > 0

    #print "eval_count", eval_count
    if converged:
        assert p2 is None
        return p1
    else:
        return None


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

    #if not FAST_VERSION:
    for j in range(MAX_ITER): #Reducing this does not make it faster. Removing prop_bisec will make it slower.

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


def bisection_3_standard(iobj, p1, p2, f1, f2, MAX_ITER):
    TH1 = 0.001
    TH3 = 0.001

    assert p1.shape[0] == 1
    assert p2.shape[0] == 1

    assert f1 < 0
    assert f1*f2 < 0, "Opposite signs required"
    for j in range(MAX_ITER):
        assert np.allclose(p1[:, 3], p2[:, 3])
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



























#########################################################
# A cleaner version, same logic.
#########################################################







def MAKE_ALPHA_LIST(step_size, max_dist, average_edge, MAX_ITER, EXTREME_ALPHA):
    alpha_list = []
    min_step_size = 0.001
    unit_of_length = average_edge
    assert step_size > min_step_size
    while step_size > min_step_size:
        step_size = step_size * 0.5
        max_step = min(MAX_ITER, int(math.floor(max_dist/math.fabs(step_size) + 0.001)) )
        assert max_step >= 2  # at least one step
        #if max_step
        #violated only at first time but the first point is already done.
        for i in range(1, max_step+1, 2): #Step size is two, to avoid aready visited points
            alpha = float(i)*step_size
            # print i, alpha/unit_of_length
            alpha_list += [alpha/unit_of_length]
            alpha_list += [-alpha/unit_of_length]
            # alpha is prepared

    if EXTREME_ALPHA:
        # This violates the max_dim condition. This is a bad idea. Set to False.
        alpha_list += [+1., -1., +1.5, -1.5, +2., -2.]

    """ Interesting observation: alpha_list[4] always finds many points. We can bring it forward, after [1] or [2] """
    """ If we use mesh normals first, the results seem better! """


    print "Alphas:"
    for i in range(len(alpha_list)):
        print "[%d]%g"%(i, alpha_list[i]),
    print
    print

    return alpha_list

# TEST IF all elements in x0_v3 & best_result_x are DUAL (CONJUGATE)
def TEST_OPPOSITE_SIGNS(iobj, x1, x2, THRESHOLD_zero_interval):
    # x2 = augment4(best_result_x[already_success, :])
    f2 = iobj.implicitFunction(x2)
    f2[np.abs(f2) < THRESHOLD_zero_interval] = 0.
    #x1 = centroids[already_success, :]
    assert x1.shape[1] == 4
    f1 = iobj.implicitFunction(x1)
    f1[np.abs(f1) < THRESHOLD_zero_interval] = 0.
    bad_accepted =  (f1*f2)[f1 * f2 > 0]
    if bad_accepted.size > 0:
        print "bad accepted f1*f2:", (f1*f2)[f1 * f2 > 0]
    #assert np.all((f1*f2)[f1 * f2 > 0])
    assert np.all(f1 * f2 <= 0)
    del f1, f2, x1, x2



def TEST_CONJUGATE(iobj, x0_v3, best_result_x, THRESHOLD_zero_interval):
    xa1 = augment4(x0_v3)
    f1_test = iobj.implicitFunction(xa1)
    xa2 = augment4(best_result_x)
    f2_test = iobj.implicitFunction(xa2)
    s = f1_test*f2_test
    s[still_nonsuccess_indices] = -1.
    assert np.all(s <= +THRESHOLD_zero_interval)  # can contain almost-zeros. Include the ==equality in zero-ness
    del s
    del f1_test
    del f2_test
    print "OK"

def set_centers_on_surface__ohtake_v3s_003(iobj, centroids, average_edge, nones_map, debug_vf=None, mesh_normals=None):
    """ see set_centers_on_surface__ohtake() """
    print "Projecting the centroids: new version v3s_003 ****"


    THRESHOLD_minimum_gradient_len = 0.000001  # kill gradients smaller than this
    THRESHOLD_zero_interval = 0.0001  # f == TH is NOT zero.
    MAX_ITER = 20
    EXTREME_ALPHA = False  # To use alpha > 1. that exceeds (violates) max_dim

    max_dist = average_edge

    x = centroids

    fc_a = iobj.implicitFunction(x)
    g_a = iobj.implicitGradient(x)[:, :3]
    glen_a = np.linalg.norm(g_a, axis=1)
    glen_a[np.abs(glen_a) < THRESHOLD_minimum_gradient_len] = 1.

    g_normalization_factors = 1. / glen_a[:, np.newaxis]
    g_direction_a = g_a * g_normalization_factors
    del g_a

    negative_f_c = fc_a < -THRESHOLD_zero_interval
    g_direction_a[negative_f_c, :] = -g_direction_a[negative_f_c, :]
    # The only difference here is the correct direction is tested "first".
    #Issue remains: The "distance-after-bisection" should be compared. But here, we just compare before the subdivision.

    signs_c = (fc_a > THRESHOLD_zero_interval)*1. - (fc_a < -THRESHOLD_zero_interval)*1.
    x0_v3 = x[:, :3]   # it's a const
    del x
    #Move opposite the direction toward center if the value is positive.

    dx0_c_grad = - g_direction_a * signs_c[:, np.newaxis]
    directions_basedon_gradient = dx0_c_grad
    del dx0_c_grad

    step_size = max_dist / 2. * 2.

    alpha_list = MAKE_ALPHA_LIST(step_size, max_dist, average_edge, MAX_ITER, EXTREME_ALPHA)

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

    assert mesh_normals is not None
    dx0c_mesh_normals = mesh_normals
    assert np.allclose(np.linalg.norm(dx0c_mesh_normals, axis=1), 1.)

    # !!!!!!!!!!!!!!!
    # Was left unnecessary testing !!
    TEST = not optimised_used()

    # todo: move definition of alpha_list here

    for it in [0, 1, 2, 3, 4,5,6]:

        # alpha_list1 = alpha_list[:10]  // Is this the default value?
        if it == 0:
            dxc = directions_basedon_gradient
            alpha_list1 = alpha_list
        elif it == 1:
            dxc = dx0c_mesh_normals.copy() #* 0.5
            print
            print "now mesh normals"
            alpha_list1 = alpha_list[:10]
        elif it == 2:
            n = directions_basedon_gradient.shape[0]
            r = 0.000001
            perturb = (np.random.rand(n, 3)*2.-1.) * r
            #set_trace()
            z = np.cross(dx0c_mesh_normals, directions_basedon_gradient+perturb, axis=1)
            z = z / np.linalg.norm(z, axis=1, keepdims=True)
            assert np.allclose(np.linalg.norm(z, axis=1), 1.)
            global z3
            z3 = z
            #set_trace()
            #dxc[:, :] = z
            dxc = z.copy()

        elif it == 3:
            n = directions_basedon_gradient.shape[0]
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
            dxc = z3  # !!missed!
        else:
            print "Error"

        counter = -1
        for alpha in alpha_list1:
            counter += 1
            x1_half = x0_v3 + (max_dist*alpha)*dxc

            active_indices = still_nonsuccess_indices
            #set_trace()

            # Todo: For those that have changed sign, check if they are closer actually.
            with Timer() as t1:
                xa4 = augment4(x1_half[active_indices, :])
                f_a = iobj.implicitFunction(xa4)
            print "Evaluating on %d points"%xa4.shape[0], "took %g sec"%(t1.interval)
            signs_a = (f_a > THRESHOLD_zero_interval)*1. + (f_a < -THRESHOLD_zero_interval)*(-1.)
            #success = signs_a * signs_c <= 0.
            success0 = signs_a * signs_c[active_indices] <= 0.
            success[:] = False
            assert np.all(success == False)
            assert np.all(success[active_indices] == False)
            success[active_indices] = success0


            assert success.ndim == 1

            new_success_indices = np.nonzero(np.logical_and(success, np.logical_not(already_success)))[0]
            #collect zeros
            #efficient verison: only from new_success_indices

            still_nonsuccess_indices = np.nonzero(np.logical_and(np.logical_not(success), np.logical_not(already_success)))[0]
            best_result_x[new_success_indices, :] = x1_half[new_success_indices, :]

            #todo: also try som ein already_success and improve by replacing those that are CLOSER.

            #for next round
            already_success = np.logical_or(success, already_success)  # Union

            print ("[%d](+%d)%d "%(counter, new_success_indices.size, still_nonsuccess_indices.size,)),

            if TEST:
                TEST_OPPOSITE_SIGNS(iobj, centroids[already_success, :], best_result_x[already_success, :], THRESHOLD_zero_interval)

            if still_nonsuccess_indices.shape[0] == 0:
                break

        if still_nonsuccess_indices.shape[0] == 0:
            break  # break the 'for' loop too

    # if still_nonsuccess_indices.shape[0] > 0:
    best_result_x[still_nonsuccess_indices, :] = x0_v3[still_nonsuccess_indices, :]  # failed to converge

    if TEST:
        # TEST IF all elements in x0_v3 & best_result_x are DUAL (CONJUGATE)
        TEST_CONJUGATE(iobj, x0_v3, best_result_x, THRESHOLD_zero_interval)

    assert np.all(centroids[:, 3] == 1.)

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

    ROOT_TOLERANCE = THRESHOLD_zero_interval  # because of the assert
    relevants_bool = np.logical_and(already_success, np.logical_not(zeros1or2))
    assert np.all(np.abs(f2[relevants_bool]) > +THRESHOLD_zero_interval)
    assert np.all(np.abs(f2[zeros1or2]) <= +THRESHOLD_zero_interval)

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
    del f2_relevants

    x_bisect = bisection_vectorized5_(iobj, x1_relevant_v4, x2_relevant_v4, ROOT_TOLERANCE)
    assert np.all(np.abs(iobj.implicitFunction(x_bisect) ) < THRESHOLD_zero_interval)
    assert x_bisect.shape[0] == np.sum(relevants_bool)
    centroids[relevants_bool, :] = x_bisect[:, :]  # x4
    print "centroids:", centroids.shape, "best_result_x:", best_result_x.shape, "relevants_bool:", np.sum(relevants_bool), "+", np.sum(np.logical_not(relevants_bool)), "zeros1or2=", np.sum(zeros1or2)
    centroids[zeros1or2, :3] = best_result_x[zeros1or2, :]  # x3

    return

# old, slow, unimporoved
# set_centers_on_surface__ohtake_v3s = set_centers_on_surface__ohtake_v3s_002

#good but unrefactored
set_centers_on_surface__ohtake_v3s = set_centers_on_surface__ohtake_v3s_002

#refactored
#set_centers_on_surface__ohtake_v3s = set_centers_on_surface__ohtake_v3s_003
