import numpy as np
VERBOSE = False
from basic_types import check_vector4_vectorized
FAST_VERSION = True
import math
from ipdb import set_trace


def set_centers_on_surface__ohtake(iobj, centroids, average_edge, nones_map):
    #nones_map = centroids[:,0]*0 < 100
    print "Projecting the centroids:"
    for i in range(centroids.shape[0]):
        print i,
        if i==873:
            set_trace()
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
