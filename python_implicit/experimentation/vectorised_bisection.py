import sys
sys.path.append("..")

#import kernprof

from basic_types import check_vector4_vectorized
#from numerical_raycast_bisection_vectorized
#from numerical_utils import numerical_raycast_bisection_vectorized
import math


import __builtin__
try:
    __builtin__.profile
except AttributeError:
    def profile(func):
        return func
    __builtin__.profile = profile


USE_KERNPROF = False

ROOT_TOLERANCE = 0.000001


def mysign_np(v):
    return np.sign(v) * (np.abs(v) > ROOT_TOLERANCE)


# Versions:
# 1 non-vectorised
# 2 old vectorised
# 3 new vectorised (x1, x2), polished
# 4 avoiding gc.

def bisection_3_standard(iobj, p1, p2, f1, f2, MAX_ITER, ROOT_TOLERANCE=0.001, MIN_DISTANCE=0.001):
    TH1 = MIN_DISTANCE  # 0.001 mm
    TH3 = ROOT_TOLERANCE  # 0.001 (function value)

    assert p1.shape[0] == 1
    assert p2.shape[0] == 1

    assert f1 < 0
    assert f1*f2 < 0, "Opposite signs required"
    for j in range(MAX_ITER):
        if np.linalg.norm(p1-p2) < TH1:
            if j == 0:  # Why was it not needed before??
                p3 = 0.5 * (p1 + p2)
            return p3, j  # Yes it happens
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


def bisection_pointwise0_6(iobj, x1_arr_, x2_arr_, ROOT_TOLERANCE=ROOT_TOLERANCE):
    # A point-wise version (version 1)
    """ x1_arr_ must be outside and x2_arr_ must be inside the object.
    Then this function finds points x=x1_arr+(lambda)*(x2_arr_-x1_arr_) where f(x)=0 using the bisection method."""
    check_vector4_vectorized(x1_arr_)
    check_vector4_vectorized(x2_arr_)
    assert x1_arr_.shape[0] == x2_arr_.shape[0]
    result_x_arr = np.zeros(x1_arr_.shape)
    EPS = 0.000001  # sign

    n = x1_arr_.shape[0]
    for i in range(n):

        x1v = x1_arr_[np.newaxis, i, :].copy()
        x2v = x2_arr_[np.newaxis, i, :].copy()

        v1_arr = iobj.implicitFunction(x1v)[0]
        v2_arr = iobj.implicitFunction(x2v)[0]

        MAX_ITER = 20
        result_x_arr[i, :], iterations = bisection_3_standard(iobj, x1v, x2v, v1_arr, v2_arr, MAX_ITER)

    return result_x_arr

#bisection_pointwise1 = bisection_pointwise0

def bisection_pointwise1(iobj, x1_arr_, x2_arr_, ROOT_TOLERANCE=ROOT_TOLERANCE):
    # A point-wise version (version 1)
    """ x1_arr_ must be outside and x2_arr_ must be inside the object.
    Then this function finds points x=x1_arr+(lambda)*(x2_arr_-x1_arr_) where f(x)=0 using the bisection method."""
    check_vector4_vectorized(x1_arr_)
    check_vector4_vectorized(x2_arr_)
    assert x1_arr_.shape[0] == x2_arr_.shape[0]
    result_x_arr = np.zeros(x1_arr_.shape)
    EPS = 0.000001  # sign

    n = x1_arr_.shape[0]
    for i in range(n):

        x1v = x1_arr_[i, :].copy()
        x2v = x2_arr_[i, :].copy()

        v1_arr = iobj.implicitFunction(x1v[np.newaxis, :])[0]
        v2_arr = iobj.implicitFunction(x2v[np.newaxis, :])[0]

        active_indices = np.arange(0, n)  # mid
        iteration = 1
        while True:
            assert mysign_np(v2_arr) * mysign_np(v1_arr) < 0 - EPS  # greater or equal
            assert v1_arr < 0-ROOT_TOLERANCE
            x_midpoint = (x1v[:] + x2v[:]) / 2.0
            x_midpoint[3] = 1
            v_mid_val = iobj.implicitFunction(x_midpoint[np.newaxis, :])[0]

            on_boundary_boolean = np.abs(v_mid_val) <= ROOT_TOLERANCE  #eq
            boolean_outside = v_mid_val < -ROOT_TOLERANCE  # gt
            boolean_inside  = v_mid_val > +ROOT_TOLERANCE  # -v_mid_val <  ROOT_TOLERANCE

            assert boolean_outside or boolean_inside == (not(on_boundary_boolean))
            if on_boundary_boolean:
                result_x_arr[i, :] = x_midpoint[:]
                break
            assert boolean_outside or boolean_inside
            if boolean_inside:
                v2_arr = v_mid_val
                x2v[:] = x_midpoint[:]
            if boolean_outside:
                v1_arr = v_mid_val
                x1v = x_midpoint

            iteration += 1

            assert x1v.shape == x2v.shape

    assert np.all(np.abs(iobj.implicitFunction(result_x_arr)) < ROOT_TOLERANCE)
    return result_x_arr




def numerical_raycast_bisection_vectorized(iobj, ray_x, ray_target, ROOT_TOLERANCE=ROOT_TOLERANCE):
    # My old vectorised implementation (version 2)
    """ ray_x must be outside and ray_target must be inside the object.
    This function the point x=ray_x+(lambda)*ray_n where f(x)=0, using the Bisection method."""
    ray_n = ray_target - ray_x
    ray_n[:, 3] = 1
    check_vector4_vectorized(ray_x)
    check_vector4_vectorized(ray_n)
    assert ray_x.shape[0] == ray_n.shape[0]
    x1_arr = ray_x  # start
    v1_arr = iobj.implicitFunction(x1_arr)
    x2_arr = x1_arr + ray_n * 1.0
    x2_arr[:,3] = 1
    v2_arr = iobj.implicitFunction(x2_arr)

    result_x_arr = np.zeros(ray_x.shape)

    EPS = 0.0000001  # sign

    n = x1_arr.shape[0]
    already_root = np.zeros((n,), dtype=np.int)
    #assert v2_arr * va > 0 - EPS  # greater or equal
    active_indices  = np.arange(0,n)  # mid
    iteration = 1
    while True:
        #print(v1_arr.shape, "*")
        #print(np.vstack((v1_arr,v2_arr)))
        #print(v2_arr * v1_arr )
        #print("aaaaaaaaaaaaaaaa")
        #print(v1_arr)
        #print(v2_arr)
        #print(v2_arr * v1_arr)
        #assert np.all(v2_arr * v1_arr < 0 - EPS)  # greater or equal
        assert np.all(mysign_np(v2_arr) * mysign_np(v1_arr) < 0 - EPS)  # greater or equal

        assert np.all(v1_arr < 0-ROOT_TOLERANCE)
        assert active_indices.shape[0] == x1_arr.shape[0]
        assert active_indices.shape[0] == x2_arr.shape[0]
        x_mid_arr = ( x1_arr + x2_arr ) / 2.0
        x_mid_arr[:,3] = 1
        v_mid_arr = iobj.implicitFunction(x_mid_arr)
        assert v_mid_arr.shape == active_indices.shape
        assert active_indices.ndim == 1

        #flipped_i = v_mid_arr
        #contains the indices
        dif = -v_mid_arr  # assuming x1 is always outside and x2 is inside
        assert dif.shape == active_indices.shape
        boolean_eq = np.abs(dif) <= ROOT_TOLERANCE
        boolean_gt =  dif >  ROOT_TOLERANCE
        boolean_lt =  dif <  -ROOT_TOLERANCE # dif <  ROOT_TOLERANCE
        boolean_neq = np.logical_not( boolean_eq )
        #logical_or
        assert np.all( np.logical_or(boolean_gt, boolean_lt) == np.logical_not(boolean_eq) )
        #print("boolean_gt", boolean_gt)  #t
        #print("boolean_lt", boolean_lt)  #f


        which_zeroed     = active_indices[ boolean_eq ] # new start = mid
        which_flippedAt1 = active_indices[ boolean_gt ] # new end = mid
        which_flippedAt2 = active_indices[ boolean_lt ]
        which_flippedAny = active_indices[ boolean_neq ]

        already_root[which_zeroed] = 1  # iteration
        result_x_arr[which_zeroed,:] = x_mid_arr[boolean_eq,:]


        #x1_arr and x2_arr should have the same size eventually. the boolean_eq should be removed from their indices.
        #the total is np.arange(n)
        v2_arr[boolean_lt] = v_mid_arr[boolean_lt]#[which_flippedAny]
        x2_arr[boolean_lt,:] = x_mid_arr[boolean_lt,:]#[which_flippedAny]   # which_flippedAt2

        #note: x2_arr is modified inplace (side-effect)


        #x1_arr and x2_arr both shrink here

        v1_arr[boolean_gt] = v_mid_arr[boolean_gt]#[which_flippedAny]
        x1_arr[boolean_gt,:] = x_mid_arr[boolean_gt,:]#[which_flippedAny]   #which_flippedAt1

        v1_arr = v1_arr[boolean_neq]
        v2_arr = v2_arr[boolean_neq]
        x1_arr = x1_arr[boolean_neq,:]
        x2_arr = x2_arr[boolean_neq,:]
        #print("active_indices = ", active_indices)
        #print("which_flippedAny = ", which_flippedAny)
        #print("active_indices[which_flippedAny] = ", active_indices[boolean_neq])
        active_indices = active_indices[boolean_neq] #which_flippedAt1 || which_flippedAt2
        iteration += 1

        assert x1_arr.shape == x2_arr.shape
        assert v1_arr.shape == v2_arr.shape
        assert active_indices.shape == v1_arr.shape

        #print(active_indices)
        #print(v1_arr, "****", v2_arr)
        #print(boolean_lt)
        #print("*******")

        if len(active_indices) == 0:
            break

    assert len(active_indices) == 0
    #result_x_arr
    v_arr = iobj.implicitFunction(result_x_arr)
    assert np.all(np.abs(v_arr) < ROOT_TOLERANCE)
    return result_x_arr

bisection_vectorized1 = numerical_raycast_bisection_vectorized


@profile
def bisection_vectorized2(iobj, x1_arr, x2_arr, ROOT_TOLERANCE=ROOT_TOLERANCE):
    # vectorised: polished. (version 3)
    """ x1_arr must be outside and x2_arr must be inside the object.
    ..."""
    check_vector4_vectorized(x1_arr)
    check_vector4_vectorized(x2_arr)
    assert x1_arr.shape[0] == x2_arr.shape[0]
    #x1_arr = x1_arr  # start
    v1_arr = iobj.implicitFunction(x1_arr)
    #x2_arr = x1_arr + ray_n * 1.0
    x2_arr[:, 3] = 1
    v2_arr = iobj.implicitFunction(x2_arr)

    result_x_arr = np.zeros(x1_arr.shape)

    EPS = 0.000001  # sign

    n = x1_arr.shape[0]
    #already_root = np.zeros((n,), dtype=np.int)
    #assert v2_arr * va > 0 - EPS  # greater or equal
    active_indices = np.arange(0, n)  # mid
    iteration = 1
    while True:
        #if not np.all(mysign_np(v2_arr) * mysign_np(v1_arr) < 0 - EPS):
        #    #print np.max(mysign_np(v2_arr) * mysign_np(v1_arr))
        #    #i = mysign_np(v2_arr) * mysign_np(v1_arr) < 0 - EPS
        #    #print np.concatenate((v2_arr[i, np.newaxis], v1_arr[i, np.newaxis]), axis=1)
        #    # #print np.max(v2_arr), np.min(v1_arr)
        #    #print np.min(v2_arr), np.max(v1_arr)
        #    #print mysign_np(np.min(v2_arr)), mysign_np(np.max(v1_arr))
        #    pass
        assert np.all(mysign_np(v2_arr) * mysign_np(v1_arr) < 0 - EPS)  # greater or equal
        assert np.all(v1_arr < 0-ROOT_TOLERANCE)
        assert active_indices.shape[0] == x1_arr.shape[0]
        assert active_indices.shape[0] == x2_arr.shape[0]
        x_mid_arr = ( x1_arr + x2_arr ) / 2.0
        x_mid_arr[:,3] = 1
        v_mid_arr = iobj.implicitFunction(x_mid_arr)
        assert v_mid_arr.shape == active_indices.shape
        assert active_indices.ndim == 1

        assert v_mid_arr.shape == active_indices.shape
        boolean_boundary = np.abs(v_mid_arr) <= ROOT_TOLERANCE  #eq
        boolean_outside = v_mid_arr < -ROOT_TOLERANCE  # gt
        boolean_inside  = v_mid_arr > +ROOT_TOLERANCE  # -v_mid_arr <  ROOT_TOLERANCE
        boolean_eitherside = np.logical_not(boolean_boundary)
        assert np.all(np.logical_or(boolean_outside, boolean_inside) == np.logical_not(boolean_boundary) )

        which_zeroed = active_indices[ boolean_boundary ] # new start = mid

        #already_root[which_zeroed] = 1  # iteration
        result_x_arr[which_zeroed,:] = x_mid_arr[boolean_boundary,:]


        #x1_arr and x2_arr should have the same size eventually. the boolean_boundary should be removed from their indices.
        #the total is np.arange(n)
        v2_arr[boolean_inside] = v_mid_arr[boolean_inside]
        x2_arr[boolean_inside,:] = x_mid_arr[boolean_inside,:]
        #note: x2_arr is modified

        #x1_arr and x2_arr both shrink here
        v1_arr[boolean_outside] = v_mid_arr[boolean_outside]
        x1_arr[boolean_outside,:] = x_mid_arr[boolean_outside,:]

        # ------ next round: --------

        #it does re-allocate
        v1_arr = v1_arr[boolean_eitherside]
        v2_arr = v2_arr[boolean_eitherside]
        x1_arr = x1_arr[boolean_eitherside,:]
        x2_arr = x2_arr[boolean_eitherside,:]
        active_indices = active_indices[boolean_eitherside]
        iteration += 1

        assert x1_arr.shape == x2_arr.shape
        assert v1_arr.shape == v2_arr.shape
        assert active_indices.shape == v1_arr.shape

        if len(active_indices) == 0:
            break

    assert len(active_indices) == 0
    v_arr = iobj.implicitFunction(result_x_arr)
    assert np.all(np.abs(v_arr) < ROOT_TOLERANCE)
    return result_x_arr


@profile
def bisection_vectorized_frozen4(iobj, x1_arr, x2_arr, ROOT_TOLERANCE=ROOT_TOLERANCE):
        # New vectorised implementation (version 4)
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
            assert np.all(mysign_np(v2_arr) * mysign_np(v1_arr) < 0 - EPS)  # greater or equal
            assert np.all(v1_arr < 0-ROOT_TOLERANCE)
            assert active_indices.shape[0] == x1_arr.shape[0]
            assert active_indices.shape[0] == x2_arr[:active_count].shape[0]

            x_mid_arr[:active_count] = ( x1_arr + x2_arr[:active_count] ) / 2.0

            v_mid_arr[:active_count] = iobj.implicitFunction(x_mid_arr[:active_count, :])
            assert active_indices.shape == (active_count,)
            assert active_indices.ndim == 1

            boolean_boundary = np.abs(v_mid_arr[:active_count]) <= ROOT_TOLERANCE  #eq
            boolean_outside = v_mid_arr[:active_count] < -ROOT_TOLERANCE  # gt
            boolean_inside  = v_mid_arr[:active_count] > +ROOT_TOLERANCE  # -v_mid_arr <  ROOT_TOLERANCE
            boolean_eitherside = np.logical_not(boolean_boundary)
            assert np.all(np.logical_or(boolean_outside, boolean_inside) == np.logical_not(boolean_boundary) )
            assert boolean_boundary.shape[0] == active_count

            which_zeroed = active_indices[ boolean_boundary ] # new start = mid
            found_count = np.sum(boolean_boundary)
            solved_count += found_count
            assert active_count-found_count+solved_count == n

            result_x_arr[which_zeroed, :3] = x_mid_arr[:active_count, :3][boolean_boundary[:], :3]

            v2_arr[boolean_inside] = v_mid_arr[:active_count][boolean_inside]
            x2_arr[:active_count][boolean_inside, :] = x_mid_arr[:active_count, :][boolean_inside,:]

            v1_arr[boolean_outside] = v_mid_arr[:active_count][boolean_outside]
            x1_arr[boolean_outside,:] = x_mid_arr[:active_count, :][boolean_outside,:]

            # ------ next round: --------
            assert active_count == active_indices.size
            active_indices = active_indices[boolean_eitherside]
            assert active_count - found_count == active_indices.size
            old_active_count = active_count
            active_count = active_count - found_count
            assert active_count == np.sum(boolean_eitherside)
            iteration += 1

            v1_arr = v1_arr[boolean_eitherside]
            v2_arr = v2_arr[boolean_eitherside]
            x1_arr = x1_arr[boolean_eitherside,:]
            x2_arr[:active_count] = x2_arr[:old_active_count][boolean_eitherside,:]

            assert active_count == v1_arr.shape[0]
            assert active_count == x1_arr.shape[0]
            assert active_count == v2_arr.shape[0]

            assert x1_arr.shape == x2_arr[:active_count].shape
            assert v1_arr.shape == v2_arr.shape
            assert active_indices.shape == v1_arr.shape
            assert active_indices.shape[0] == active_count

            del old_active_count

            if len(active_indices) == 0:
                break

        assert len(active_indices) == 0
        v_arr = iobj.implicitFunction(result_x_arr)
        assert np.all(np.abs(v_arr) < ROOT_TOLERANCE)
        return result_x_arr



@profile
def bisection_vectorized5(iobj, x1_arr, x2_arr, ROOT_TOLERANCE=ROOT_TOLERANCE):
    # New vectorised implementation (version 4)
    """ x1_arr must be outside and x2_arr must be inside the object.
    This function the point x=x1_arr+(lambda)*(x2-x1) where f(x)=0, using the Bisection method.
    Note: inplace operations will be done on x1_arr, x2_arr."""
    check_vector4_vectorized(x1_arr)
    check_vector4_vectorized(x2_arr)
    assert x1_arr.shape[0] == x2_arr.shape[0]
    #x1_arr = x1_arr  # start
    v1_arr = iobj.implicitFunction(x1_arr)
    #x2_arr = x1_arr + ray_n * 1.0
    x2_arr[:, 3] = 1
    v2_arr = iobj.implicitFunction(x2_arr)

    result_x_arr = np.ones(x1_arr.shape)

    EPS = 0.000001  # sign


    n = x1_arr.shape[0]
    #already_root = np.zeros((n,), dtype=np.int)
    #assert v2_arr * va > 0 - EPS  # greater or equal
    active_indices = np.arange(0, n)  # mid
    active_count = n
    solved_count = 0
    #invariants...

    x_mid_arr = np.ones((active_count, 4))
    v_mid_arr = np.zeros((active_count,))

    iteration = 1
    while True:
        #if not np.all(mysign_np(v2_arr) * mysign_np(v1_arr) < 0 - EPS):
        #    #print np.max(mysign_np(v2_arr) * mysign_np(v1_arr))
        #    #i = mysign_np(v2_arr) * mysign_np(v1_arr) < 0 - EPS
        #    #print np.concatenate((v2_arr[i, np.newaxis], v1_arr[i, np.newaxis]), axis=1)
        #    # #print np.max(v2_arr), np.min(v1_arr)
        #    #print np.min(v2_arr), np.max(v1_arr)
        #    #print mysign_np(np.min(v2_arr)), mysign_np(np.max(v1_arr))
        #    pass
        assert np.all(mysign_np(v2_arr[:active_count]) * mysign_np(v1_arr[:active_count]) < 0 - EPS)  # greater or equal
        assert np.all(v1_arr[:active_count] < 0-ROOT_TOLERANCE)
        assert active_indices.shape[0] == x1_arr[:active_count].shape[0]
        assert active_indices.shape[0] == x2_arr[:active_count].shape[0]

        x_mid_arr[:active_count] = ( x1_arr[:active_count] + x2_arr[:active_count] ) / 2.0
        #assert x_mid_arr.shape[0] == active_count
        #x_mid_arr[:active_count, 3] = 1

        #assert np.all(x_mid_arr[:active_count, 3] == 1)
        v_mid_arr[:active_count] = iobj.implicitFunction(x_mid_arr[:active_count, :])
        #assert v_mid_arr.shape[0] == active_count
        #assert v_mid_arr.shape == active_indices.shape
        #print active_indices.shape, active_count
        assert active_indices.shape == (active_count,)
        assert active_indices.ndim == 1

        ####assert v_mid_arr.shape == active_indices.shape
        abs_ = np.abs(v_mid_arr[:active_count])
        indices_boundary = np.nonzero(abs_ <= ROOT_TOLERANCE)[0]  #eq
        indices_outside = np.nonzero(v_mid_arr[:active_count] < -ROOT_TOLERANCE)[0]  # gt
        indices_inside  = np.nonzero(v_mid_arr[:active_count] > +ROOT_TOLERANCE)[0]  # -v_mid_arr <  ROOT_TOLERANCE
        indices_eitherside = np.nonzero(abs_ > ROOT_TOLERANCE)[0]
        #assert np.all(np.logical_or(indices_outside, indices_inside) == np.logical_not(indices_boundary) )
        #assert np.sorted(np.hstack((indices_outside, indices_inside)) == ?np.logical_not(indices_boundary) )
        #assert indices_boundary.shape[0] == active_count
        assert indices_boundary.size + indices_inside.size + indices_outside.size == active_count
        #print len(indices_eitherside), indices_eitherside.size, indices_eitherside.shape
        assert indices_eitherside.size + indices_boundary.size == active_count

        which_zeroed = active_indices[ indices_boundary ] # new start = mid
        found_count = indices_boundary.shape[0]
        solved_count += found_count
        assert active_count-found_count+solved_count == n
        #assert which_zeroed.shape[0] == active_count

        #already_root[which_zeroed] = 1  # iteration
        result_x_arr[which_zeroed] = x_mid_arr[indices_boundary]
        assert np.all(indices_boundary < active_count)
        #todo: x_mid_arr[indices_boundary[:active_count], :]

        #x1_arr and x2_arr should have the same size eventually. the boolean_boundary should be removed from their indices.
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
        iteration += 1

        #assert boolean_eitherside.size == old_active_count
        assert np.all(indices_eitherside < old_active_count)
        v1_arr[:active_count] = v1_arr[indices_eitherside]
        v2_arr[:active_count] = v2_arr[indices_eitherside]
        x1_arr[:active_count] = x1_arr[indices_eitherside]
        x2_arr[:active_count] = x2_arr[indices_eitherside]


        #assert active_count == v1_arr[:active_count].shape[0]
        #assert active_count == x1_arr[:active_count].shape[0]
        #assert active_count == v2_arr[:active_count].shape[0]
        #assert active_count == x2_arr.shape[0]

        #assert x1_arr[:active_count].shape == x2_arr[:active_count].shape
        #assert v1_arr[:active_count].shape == v2_arr[:active_count].shape
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



import numpy as np
global STEPSIZE

import sys

def flush():
    sys.stdout.flush()


def getifunc():
    from example_objects import make_example_vectorized
    ifunc = make_example_vectorized("cube_example")
    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-3/5., +5/5., 0.2*1.5/1.5  *2. /2.)
    return ifunc, (RANGE_MIN, RANGE_MAX, STEPSIZE)


def prepare():
    ifunc, (RANGE_MIN, RANGE_MAX, STEPSIZE) = getifunc()

    from stl_tests import make_mc_values_grid
    gridvals = make_mc_values_grid(ifunc, RANGE_MIN, RANGE_MAX, STEPSIZE, old=False)
    verts, facets = vtk_mc(gridvals, (RANGE_MIN, RANGE_MAX, STEPSIZE))
    return verts, facets

# verts, f_ = prepare(*getifunc())

ifunc, (RANGE_MIN, RANGE_MAX, STEPSIZE) = getifunc()

global xi
global xo
global maxn_pointwise
maxn_pointwise = 1000000*3 #/1000
n = 1000000*3 # / 1000
x = (np.random.rand(n, 4)*2-1)*(RANGE_MAX-RANGE_MIN)+RANGE_MIN
x[:, 3] = 1

f = ifunc.implicitFunction(x)

xi = x[f > ROOT_TOLERANCE*2, :]
xo = x[f < -ROOT_TOLERANCE*2, :]
print xi.shape, xo.shape
n_min = min(xi.shape[0], xo.shape[0])


#xio = xi[:n_min, :], xo[:n_min, :]
xi0 = xi[:n_min, :]
xo0 = xo[:n_min, :]
del xi, xo

#bisection_vectorized2

def testout(ifunc, q):
    f = ifunc.implicitFunction(q)
    #print np.max(np.fabs(f)), ROOT_TOLERANCE
    return True






def test2():
    global q2
    print "test2", ; flush()
    q2 = bisection_vectorized1(ifunc, xo, xi)
    assert testout(ifunc, q2)
    print "()", ; flush()

def test1():
    global maxn_pointwise
    if xo.shape[0]>maxn_pointwise:
        return

    global q1
    #print "test1", ; flush()
    q1 = bisection_pointwise1(ifunc, xo, xi)
    #print "()" ; flush()
    assert testout(ifunc, q1)


def test3():
    global q3
    #print "test3", ; flush()
    q3 = bisection_vectorized2(ifunc, xo, xi)
    #print "()" ; flush()
    assert testout(ifunc, q3)

    TEST_RESULTS = False
    if TEST_RESULTS:
        global q1
        global q2
        e1 = np.sum(np.fabs(q2-q1))
        e2 = np.sum(np.fabs(q3-q2))
        TOL = 0.00000001
        if not (e2 <TOL and e2 < TOL):
            print "error:",
            print e1,e2
        assert e1 < TOL
        assert e2 < TOL


def test5():
    global q5
    print "test5", ; flush()
    q5 = bisection_vectorized5(ifunc, xo, xi)
    print "()", ; flush()
    assert testout(ifunc, q5)


def test4():
    global q4
    q4 = bisection_vectorized_frozen4(ifunc, xo, xi)
    assert testout(ifunc, q4)


def test6():
    #Back to point-wise
    global q6
    print "test6", ; flush()
    q6 = bisection_pointwise0_6(ifunc, xo, xi)
    print "()", ; flush()
    assert testout(ifunc, q6)



def optimised_used():
    global _optimised_used
    _optimised_used = True
    def side_effect():
        global _optimised_used
        _optimised_used = False
        return True
    assert side_effect()
    #print "optimisation", _optimised_used
    return  _optimised_used


import matplotlib.pyplot as plt


def experiment1():
    #test3()
    #test2()
    # 20000
    # Best results: test5(): 4.5 (sec) for 300'000 points
    na = [1,2,5,10, 100,200,400,600, 800,  1000, 2000, 10000, ]  # 100000, 200000, 300000] #, 1000000]
    if USE_KERNPROF:
        na = [10000]
    global n_min
    na = filter(lambda e: e <= n_min, na)

    test_scripts = ['test1()', 'test2()', 'test3()', 'test4()', 'test5()', 'test6()']
    sty = {2: "r*-", 1: "bs-", 0: ".k-", 3: "m^-", 4: "r^--", 5: ".k:"}
    lbl = {2: 'numpy vec. (new)', 1: 'numpy vec. (old)', 0: 'point-wise', 3: 'vec-inplace (old)', 4: 'vec-inplace', 5: 'pointwise actually used'}  # vec-no-alloc
    COUNT = len(test_scripts)-1

    tl = []
    for ei in range(len(na)):
        n = na[ei]
        #x = (np.random.rand(n, 4)*2-1)*(RANGE_MAX-RANGE_MIN)+RANGE_MIN
        #x[:, 3] = 1
        global xi
        global xo
        assert xi0.shape[0] >= n
        xi = xi0[:n, :].copy()
        xo = xo0[:n, :].copy()
        print "init:", ei,

        repeats = 1 # max(int(math.ceil(10/n)), 2)
        import timeit
        tt = ()
        for test_i in range(COUNT):
            #t1 = timeit.timeit(test_scripts[0], "from __main__ import test3", number=repeats)
            #t2 = timeit.timeit(test_scripts[1], "from __main__ import test2", number=repeats)

            #test5()
            #exit()

            if USE_KERNPROF:
                test3()
                test5()
                exit()

            if USE_KERNPROF:
                qq1 = test3()
                qq2 = test4() #python -O -m kernprof -v -l   vectorised_bisection.py
                exit()


            t1 = timeit.timeit(test_scripts[test_i], "from __main__ import test"+str(test_i+1), number=repeats)
            tt += (t1/repeats,)
            #print tt
        #tl.append((t1/repeats, t2/repeats))
        tl.append(tt)
        print ".",
    ta = np.array(tl)
    #print ta.shape
    print ta * 1000  # in msec
    print "speedup ratios:"
    for test_i in range(1, COUNT):
        print lbl[test_i], "ratio = ", ta[:, 0]/ta[:, test_i]
        # print "ratio =", ta[:, 0]/ta[:, test_i]
    print "abolute time:",
    test_i = 5-1
    print lbl[test_i], "(sec) = ", ta[:, test_i]

    naa = np.array(na)


    if optimised_used():
        asymp_ratio1 = 0.00001
        asymp_ratio2 = 0.00001
        asymp_bias = 0.01/2.
    else:
        asymp_ratio1 = 0.00001*2
        asymp_ratio2 = 0.00001*2
        asymp_bias = 0.01

    #g0 = plt.loglog(naa, naa/asymp[0]*asymp[1]/100, "b:", label='asympt')
    #sty = {'point-wise':"r*-", 'numpy vectorised':"bs-" }
    for j in range(COUNT):
        #k = sty.keys()[j]
        #plt.loglog(naa, ta[:, j], k, label=sty[k])
        plt.loglog(naa, ta[:, j], sty[j], label=lbl[j])
    #g1 = plt.loglog(naa, ta[:, 1], "r*-", label='point-wise')
    #g2 = plt.loglog(naa, ta[:, 0], "bs-", label='numpy vectorised')
    g0 = plt.loglog(naa, naa * asymp_ratio1, "b:", label='asympt')
    g01 = plt.loglog(naa, naa * asymp_ratio2+ asymp_bias, "b:", label='asympt')
    plt.legend(loc='upper left', numpoints = 1)
    plt.show()

    #g1 = plt.plot(naa, ta[:, 1], "r*-", label='point-wise')

    #g2 = plt.plot(naa, ta[:, 0], "bs-", label='numpy vectorised')
    for ti in range(1, COUNT):
        g2 = plt.plot(naa, ta[:, ti], sty[ti], label=lbl[ti])
    plt.plot(naa, naa * asymp_ratio1, "b:", label='asympt')
    g0 = plt.plot(naa, naa * asymp_ratio2 + asymp_bias, "b:", label='asympt')
    plt.legend(loc='upper left', numpoints = 1)
    plt.show()

if __name__ == "__main__":
    experiment1()

#python -O -m kernprof -v -l vectorised_bisection.py
