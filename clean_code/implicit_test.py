import unittest

import numpy as np

from basic_functions import make_random_vector3, almost_equal1, make_vector3, check_vector3

from implicit_config import TOLERANCE

from primitives import UnitSphere, UnitCube1

import example_objects

from basic_functions import check_scalar_vectorized, repeat_vect3, make_random_vector3_vectorized, check_vector3_vectorized

# import simple_blend
import vector3

NUMERICAL_GRADIENT_TOLERANCE = 0.0001  # 0.00001   # 0.001
assert NUMERICAL_GRADIENT_TOLERANCE > 0.0000059
# assert NUMERICAL_GRADIENT_TOLERANCE > 0.00000001


def almost_equal3(a, b, TOLERANCE):
    assert not np.any(np.isnan(a.ravel()))
    assert not np.any(np.isinf(b.ravel()))
    assert not issubclass(a.dtype.type, np.integer)
    check_vector3(a)
    check_vector3(b)
    return np.sum(np.abs(a - b)) < TOLERANCE


def weights(k, x0, xs):
    """Calculate weights for the finite difference approximation.
    Arguments:
    k - The k-th derivative will be approximated.
    x0 - The point at which to approximate the derivative.
    xs - The grid points at which the function's value is known.
    This uses the algorithm described in:
    B. Fornberg, "Calculation of weights in finite difference formulas",
    SIAM Review 40 (1998), pp. 685-691.
    """
    # Size of the system.
    n = xs.size
    assert k < n, "need more grid points to calculate this derivative"
    # Measure points relative to x0.
    xs = xs - x0
    # Weight matrix to calculate successive finite difference
    # approximations.
    w = np.zeros((k+1, n))
    # Before starting, we want to pre-compute certain reusable
    # quantities.
    product_ratio = np.ones(n)
    for j in range(n):
        for i in range(j):
            product_ratio[j] *= xs[j] - xs[i]
    product_ratio[1:] = product_ratio[:n-1] / product_ratio[1:]
    # Each iteration of this loop will produce weights for the
    # derivatives that use one point more than the previous iteration.
    for j in range(n - k):
        # 0-th derivative approximation.
        if j == 0:
            # The approximation to the 0th derivative given only one
            # function value is trivially to just use that value.
            w[0, 0] = 1
        else:
            w[0, j] = - xs[j-1] * w[0, j-1] * product_ratio[j]
            for i in range(1, j+1):
                w[0, j-i] = xs[j] * w[0, j-i] / (xs[j] - xs[j-i])
        for m in range(1, k+1):
            # Generate weights for each derivative using the
            # previous one.
            # m is the derivative we are currently working on,
            # and l is the number of points used in this round,
            # minus one.
            l = j+m
            w[m, l] = (m*w[m-1, l-1] - xs[l-1]*w[m, l-1]) \
                        * product_ratio[l]
            for i in range(1, l+1):
                w[m, l-i] = (xs[l] * w[m, l-i] - m*w[m-1, l-i]) \
                            / (xs[l] - xs[l-i])
    return w[k, :]


def numerical_gradient_vectorized_v2(iobj, pos0, delta_t=0.01/100., order=5):
    # not useful here because we are working with a vector whose dimension are(3,1)
    """ A proper vectorized implementation. See numerical_gradient() """
    # Note: 0.1 is not enough for delta_t
    check_vector3_vectorized(pos0)
    assert issubclass(type(iobj), vector3.ImplicitFunctionVectorized)
    assert pos0.ndim == 2
    if pos0.shape[0] == 0:
        return np.zeros((0, 3))

    m = order  # sample points: -m,...,-1,0,1,2,...,+m

    sample_points = range(-m, m+1)
    n = m*2+1

    x0 = 0
    findiff_weights = weights(k=1, x0=x0, xs=np.array(sample_points) * delta_t)
    del x0

    assert n < 20
    pos0_3 = pos0[:, np.newaxis, :]
    pos = np.tile(pos0_3, (1, 3*n, 1))
    assert not issubclass(pos.dtype.type, np.integer)

    dx = make_vector3(1, 0, 0)[np.newaxis, np.newaxis, :]
    dy = make_vector3(0, 1, 0)[np.newaxis, np.newaxis, :]
    dz = make_vector3(0, 0, 1)[np.newaxis, np.newaxis, :]
    dxyz = [dx, dy, dz]

    ci = 0
    for d in range(3):
        dd = dxyz[d]
        for i in sample_points:
            pos[:, ci, :] = pos[:, ci, :] + (dd * (delta_t * float(i)))

            assert ci < 3*n
            ci += 1

    vsize = pos0.shape[0]

    v = iobj.implicitFunction(pos.reshape((vsize * 3*n), 3))    # v .shape: (3,11)
    v3 = np.reshape(v, (vsize, 3, n), order='C')  # v3 .shape: (11,)

    if True:
        Lipchitz_B = 50  # Lipschitz constant
        Lipchitz_beta = 1  # order. Keep it 1

        b_h_beta = Lipchitz_B*(np.abs(delta_t)**Lipchitz_beta)

        d0 = np.abs(np.diff(v3, axis=1+1))
        nonsmooth_ness = d0 / (np.abs(delta_t)**Lipchitz_beta)

        del d0, b_h_beta, Lipchitz_beta, Lipchitz_B

        # print nonsmooth_ness.shape
        # set_trace()
        if(np.max(np.ravel(nonsmooth_ness))) > 100*10:
            print "warning: nonsmooth ",
        del nonsmooth_ness

    if False:

        d = np.abs(np.diff(v3, n=1, axis=1+1)) / np.abs(delta_t)
        d = d - np.tile(np.mean(d, axis=1+1, keepdims=True), (1, 1, d.shape[1+1]))
        d = np.abs(d) / np.abs(delta_t)
        d = d - np.tile(np.mean(d, axis=1+1, keepdims=True), (1, 1, d.shape[1+1]))
        del d

    """ Calculating the numerical derivative using finite difference (convolution with weights) """
    # convolusion
    grad_cnv = np.dot(v3, findiff_weights)  # "sum product over the last axis of a and the second-to-last of b"

    assert not np.any(np.isnan(grad_cnv), axis=None)
    return grad_cnv


def numerical_gradient(iobj, pos0, delta_t=0.01/10.0/10.0, order=5):

    check_vector3(pos0)

    assert issubclass(type(iobj), vector3.ImplicitFunction)

    m = order  # sample points: -m,...,-1,0,1,2,...,+m

    sample_points = range(-m, m+1)
    n = m*2+1

    x0 = 0
    findiff_weights = weights(k=1, x0=x0, xs=np.array(sample_points) * delta_t)

    pos = repeat_vect3(1, pos0)
    pos3 = np.tile(pos, (3*n, 1))

    assert not issubclass(pos.dtype.type, np.integer)

    dx = repeat_vect3(1, make_vector3(1, 0, 0))
    dy = repeat_vect3(1, make_vector3(0, 1, 0))
    dz = repeat_vect3(1, make_vector3(0, 0, 1))
    dxyz = [dx, dy, dz]

    ci = 0
    for d in range(3):
        for i in sample_points:
            dd = dxyz[d]
            pos3[ci, :] = pos3[ci, :] + (dd * delta_t * float(i))
            ci += 1

    v = iobj.implicitFunction(pos3)

    v3 = np.reshape(v, (3, n), order='C')

    Lipchitz_beta = 1  # order. Keep it 1

    d0 = np.abs(np.diff(v3, axis=1))

    nonsmooth_ness = d0 / (np.abs(delta_t)**Lipchitz_beta)

    d = np.abs(np.diff(v3, n=1, axis=1)) / np.abs(delta_t)
    d = d - np.tile(np.mean(d, axis=1, keepdims=True), (1, d.shape[1]))
    d = np.abs(d) / np.abs(delta_t)
    d = d - np.tile(np.mean(d, axis=1, keepdims=True), (1, d.shape[1]))

    if(np.max(np.ravel(nonsmooth_ness))) > 100*10:
        print "warning: nonsmooth ",

    """ Calculating the numerical derivative using finite difference (convolution with weights) """
    # convolusion
    grad_cnv = np.dot(v3, findiff_weights)

    # Detecting sharp edges (non-smooth points, i.e. corners and edges and ridges)
    if np.max(np.abs(grad_cnv)) > 100:
        pass

    return grad_cnv.reshape(1, 3)


class ImplicitFunctionTests(unittest.TestCase):

    def test_simple_sphere(self):
        """ Testing for correct evaluations of the implicit function in sample points """
        """ Should return true/false if point correctly lies in volume of the sphere """
        my_sphere = UnitSphere()
        for i in range(100):
            my_point = np.array([np.random.randn(), np.random.randn(), np.random.randn()])
            flag_outside = np.linalg.norm(my_point) - 1 > 0
            flag0_on_surface = almost_equal1(1 - np.linalg.norm(my_point), 0, TOLERANCE)
            v = my_sphere.implicitFunction(my_point)
            my_func_ouput_outside = v < 0
            my_func_ouput0_on_surface = almost_equal1(v, 0, TOLERANCE)
            self.assertEqual(flag_outside, my_func_ouput_outside, 'A point is wrongly classified')
            self.assertEqual(flag0_on_surface, my_func_ouput0_on_surface, 'A point on surface is wrongly classified')

    def test_unit_cube1(self):
        """  Points inside or outside the unit cube [0,1]^3 """
        count_inside = 0
        count_outside = 0

        my_cube = UnitCube1()
        for i in range(100):
            my_point = np.array([np.random.randn(), np.random.randn(), np.random.randn()])*0.5 + 0.5

            is_inside = (
                (my_point[0] >= 0-0.5) and (my_point[0] <= 1.0-0.5) and
                (my_point[1] >= 0-0.5) and (my_point[1] <= 1.0-0.5) and
                (my_point[2] >= 0-0.5) and (my_point[2] <= 1.0-0.5))

            if is_inside:
                count_inside += 1
            else:
                count_outside += 1

            v2 = my_cube.implicitFunction((repeat_vect3(1, my_point)).reshape(3))

            self.assertEqual(is_inside, v2 >= 0, 'A point on surface is wrongly classified')

        print "point count: inside, outside = ", count_inside, ",", count_outside

    def test_correctly_classifies_random_points(self):
        """ Testing for correct classification of random sample points """
        for i in range(50):
            v = make_random_vector3(1, 1)
            my_sphere = UnitSphere()
            self.assertTrue(almost_equal1(my_sphere.implicitFunction(v), 0.0, TOLERANCE))

    def test_gradient_function(self):
        """ Testing the gradient function of the implicit function"""
        pass

    def test_hessian_matrix(self):
        pass


class ImplicitObjectTests(unittest.TestCase):
    pass


def vectors_parallel_and_direction(v1, v2):

    assert v1.shape == (3,), str(v1.shape)
    assert v2.shape == (3,), str(v2.shape)
    cross_prod = np.cross(v1, v2)
    inner_prod = np.dot(np.transpose(v1), v2)
    are_parallel1 = almost_equal3(cross_prod, make_vector3(0, 0, 0), TOLERANCE)
    are_directed = (inner_prod > 0)

    are_parallel2 = np.allclose(cross_prod, np.array([0, 0, 0]), atol=TOLERANCE)
    return (are_parallel1 and are_parallel2, are_directed)


class ImplicitFunctionVectorizedTests(unittest.TestCase):

    def ellipsoid_point_and_gradient_vectorized(self, m, xa, correctGrad, center=None):
        """ Checks if the point x is on the surface, and if the gradient is correct."""

        e = vector3.Ellipsoid(m)
        msg = "Ellipsoid(m): "+str(e)

        va = e.implicitFunction(xa)
        ga = e.implicitGradient(xa)
        check_vector3_vectorized(ga)
        N = xa.shape[0]
        check_scalar_vectorized(va, N)
        assert ga.ndim == 2
        assert ga.shape == (N, 3)
        assert correctGrad.shape == (N, 3)

        correctScalar = 0
        less_a = np.less(np.abs(va - correctScalar), TOLERANCE)

        if not np.all(less_a):
            print("Some error:")
            print(xa)
            print(va)
            print(ga)
            print(e)

        self.assertTrue(np.all(less_a), ("Implicit Function's scalar value incorrect"))

        for i in range(ga.shape[0]):
            (are_parallel, are_directed) = vectors_parallel_and_direction(ga[i, :], correctGrad[i, :])
            self.assertTrue(are_parallel, "Incorrect gradient: not parallel "+msg)
            self.assertTrue(are_directed, "parallel but opposite directions "+msg)

    def test_ellipsoid_certain_points(self):
        """ Ellipsoid using vectorized calculations """
        N = 10
        _x = make_vector3(0, 0, 1)
        xa = repeat_vect3(N, _x)

        m = np.eye(4)  # makeMatrix4( 1,0,0,0, 0,1,0,0,  0,0,1,0 )
        m[0, 0] = 1
        ga2 = repeat_vect3(N, make_vector3(0, 0, -1))
        self.ellipsoid_point_and_gradient_vectorized(m, xa, ga2)

        m = np.eye(4)
        self.ellipsoid_point_and_gradient_vectorized(m, repeat_vect3(N, make_vector3(0, 1, 0)), repeat_vect3(N, make_vector3(0, -1, 0)))

        m = np.eye(4)
        self.ellipsoid_point_and_gradient_vectorized(m, repeat_vect3(N, make_vector3(1, 0, 0)), repeat_vect3(N, make_vector3(-1, 0, 0)))

        xa = repeat_vect3(N, make_vector3(0, 0, 2))
        m = np.eye(4)
        m[2, 2] = 2
        self.ellipsoid_point_and_gradient_vectorized(m, xa, repeat_vect3(N, make_vector3(0, 0, -1)))

        xa = repeat_vect3(N, make_vector3(2, 0, 0))
        m = np.eye(4)
        m[0, 0] = 2
        self.ellipsoid_point_and_gradient_vectorized(m, xa, repeat_vect3(N, make_vector3(-1, 0, 0)))

        xa = repeat_vect3(N, make_vector3(0, 2, 0))
        m = np.eye(4)
        m[0, 0] = 2
        m[1, 1] = 2
        m[2, 2] = 2
        self.ellipsoid_point_and_gradient_vectorized(m, xa, repeat_vect3(N, make_vector3(0, -2, 0)))

    def test_ellipsoid_random_points(self):
        """Testing hundreds of random points on a sphere of size RADIUS"""
        for i in range(0, 30):

            RADIUS = 3
            POW = 4  # higher POW will get points more toward parallel to axes
            N = 500

            rcenter = make_random_vector3(1000, 1.0)
            centers_a = repeat_vect3(N, rcenter)
            r0 = make_random_vector3_vectorized(N, RADIUS, POW)
            r = r0 + centers_a

            assert r.shape[1] == 3
            xa = r

            m = np.eye(4) * RADIUS
            m[0:3, 3] = rcenter
            m[3, 3] = 1

            expected_grad = -r0
            check_vector3_vectorized(expected_grad)

            self.ellipsoid_point_and_gradient_vectorized(m, xa, expected_grad)

    def test_gradients_using_numerical_gradient(self):
        def blend2():
            m1 = np.eye(4) * 1
            m1[0:3, 3] = [0, 1, 0]
            m1[3, 3] = 1

            m2 = np.eye(4) * 2
            m2[0:3, 3] = [2.5, 0, 0]
            m2[3, 3] = 1

            iobj_v = vector3.SimpleBlend(vector3.Ellipsoid(m1), vector3.Ellipsoid(m2))
            return iobj_v

        iobj_v = blend2()
        self.check_gradient_function(iobj_v, objname="blend2")

        examples_list = example_objects.get_all_examples([2])
        for example_name in examples_list:
            print("example_name = ", example_name)
            iobj = example_objects.make_example_vectorized(example_name)
            self.check_gradient_function(iobj, objname=example_name)

        """ numerical """
        x = make_vector3(0, 2, 0)
        g2 = numerical_gradient(iobj_v, x)
        g = iobj_v.implicitGradient(repeat_vect3(1, x))
        np.set_printoptions(formatter={'all': lambda x: ''+("%2.19f" % (x,))})
        err = np.sum(np.abs(g - g2), axis=1)
        err_max = np.max(np.abs(g - g2), axis=1)
        err_rel = np.sum(np.abs(g - g2), axis=1) / np.mean(np.abs(g))
        print(err, err_max, err_rel)
        # print(err)
        self.assertTrue(np.all(err < NUMERICAL_GRADIENT_TOLERANCE))

    def check_gradient_function(self, iobj, tolerance=NUMERICAL_GRADIENT_TOLERANCE, objname=None):
        """testing the gradient on 100 random points """

        assert (-1)*-tolerance > 0, "should be a number"

        if iobj is not None:
            from vector3 import ImplicitFunction
            assert issubclass(type(iobj), ImplicitFunction)

        for i in range(100):
            radius = 10  # 10/5
            x = make_random_vector3(radius, 1, type="randn")  # [0:3]
            x = x * np.random.randn() * 5
            self.check_gradient_function_point1(iobj, x, tolerance, objname)

        x = make_vector3(0, 2, 0)
        self.check_gradient_function_point1(iobj, x, tolerance, objname)

    def check_two_vectors(self, a, b, tolerance, description):
        err = np.sum(np.abs(a - b), axis=1)
        err_max = np.max(np.abs(a - b), axis=1)
        err_rel = np.sum(np.abs(a - b), axis=1) / np.mean(np.abs(a))
        if err >= tolerance:
            print("err: ", err)
            print("err(max): ", err_max, "  err_rel:", err_rel)
            print("a,g_numerical: ", a, b)

        err_n = np.sum(np.abs(a - b), axis=1)
        if np.any(err_n >= tolerance):
            print(np.any(err_n < tolerance), "aaa")
            # print(err_n.shape, "*******")
            print(a, b, a - b)

    #    self.assertTrue(np.all(err_n < tolerance), msg="absolute error exceeds: %f %s" % (np.max(err_n), description,))
    #    self.assertTrue(np.all(err < tolerance), msg=" %s" % (description,))

    def check_gradient_function_point1(self, iobj, x, tolerance, objname):
        """Tests the gradient using numerical method, verify with analytical and vectorized-analytical"""

        assert iobj is not None

        if iobj is not None:
            g_numer_vec = numerical_gradient(iobj, x)

            from vector3 import ImplicitFunction
            assert issubclass(type(iobj), ImplicitFunction)
            g_vec = iobj.implicitGradient(repeat_vect3(1, x))

            np.set_printoptions(formatter={'all': lambda x: ''+("%2.7f" % (x,))})
            self.check_two_vectors(g_vec, g_numer_vec, tolerance, "vec numerical versus analytical gradients (%s)" % (objname,))


if __name__ == '__main__':
    unittest.main()
