import unittest

import numpy as np

from basic_functions import make_random_vector3, almost_equal1, make_vector3, check_vector3

from implicit_config import TOLERANCE

from basic_functions import check_scalar_vectorized, repeat_vect3, make_random_vector3_vectorized, check_vector3_vectorized

# import simple_blend
import vector3

import numerical_utils

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


class ImplicitFunctionTests(unittest.TestCase):

    def test_simple_sphere(self):
        """ Testing for correct evaluations of the implicit function in sample points """
        """ Should return true/false if point correctly lies in volume of the sphere """
        my_sphere = vector3.UnitSphere()
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

        my_cube = vector3.UnitCube1()
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

            v2 = my_cube.implicitFunction(repeat_vect3(1, my_point))

            self.assertEqual(is_inside, v2 >= 0, 'A point on surface is wrongly classified')

        print "point count: inside, outside = ", count_inside, ",", count_outside

    def test_correctly_classifies_random_points(self):
        """ Testing for correct classification of random sample points """
        for i in range(50):
            v = make_random_vector3(1, 1)
            my_sphere = vector3.UnitSphere()
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

    are_parallel1 = almost_equal3(make_vector3(cross_prod), make_vector3(0, 0, 0), TOLERANCE)
    are_directed = (inner_prod > 0)

    are_parallel2 = np.allclose(cross_prod, np.array([0, 0, 0]), atol=TOLERANCE)
    return (are_parallel1 and are_parallel2, are_directed)


class EllipsoidTests(unittest.TestCase):

    def ellipsoid_point_and_gradient(self, m, x, correctGrad, center=None):
        """ Checks if the point x is on the surface, and if the gradient is correct."""
        e = vector3.Ellipsoid(m)
        v = e.implicitFunction(x)
        g = e.implicitGradient(x)
        check_vector3(g)
        correctScalar = 0
        assert np.abs(v - correctScalar) < TOLERANCE, ("Implicit Function's scalar value incorrect: %2.20f" % (v,))

        msg = "vector3.Ellipsoid(m): "+str(e)

        (are_parallel, are_directed) = vectors_parallel_and_direction(g, correctGrad)
        self.assertTrue(are_parallel, "Incorrect gradient: not parallel "+msg)
        self.assertTrue(are_directed, "parallel but opposite directions "+msg)

    def test_ellipsoid_certain_points(self):
        """Ellipsoid & unit sphere hard coded values"""
        x = make_vector3(0, 0, 1)
        m = np.eye(4)  # makeMatrix4( 1,0,0,0, 0,1,0,0,  0,0,1,0 )
        m[0, 0] = 1
        self.ellipsoid_point_and_gradient(m, x, make_vector3(0, 0, -1))

        m = np.eye(4)
        self.ellipsoid_point_and_gradient(m, make_vector3(0, 1, 0), make_vector3(0, -1, 0))

        m = np.eye(4)
        self.ellipsoid_point_and_gradient(m, make_vector3(1, 0, 0), make_vector3(-1, 0, 0))

        x = make_vector3(0, 0, 2)
        m = np.eye(4)
        m[2, 2] = 2
        self.ellipsoid_point_and_gradient(m, x, make_vector3(0, 0, -1))

        x = make_vector3(2, 0, 0)
        m = np.eye(4)
        m[0, 0] = 2
        self.ellipsoid_point_and_gradient(m, x, make_vector3(-1, 0, 0))

        x = make_vector3(0, 2, 0)
        m = np.eye(4)
        m[0, 0] = 2
        m[1, 1] = 2
        m[2, 2] = 2
        self.ellipsoid_point_and_gradient(m, x, make_vector3(0, -2, 0))

    def test_ellipsoid_random_points(self):
        """Testing hundreds of random points on a sphere of size RADIUS=3"""
        for i in range(0, 100):

            RADIUS = 3
            POW = 4  # higher POW will get points more toward parallel to axes

            rcenter = make_random_vector3(1000, 1.0)
            r0 = make_random_vector3(RADIUS, POW)
            r = r0 + rcenter
            assert r.shape[0] == 3
            x = make_vector3(r[0], r[1], r[2])

            m = np.eye(4) * RADIUS
            m[0:3, 3] = rcenter[0:3]
            m[3, 3] = 1

            expected_grad = make_vector3(-r0[0], -r0[1], -r0[2])
            check_vector3(expected_grad)

            self.ellipsoid_point_and_gradient(m, x, expected_grad)


class ImplicitFunctionVectorizedTests(unittest.TestCase):

    def ellipsoid_point_and_gradient_vectorized(self, m, xa, correctGrad, center=None):
        """ Checks if the point x is on the surface, and if the gradient is correct."""

        e = vector3.Ellipsoid(m)
        msg = "vector3.Ellipsoid(m): "+str(e)

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
        g2 = numerical_utils.numerical_gradient(iobj_v, x, is_vectorized=True)
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

        self.assertTrue(np.all(err_n < tolerance), msg="absolute error exceeds: %f %s" % (np.max(err_n), description,))
        self.assertTrue(np.all(err < tolerance), msg=" %s" % (description,))

    def check_gradient_function_point1(self, iobj, x, tolerance, objname):
        """Tests the gradient using numerical method, verify with analytical and vectorized-analytical"""

        assert iobj is not None

        if iobj is not None:
            g_numer_vec = numerical_utils.numerical_gradient(iobj, x, is_vectorized=True)

            from vector3 import ImplicitFunction
            assert issubclass(type(iobj), ImplicitFunction)
            g_vec = iobj.implicitGradient(repeat_vect3(1, x))

            np.set_printoptions(formatter={'all': lambda x: ''+("%2.7f" % (x,))})

            self.check_two_vectors(g_vec, g_numer_vec, tolerance, "vec numerical versus analytical gradients (%s)" % (objname,))

import example_objects


class Examples(unittest.TestCase):

    def test_examples(self):
        """ dummy test. not necessary"""
        # iobj = example_objects.bowl_hole()
        # iobj = example_objects.cube1()
        # iobj = example_objects.csg_example1()
        # iobj = example_objects.dice()
        # iobj = example_objects.rdice()

        usable_examples = example_objects.get_all_examples([2])

        # ["sphere_example", "ell_example1", "blend_example2", "cube_example", "blend_example2_discs", "blend_example1",
        # "bowl_15_holes", "first_csg", "french_fries", "rdice_vec", "rcube_vec", "screw1", "screw2", "udice_vec", "rods",
        # "cyl1", "cyl2", "cyl3", "cyl4", "cube_with_cylinders", "union_of_two_cubes", 'crisp_cube_sphere', 'cube_modified']
        print(usable_examples)

        iobj = example_objects.make_example_vectorized('csg_example1')

        """ Implicit object is not defined. Now going for visualisation. """

        x = make_vector3(0.5, 0.5, 0.5)
        #   x = make_vector4( 1, 1, 1 )
        g = iobj.implicitGradient(x)
        v = iobj.implicitFunction(x)



if __name__ == '__main__':
    unittest.main()
