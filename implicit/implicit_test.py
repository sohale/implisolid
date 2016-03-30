import unittest

import numpy as np

from basic_types import check_vector4, make_vector4, make_vector4_numpy, make_random_vector, almost_equal1, almost_equal4

from implicit_config import TOLERANCE

from primitives import UnitSphere
from primitives import UnitCube1

#from ellipsoid import Ellipsoid
import nonvec

#vectorizes tests
#import implicit_vectorized
#import simple_blend_nonvec
import ellipsoid_vectorized

from basic_types import check_vector4_vectorized, check_scalar_vectorized, make_random_vector_vectorized, repeat_vect4

#import simple_blend
import vectorized

import numerical_utils

#import example_objects

NUMERICAL_GRADIENT_TOLERANCE = 0.0001 # 0.00001   # 0.001
assert NUMERICAL_GRADIENT_TOLERANCE > 0.0000059
assert NUMERICAL_GRADIENT_TOLERANCE > 0.00000001


class ImplicitFunctionTests(unittest.TestCase):

    # def setUp(self):
    #     print ( "Setting up ImplicitTests: _start" )
    #
    #     print ( "Setting up ImplicitTests: _end" )
    #
    # def tearDown(self):
    #     print ( "Tearing down ImplicitTests: _start")
    #     # set up objects or init
    #     print ("Tearing down ImplicitTests: _end")

    def test_simple_sphere(self):
        """ Testing for correct evaluations of the implicit function in sample points """
        """ Should return true/false if point correctly lies in volume of the sphere """
        my_sphere = UnitSphere()
        for i in range(100):
            my_point = np.array([np.random.randn(), np.random.randn(), np.random.randn(), 1])
            flag_outside = np.linalg.norm(my_point[0:3]) - 1 > 0
            flag0_on_surface = almost_equal1(1 - np.linalg.norm(my_point[0:3]), 0, TOLERANCE)
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
        my_cube_vec = vectorized.UnitCube1()
        for i in range(100):
            my_point = np.array([np.random.randn(), np.random.randn(), np.random.randn(), 1])*0.5 + 0.5
            my_point[3] = 1
            #print(my_point)
            is_inside = (
                (my_point[0] >= 0-0.5) and (my_point[0] <= 1.0-0.5) and
                (my_point[1] >= 0-0.5) and (my_point[1] <= 1.0-0.5) and
                (my_point[2] >= 0-0.5) and (my_point[2] <= 1.0-0.5))
            #if is_inside:
            #    print (my_point[0:3])
            if is_inside:
                count_inside += 1
            else:
                count_outside += 1

            v = my_cube.implicitFunction(my_point)
            #print( is_inside, v )

            self.assertEqual(is_inside, v >= 0, 'A point on surface is wrongly classified')

            v2 = my_cube_vec.implicitFunction( repeat_vect4(1,my_point) )
            #print(v2, my_point)

            self.assertEqual(is_inside, v2 >= 0, 'A point on surface is wrongly classified')


        print("point count: inside, outside = ", count_inside, ",", count_outside)

    def test_correctly_classifies_random_points(self):
        """ Testing for correct classification of random sample points """
        for i in range(50):
            v = make_random_vector(1, 1)
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
    #check_vector3(v1)
    #check_vector3(v2)
    #print(v1.shape)
    assert v1.shape == (3,) , str(v1.shape)
    assert v2.shape == (3,) , str(v2.shape)
    cross_prod = np.cross(v1[0:3], v2[0:3])
    inner_prod = np.dot(np.transpose(v1[0:3]), v2[0:3])

    are_parallel1 = almost_equal4(make_vector4_numpy(cross_prod), make_vector4(0, 0, 0), TOLERANCE)
    are_directed = (inner_prod > 0)  #no tolerance is needed

    are_parallel2 = np.allclose(cross_prod, np.array([0, 0, 0]), atol=TOLERANCE)
    return (are_parallel1 and are_parallel2, are_directed)

class EllipsoidTests(unittest.TestCase):

    def ellipsoid_point_and_gradient(self, m, x, correctGrad, center=None):
        """ Checks if the point x is on the surface, and if the gradient is correct."""
        e = nonvec.Ellipsoid(m)
        v = e.implicitFunction(x)
        g = e.implicitGradient(x)
        check_vector4(g)
        correctScalar = 0
        assert np.abs(v - correctScalar) < TOLERANCE, ("Implicit Function's scalar value incorrect: %2.20f" % (v,))
        #assert np.allclose( g, correctGrad , atol=TOLERANCE ) , "Incorrect gradient"
        msg = "nonvec.Ellipsoid(m): "+str(e)

        (are_parallel,are_directed) = vectors_parallel_and_direction(g[0:3], correctGrad[0:3])
        self.assertTrue(are_parallel, "Incorrect gradient: not parallel "+msg)
        self.assertTrue( are_directed, "parallel but opposite directions "+msg)

    def test_ellipsoid_certain_points(self):
        """Ellipsoid & unit sphere hard coded values"""
        x = make_vector4(0, 0, 1)
        m = np.eye(4)  # makeMatrix4( 1,0,0,0, 0,1,0,0,  0,0,1,0 )
        m[0, 0] = 1
        self.ellipsoid_point_and_gradient(m, x, make_vector4(0, 0, -1))

        m = np.eye(4)
        self.ellipsoid_point_and_gradient(m, make_vector4(0, 1, 0), make_vector4(0, -1, 0))

        m = np.eye(4)
        self.ellipsoid_point_and_gradient(m, make_vector4(1, 0, 0), make_vector4(-1, 0, 0))

        x = make_vector4(0, 0, 2)
        m = np.eye(4)
        m[2, 2] = 2
        self.ellipsoid_point_and_gradient(m, x, make_vector4(0, 0, -1))

        x = make_vector4(2, 0, 0)
        m = np.eye(4)
        m[0, 0] = 2
        self.ellipsoid_point_and_gradient(m, x, make_vector4(-1, 0, 0))

        x = make_vector4(0, 2, 0)
        m = np.eye(4)
        m[0, 0] = 2
        m[1, 1] = 2
        m[2, 2] = 2
        self.ellipsoid_point_and_gradient(m, x, make_vector4(0, -2, 0))

    def test_ellipsoid_random_points(self):
        """Testing hundreds of random points on a sphere of size RADIUS=3"""
        for i in range(0, 100):

            RADIUS = 3
            POW = 4  # higher POW will get points more toward parallel to axes

            rcenter = make_random_vector(1000, 1.0)[0:3]
            r0 = make_random_vector(RADIUS, POW)[0:3]
            r = r0 + rcenter
            assert r.shape[0] == 3
            x = make_vector4(r[0], r[1], r[2])

            m = np.eye(4) * RADIUS
            m[0:3, 3] = rcenter[0:3]
            m[3, 3] = 1

            expected_grad = make_vector4(-r0[0], -r0[1], -r0[2])
            check_vector4(expected_grad)

            self.ellipsoid_point_and_gradient(m, x,  expected_grad)


class ImplicitFunctionVectorizedTests(unittest.TestCase):

    def ellipsoid_point_and_gradient_vectorized(self, m, xa, correctGrad, center=None):
        """ Checks if the point x is on the surface, and if the gradient is correct."""

        e = vectorized.Ellipsoid(m)
        msg = "vectorized.Ellipsoid(m): "+str(e)

        va = e.implicitFunction(xa)
        ga = e.implicitGradient(xa)
        check_vector4_vectorized(ga)
        N = xa.shape[0]
        check_scalar_vectorized(va, N)
        assert ga.ndim == 2
        assert ga.shape == (N, 4)
        assert correctGrad.shape == (N, 4)

        correctScalar = 0
        less_a = np.less(np.abs(va - correctScalar), TOLERANCE)
        #print(va - correctScalar)

        if not np.all(less_a):
            print("Some error:")
            print(xa)
            print(va)
            print(ga)
            print(e)

        self.assertTrue(np.all(less_a), ("Implicit Function's scalar value incorrect"))

        for i in range(ga.shape[0]):
            (are_parallel,are_directed) = vectors_parallel_and_direction(ga[i, 0:3], correctGrad[i, 0:3])
            self.assertTrue(are_parallel, "Incorrect gradient: not parallel "+msg)
            self.assertTrue( are_directed, "parallel but opposite directions "+msg)

    def test_ellipsoid_certain_points(self):
        """ Ellipsoid using vectorized calculations """
        N = 10
        _x = make_vector4(0, 0, 1)
        xa = repeat_vect4(N, _x)

        m = np.eye(4)  # makeMatrix4( 1,0,0,0, 0,1,0,0,  0,0,1,0 )
        m[0, 0] = 1
        ga2 = repeat_vect4(N, make_vector4(0, 0, -1))
        self.ellipsoid_point_and_gradient_vectorized(m, xa, ga2)

        m = np.eye(4)
        self.ellipsoid_point_and_gradient_vectorized(m, repeat_vect4(N, make_vector4(0, 1, 0)), repeat_vect4(N, make_vector4(0, -1, 0)))

        m = np.eye(4)
        self.ellipsoid_point_and_gradient_vectorized(m, repeat_vect4(N, make_vector4(1, 0, 0)), repeat_vect4(N, make_vector4(-1, 0, 0)))

        xa = repeat_vect4(N, make_vector4(0, 0, 2))
        m = np.eye(4)
        m[2, 2] = 2
        self.ellipsoid_point_and_gradient_vectorized(m, xa, repeat_vect4(N, make_vector4(0, 0, -1)))

        xa = repeat_vect4(N, make_vector4(2, 0, 0))
        m = np.eye(4)
        m[0, 0] = 2
        self.ellipsoid_point_and_gradient_vectorized(m, xa, repeat_vect4(N, make_vector4(-1, 0, 0)))

        xa = repeat_vect4(N, make_vector4(0, 2, 0))
        m = np.eye(4)
        m[0, 0] = 2
        m[1, 1] = 2
        m[2, 2] = 2
        self.ellipsoid_point_and_gradient_vectorized(m, xa, repeat_vect4(N, make_vector4(0, -2, 0)))

    def test_ellipsoid_random_points(self):
        """Testing hundreds of random points on a sphere of size RADIUS"""
        for i in range(0, 30):

            RADIUS = 3
            POW = 4  # higher POW will get points more toward parallel to axes
            N = 500

            #rcenter = make_random_vector_vectorized(N, 1000, 1.0)
            rcenter = make_random_vector(1000, 1.0)
            centers_a = repeat_vect4(N, rcenter)
            r0 = make_random_vector_vectorized(N, RADIUS, POW)
            r = r0 + centers_a
            r[:, 3] = 1
            assert r.shape[1] == 4  # generates vector4 , but the non vectorized version generates vector 3
            xa = r  # make_vector4( r[0], r[1], r[2] )

            m = np.eye(4) * RADIUS
            m[0:3, 3] = rcenter[0:3]
            m[3, 3] = 1

            expected_grad = -r0  # r - centers_a #repeat_vect4(N, make_vector4( r0[0], r0[1], r0[2] ) )
            expected_grad[:,3] = 1
            check_vector4_vectorized(expected_grad)

            self.ellipsoid_point_and_gradient_vectorized(m, xa,  expected_grad)

    def test_gradients_using_numerical_gradient(self):
        def blend2():
            m1 = np.eye(4) * 1
            m1[0:3, 3] = [0, 1, 0]
            m1[3, 3] = 1

            m2 = np.eye(4) * 2
            m2[0:3, 3] = [2.5, 0, 0]
            m2[3, 3] = 1

            iobj_v = vectorized.SimpleBlend(vectorized.Ellipsoid(m1), vectorized.Ellipsoid(m2))
            return iobj_v

        iobj_v = blend2()
        self.check_gradient_function(iobj_v, objname="blend2")

        #examples_list = ["sphere_example", "ell_example1", "blend_example2"]
        examples_list = example_objects.get_all_examples([3])
        for example_name in examples_list:
            (iobj_v, io_) = example_objects.make_example_pair(example_name)
            self.check_gradient_function(iobj_v,  iobj_nonvec=io_, objname=example_name)


        #examples_list = example_objects.get_all_examples([3])
        #print(examples_list)

        examples_list = example_objects.get_all_examples([3])
        #examples_list = ["sphere_example", "ell_example1", "blend_example2"]
        for example_name in examples_list:
            (iobj_v, io_) = example_objects.make_example_pair(example_name)
            self.check_equal_vec_nonvec(iobj_v, io_, objname=example_name)

        examples_list = example_objects.get_all_examples([1])
        for example_name in examples_list:
            io_ = example_objects.make_example_nonvec(example_name)
            self.check_gradient_function(None,  iobj_nonvec=io_, objname=example_name)

        examples_list = example_objects.get_all_examples([2])
        for example_name in examples_list:
            print("example_name = ",example_name)
            iobj_vev = example_objects.make_example_vectorized(example_name)
            self.check_gradient_function(iobj_vev, iobj_nonvec=None, objname=example_name)


        """ numerical """
        x = make_vector4(0, 2, 0)
        g2 = numerical_utils.numerical_gradient(iobj_v, x, is_vectorized=True)
        g = iobj_v.implicitGradient(repeat_vect4(1, x))
        np.set_printoptions(formatter={'all': lambda x: ''+("%2.19f" % (x,))})
        err = np.sum(np.abs(g - g2), axis=1)
        err_max = np.max(np.abs(g - g2), axis=1)
        err_rel = np.sum(np.abs(g - g2), axis=1) / np.mean(np.abs(g))
        print(err, err_max, err_rel)
        #print(err)
        self.assertTrue(np.all(err < NUMERICAL_GRADIENT_TOLERANCE))

    def check_equal_vec_nonvec(self, vec_obj, nonvec_obj, objname):
        for i in range(100):
            radius = 10
            x = make_random_vector(radius, 1, type="randn")  # [0:3]
            x_v = repeat_vect4(1,x)
            v1 = vec_obj.implicitFunction(x_v)
            v2 = nonvec_obj.implicitFunction(x)
            g1 = vec_obj.implicitGradient(x_v)
            g2 = nonvec_obj.implicitGradient(x)
            #h1 = vec_obj.hessianMatrix(x_v)
            #h2 = nonvec_obj.hessianMatrix(x)

            err = np.sum(np.abs(g1 - g2[np.newaxis,:]), axis=1)
            self.assertTrue(np.all(err < 0.000001))

            self.check_two_vectors(g1, g2, NUMERICAL_GRADIENT_TOLERANCE, "gradient(vect versus nonvect) %s"%(objname,))


    def check_gradient_function(self, iobj_vec, tolerance=NUMERICAL_GRADIENT_TOLERANCE, iobj_nonvec=None, objname=None):
        """testing the gradient on 100 random points """

        assert (-1)*-tolerance > 0, "should be a number"

        if not iobj_vec is None:
            from vectorized import ImplicitFunctionVectorized
            assert issubclass(type(iobj_vec), ImplicitFunctionVectorized)
            import vectorized
            assert vectorized.is_implicit_type(iobj_vec)
        if not iobj_nonvec is None:
            from nonvec import ImplicitFunctionPointwise
            assert issubclass(type(iobj_nonvec), ImplicitFunctionPointwise)
            import nonvec
            assert nonvec.is_implicit_type(iobj_nonvec)

        for i in range(100):
            radius = 10  # 10/5
            x = make_random_vector(radius, 1, type="randn")  # [0:3]
            x = x * np.random.randn() * 5
            x[3]=1
            self.check_gradient_function_point1(iobj_vec, x, tolerance, iobj_nonvec, objname)

        x = make_vector4(0, 2, 0)
        self.check_gradient_function_point1(iobj_vec, x, tolerance, iobj_nonvec, objname)

    def check_two_vectors(self, a, b, tolerance, description):
        err = np.sum(np.abs(a - b), axis=1)
        err_max = np.max(np.abs(a - b), axis=1)
        err_rel = np.sum(np.abs(a - b), axis=1) / np.mean(np.abs(a))
        if err >= tolerance:
            #print("x: ", x)
            print("err: ", err)
            print("err(max): ", err_max, "  err_rel:", err_rel)
            print("a,g_numerical: ", a, b)

            #print(">")
            #if not iobj_nonvec is None:
            #    g_nonvec = iobj_nonvec.implicitGradient(x)
            #    print("nonvec: ", g_nonvec, "vec: ", a, " numeric: ", b)
            #    print("difference: ", a - g_nonvec)
            #print("<")

        err_n = np.sum(np.abs(a - b), axis=1)
        if np.any(err_n >= tolerance):
            print(np.any(err_n < tolerance), "aaa")
            #print(err_n.shape, "*******")
            print(a,b, a - b)
        #assert err_n < tolerance
        self.assertTrue(np.all(err_n < tolerance), msg="absolute error exceeds: %f %s"%(np.max(err_n),description,))
        self.assertTrue(np.all(err < tolerance), msg=" %s"%(description,))


    def check_gradient_function_point1(self, iobj_vec, x, tolerance, iobj_nonvec, objname):
        """Tests the gradient using numerical method, verify with analytical and vectorized-analytical"""

        assert not iobj_vec is None or not iobj_nonvec is None

        if not iobj_vec is None:
            #try:
            g_numer_vec = numerical_utils.numerical_gradient(iobj_vec, x, is_vectorized=True)
            #if g_numer_vec == "sharp-edge":
            #    return
        if not iobj_nonvec is None:
            g_numer_nonvec = numerical_utils.numerical_gradient(iobj_nonvec, x, is_vectorized=False)

        if not iobj_vec is None:
            from vectorized import ImplicitFunctionVectorized
            assert issubclass(type(iobj_vec), ImplicitFunctionVectorized)
            g_vec = iobj_vec.implicitGradient(repeat_vect4(1, x))

        if not iobj_nonvec is None:
            from nonvec import ImplicitFunctionPointwise
            assert issubclass(type(iobj_nonvec), ImplicitFunctionPointwise)

            g_nonvec = iobj_nonvec.implicitGradient(  x )
            self.check_two_vectors(g_nonvec, g_numer_nonvec, tolerance, "nonvec numerical versus analytical gradients (%s)"%(objname,))
            #err_n = np.sum(np.abs(g_nonvec - g_numer_vec), axis=1)
            #if(err_n >= tolerance):
            #    print(g_nonvec,g_numer_vec, g_nonvec - g_numer_vec)
            #assert err_n < tolerance

        if not iobj_vec is None and not iobj_nonvec is None:
            #err_n2 = np.sum(np.abs(g_nonvec - g_vec), axis=1)
            #assert err_n < tolerance
            self.check_two_vectors(g_nonvec, g_vec, tolerance, "nonvec versus vec gradients (%s)"%(objname,))

        np.set_printoptions(formatter={'all': lambda x: ''+("%2.7f" % (x,))})
        #print(g_vec, " =?= ", g_numer_vec)

        if not iobj_vec is None:
            self.check_two_vectors(g_vec, g_numer_vec, tolerance, "vec numerical versus analytical gradients (%s)"%(objname,))

        if not iobj_vec is None and not iobj_nonvec is None:
            self.check_two_vectors( g_numer_nonvec,  g_numer_vec, tolerance, "nonvec versus vec numerical gradients (%s)"%(objname,))

        #self.assertTrue(np.all(err < tolerance))
        #except:
        #    print("point skipped because on a sharp edge")



import example_objects



class Examples(unittest.TestCase):
    def test_creation_of_all_Examples(self):
        example_objects.test_creation_of_all_Examples()

    def test_examples(self):
        """ dummy test. not necessary"""
        #iobj = example_objects.bowl_hole()
        #iobj = example_objects.cube1()
        #iobj = example_objects.csg_example1()
        #iobj = example_objects.dice()
        #iobj = example_objects.rdice()
        #iobj = example_objects.make_example_nonvec("rdice")

        usable_examples = example_objects.get_all_examples([1])

        #['rods', 'bowl_hole', 'dice', 'csg_example1', 'rdice', 'cube1']
        print(usable_examples)
        #iobj = example_objects.make_example_nonvec(usable_examples[1])
        #iobj = example_objects.make_example_nonvec('rods')
        iobj = example_objects.make_example_nonvec('csg_example1')

        """ Implicit object is not defined. Now going for visualisation. """

        x = make_vector4(0.5, 0.5, 0.5)
        #   x = make_vector4( 1, 1, 1 )
        g = iobj.implicitGradient(x)
        v = iobj.implicitFunction(x)

        #print(v)
        #print(g)

if __name__ == '__main__':
    unittest.main()
