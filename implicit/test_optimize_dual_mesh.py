import unittest
import numpy as np
import timeit
from optimize_dual_mesh import MeshOptimizer
import random

class BisectionRootTests(unittest.TestCase):
    """ Tests for the bisection method of the Mesh Optimizer class """
    def setUp(self):
        pass
    def test_find_bisection_root_vector(self):

        """
        It should test that the bisection function provides
        a fast and approximate solution for a scalar unction with
        vector input

        Function : Unit-Sphere 1 - (x^2 + y^2 + z^2) = 0
        PointA   : [0,0,0]
        PointA   : [2,0,0]

        """
        known_function = lambda x: 1 - np.dot(x[:-1],x[:-1])
        new_opt = MeshOptimizer()
        new_opt.register_function(known_function)
        result = new_opt.find_bisection_root(np.array([0,0,0,1]),np.array([2,0,0,1]))
        result_expected = np.array([1,0,0,1],dtype=np.float32)
        print result,result_expected
        self.assertEqual(np.linalg.norm(result_expected - result),0,"It should return the point [1,0,0]")

        def tearDown(self):
            pass

class MeshOptimizerTests(unittest.TestCase):
    """ Tests for the Mesh Optimizer class  """
    def setUp(self):
        pass

    def test_should_return_a_vector(self):
        """ The vector returned by optimize_centroid should be a 4 x 1 vector """
        # new_opt = MeshOptimizer()
        # new_opt.load_example('bowl_15_holes')
        # centroid = random.choice(new_opt.centroids)
        # result = new_opt.optimize_centroid(centroid)
        # self.assertEqual(result.shape, (4,))
        pass

    def test_implicit_on_object(self):
        """
        It should test that the bisection function provides
        a fast and approximate solution for a scalar function with
        scalar input
        """
        new_opt = MeshOptimizer()
        new_opt.load_example('bowl_15_holes')
        func_evals = np.ones(len(new_opt.vertices))
        vertices = np.c_[new_opt.vertices,func_evals]
        for i in range(len(new_opt.vertices)):
            func_evals[i] =  new_opt.function(vertices[i][:].reshape(1,4))
            print func_evals[i]
        self.assertEqual(np.linalg.norm(func_evals),0,"All the initial vertices should be on the implicit")
    def test_normal_has_correct_dimensions(self):
        new_opt = MeshOptimizer()
        new_opt.load_example('bowl_15_holes')
        point = np.array([1, 1, 1, 1])
        self.assertEqual(new_opt.getNormalVectorAtPoint(point).shape,(3,3))
    def tearDown(self):
        pass


if __name__== "__main__":
    unittest.main()
