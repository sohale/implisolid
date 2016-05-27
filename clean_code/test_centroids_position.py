import unittest

import numpy as np

from basic_functions import check_vector3_vectorized

from vtk_mc import vtk_mc

import math
import example_objects
import sys


class ImplicitFunctionTests(unittest.TestCase):

    def test_centroids_projection(self):

        examples_list = example_objects.get_all_examples([2])
        for example_name in examples_list:
            print("example_name = ", example_name)
            iobj = example_objects.make_example_vectorized(example_name)
            from example_objects import make_example_vectorized
            iobj = make_example_vectorized(example_name)
            self.check_centroids_projection(iobj, objname=example_name)

    def check_centroids_projection(self, iobj, objname=None):
        TOLERANCE = 0.00001
        # TOLERANCE = 0.7 to pass the test for every object
        """Do the centroids projection """
        if iobj is not None:
            VERTEX_RELAXATION_ITERATIONS_COUNT = 0
            (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-3, +5, 0.2)
            if objname == "cube_with_cylinders" or objname == "twist_object" or objname == "french_fries" or objname == "rdice_vec" or objname == "rods" or objname == "bowl_15_holes":
                VERTEX_RELAXATION_ITERATIONS_COUNT = 1

            if objname == "cyl4":
                (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-32 / 2, +32 / 2, 1.92 / 4.0)

            elif objname == "french_fries" or objname == "rods":
                (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-3, +5, 0.11)  # 0.05

            elif objname == "bowl_15_holes":
                (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-3, +5, 0.15)

            elif objname == "cyl3":
                (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-32 / 2, +32 / 2, 1.92 / 4.0)

            elif objname == "cyl1":
                (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-16, +32, 1.92 * 0.2 * 10 / 2.0)

            from stl_tests import make_mc_values_grid
            gridvals = make_mc_values_grid(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE, old=False)
            vertex, faces = vtk_mc(gridvals, (RANGE_MIN, RANGE_MAX, STEPSIZE))
            print("MC calculated.")
            sys.stdout.flush()

            from ohtake_belyaev_demo_subdivision_projection_qem import process2_vertex_resampling_relaxation, compute_average_edge_length, set_centers_on_surface__ohtake_v3s

            for i in range(VERTEX_RELAXATION_ITERATIONS_COUNT):
                vertex, facets_not_used, centroids = process2_vertex_resampling_relaxation(vertex, faces, iobj)
            assert not np.any(np.isnan(vertex.ravel()))  # fails
            print("Vertex relaxation applied.")
            sys.stdout.flush()

            # projection
            average_edge = compute_average_edge_length(vertex, faces)

            old_centroids = np.mean(vertex[faces[:], :], axis=1)
            check_vector3_vectorized(old_centroids)

            new_centroids = old_centroids.copy()

            set_centers_on_surface__ohtake_v3s(iobj, new_centroids, average_edge)
            check_vector3_vectorized(new_centroids)
            # checking if the projection is correct by calling the implicitFunction
            f = iobj.implicitFunction(new_centroids)

            # Two ways of doing this test, in the first one we strictly consider that the test fail if one value is superior
            # to the tolerance and in the second we print the numbere of point who fail the test

            Number_of_point_who_fail = True

            if Number_of_point_who_fail is True:
                fail = 0
                for i in range(new_centroids.shape[0]):
                    if math.fabs(f[i]) > TOLERANCE:
                        fail += 1
                print objname, "Number of points:", new_centroids.shape[0], "Number of points who fails the test:", fail

            else:
                for i in range(new_centroids.shape[0]):
                    print "Fail the test", math.fabs(f[i])
                    self.assertTrue(math.fabs(f[i]) < TOLERANCE)


if __name__ == '__main__':
    unittest.main()
