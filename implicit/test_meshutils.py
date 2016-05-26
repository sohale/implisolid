import unittest

#from lettuce import *
#from ipdb import set_trace
#import ipdb


def check_if_imported(module_name):
    from sys import modules
    try:
        module = modules[module_name]
        return True
    except KeyError:
        #__import__('m')
        return False


def check_imdb_import():
    assert not check_if_imported("ipdb"), "Should import `imdb` in unittest"


def show_imdb_import():
    print "*is imported." if check_if_imported("ipdb") else "*not imported"


show_imdb_import()
from mesh_utils import make_sparse_neighbour_faces_of_vertex_csr
show_imdb_import()
from mesh_examples_for_tests import testcase_square, testcase_cube
show_imdb_import()

import numpy as np


class A(unittest.TestCase):

    def test_make_sparse_neighbour_faces_of_vertex_csr(self):
        _v, faces = testcase_square()
        sparse_matrix = make_sparse_neighbour_faces_of_vertex_csr(faces)
        self.assertTrue(sparse_matrix.shape[0] == 4)

        # Make sure we are testing the right thing
        fcs1 = np.array([[0, 1, 2], [0, 2, 3]])
        assert np.allclose(faces, fcs1)
        self.assertTrue(1)
        correct_fov = np.array([[1, 1], [1, 0], [1, 1], [0, 1]])
        self.assertTrue(np.allclose(correct_fov, sparse_matrix.todense()))

        self.full_vof_test(faces, sparse_matrix)

    def full_vof_test(self, faces, vof_matrix_sparse):

        self.assertTrue(issubclass(vof_matrix_sparse.dtype.type, np.integer))

        VERTS_DIM = 0
        FACES_DIM = 1
        faces_per_vert = np.sum(vof_matrix_sparse.todense(), axis=VERTS_DIM)
        self.assertTrue(np.allclose(faces_per_vert, 3))

        # Faces that have that vertex. Not 3.
        verts_per_face = np.sum(vof_matrix_sparse.todense(), axis=FACES_DIM)
        # print faces_per_vert,
        # print verts_per_face
        num_triangles = faces.shape[0]
        self.assertTrue((np.sum(verts_per_face)+0.) / 3. == num_triangles)

        # Master test:
        recreate_faces = np.nonzero(vof_matrix_sparse.T > 0)[1].reshape((num_triangles, 3))
        # not tested
        self.assertTrue(np.allclose(recreate_faces, faces))  # Untimate test!


if __name__ == '__main__':
    check_imdb_import()
    unittest.main()
