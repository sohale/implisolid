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
        nverts = 4
        self.assertTrue(sparse_matrix.shape[0] == nverts)

        # Make sure we are testing the right thing
        fcs1 = np.array([[0, 1, 2], [0, 2, 3]])
        assert np.allclose(faces, fcs1)
        self.assertTrue(1)
        correct_fov = np.array([[1, 1], [1, 0], [1, 1], [0, 1]])
        # Is there any way to avoid todense()?
        self.assertTrue(np.allclose(correct_fov, sparse_matrix.todense()))

        self.full_vof_test(faces, sparse_matrix)

    def full_vof_test(self, faces, vof_matrix_sparse):

        self.assertTrue(issubclass(vof_matrix_sparse.dtype.type, np.integer))

        VERTS_DIM = 0
        FACES_DIM = 1
        faces_per_vert = vof_matrix_sparse.sum(axis=VERTS_DIM)
        self.assertTrue(np.allclose(faces_per_vert, 3))

        # Faces that have that vertex. Not 3.
        verts_per_face = vof_matrix_sparse.sum(axis=FACES_DIM)
        # print faces_per_vert,
        # print verts_per_face
        num_triangles = faces.shape[0]
        self.assertTrue((np.sum(verts_per_face)+0.) / 3. == num_triangles)

        # Master test:
        recreate_faces0 = np.nonzero(vof_matrix_sparse.T > 0)[1].reshape((num_triangles, 3))

        # We need row numbers (i.e. vertex indices) of the `fov` matrix.
        # But we are interested in columns of vof_matrix_sparse.
        # cols_a contains the column numbers of all non-zero elements of fov.T
        # Doing .T is necessary because nonzero() is row-wise.

        # print vof_matrix_sparse.todense()

        (rows_a, cols_a) = np.nonzero(vof_matrix_sparse.T > 0)
        # print (rows_a, cols_a)
        # cols_a: v1,v2,v3, v1,v2,v3, v1,v2,v3, ...
        recreate_faces = cols_a.reshape((num_triangles, 3), order='C')

        recreate_faces.sort(axis=1)
        faces_c = faces
        faces_c.sort(axis=1)

        self.assertTrue(np.allclose(recreate_faces, faces_c))  # Untimate test!

    def test_2(self):
        _v, faces = testcase_cube()
        sparse_matrix = make_sparse_neighbour_faces_of_vertex_csr(faces)
        nverts = _v.shape[0]
        self.assertTrue(sparse_matrix.shape[0] == nverts)

        #print faces
        #print sparse_matrix.T.todense()
        self.full_vof_test(faces, sparse_matrix)

if __name__ == '__main__':
    check_imdb_import()
    unittest.main()
