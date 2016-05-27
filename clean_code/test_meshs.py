import unittest

import numpy as np

from vtk_mc import vtk_mc

import example_objects
import sys

B = 1000000L
CHECK_PAIRED = True

# test the correctness of all meshs, work for all objects except screw1 and screw2


def check_mesh(faces):
    # based on the function coming from ohtake_belyaev_5.py
    """ Checks if the output of the Marching Cubes and subdivision are correct. Correction of the mesh: closedness, each edge appears exactly twice, etc"""
    from mesh_utils import make_edge_lookup_sparse
    check_faces(faces)
    check_faces3(faces)

    # The following does a series of `assert`s which check the correctness of the mesh
    (edges_of_faces, faces_of_edges, vertpairs_of_edges) = \
        make_edge_lookup_sparse(faces)


def check_faces3(facets):
    REMOVE_REPEATED_EDGES(facets)
    pass


def get_edge_code_triples_of_mesh(facets):
    """ Returns an array of (F)x(3), containing the 'edge codes' of sides of the faces of a mesh.
    There are 3 sides for each face.
    An 'edge code' is a long integer (int64) v1+B*v2 where v1,v2 are the indices of the ends (vertices) of the edge, where v1<v2."""
    e0 = facets[:, np.newaxis, [0, 1]]
    e1 = facets[:, np.newaxis, [1, 2]]
    e2 = facets[:, np.newaxis, [2, 0]]   # np view
    e012 = np.concatenate((e0, e1, e2), axis=1)  # n x 3 x 2
    assert e012.base is None  # make sure it's not a view of faces
    e012.sort(axis=2)
    BB = np.array([1L, B], dtype=np.int64)
    all_edges_triples = np.dot(e012, BB)  # n x 3
    assert all_edges_triples.dtype == np.int64
    assert all_edges_triples.size == 0 or np.min(all_edges_triples) >= 0
    assert np.max(facets, axis=None) < B
    assert all_edges_triples.size == 0 or np.min(all_edges_triples) >= 0
    assert all_edges_triples.shape == (facets.shape[0], 3)
    return all_edges_triples, e012


def REMOVE_REPEATED_EDGES(faces):
    # see check_faces()

    if True:
        # >begin of refactorable region
        f0 = faces[:, np.newaxis, 0:2]
        f1 = faces[:, np.newaxis, 1:3]
        f2 = faces[:, np.newaxis, [0, 2]]
        f0 = f0.copy()
        f0.sort(axis=2)  # changes the order in faces!
        f1 = f1.copy()
        f1.sort(axis=2)
        f2 = f2.copy()
        f2.sort(axis=2)
        fe3 = np.concatenate((f0, f1, f2), axis=1)  # shape==(:,3,2)

        BB = np.array([1L, B], dtype=np.int64)
        edg = np.dot(fe3, BB)   # fx3
        assert edg.dtype == np.int64
        assert edg.size == 0 or np.min(edg) >= 0
        assert np.max(faces, axis=None) < B
        # <end of refactorable region

    edg_1, fe3_1 = get_edge_code_triples_of_mesh(faces)
    assert np.allclose(edg_1, edg)
    assert np.allclose(fe3_1, fe3)

    # Sort edges to detect repeated edges. Each edge should appear exactly twice.
    # was q
    edg_sorted = edg.ravel()
    sort_idx = edg_sorted.argsort()
    assert sort_idx.shape == (faces.shape[0]*3,)

    q_unsorted = edg_sorted.copy()
    edg_sorted.sort()
    assert np.all(q_unsorted[sort_idx] == edg_sorted)

    assert np.all(np.diff(edg_sorted)[0::2] == 0), "some edges are not exactly repeated once: not manifold."
    assert np.all(np.diff(edg_sorted)[1::2] != 0), "some edges are not exactly repeated once: not manifold."

    assert edg_sorted.shape[0] % 3 == 0

    # back to pairs
    return

    # idx_xy is a tuple
    if not np.all(np.diff(edg_sorted)[::2] == 0):
        pass


def face_triplet_codes(faces):
    """ returns an array of (f) x 1, with codes such as '169000044000031' of type long (np.int64), one single code for each unique face. """
    f3sides = faces.copy()
    f3sides.sort(axis=1)

    BBB = np.array([[1L, B, B*B]], dtype=np.int64).transpose().ravel()  # (3,)

    d = np.dot(f3sides, BBB)
    assert d.dtype == np.int64
    assert d.size == 0 or np.min(d) >= 0
    assert np.max(faces.ravel()) < B
    del f3sides

    return d


def check_face_triplets(faces):
    # unique faces
    d = face_triplet_codes(faces)
    face_triplet_ids = d.ravel().copy()

    face_order = face_triplet_ids.argsort()
    face_triplet_ids.sort()
    # Check there is no repeated faces (with exact same vertex list):
    diff0 = (np.diff(face_triplet_ids) == 0)
    if np.sum(diff0) != 0:
        nonz = np.nonzero(diff0)[0]
        diff01 = diff0.copy()
        diff01[nonz+1] = True

    # diff0 versus diff01: diff01 is for print only: to print both sides (Elements) of each "diff==0"
    # nonz = np.nonzero(diff01)[0]  #both of them
    nonz = np.nonzero(diff0)[0]  # only the redundants
    bad_faces = nonz
    # but some repeated ones may remain
    # howver the original idx (of redundant faces) are = face_order[nonz]
    original_indices = face_order[bad_faces]

    assert np.sum(diff0) == original_indices.size, "Number of redundant faces"

    return original_indices


def check_faces(faces):
    # print("------ check_faces(faces)")

    redundant_faces = check_face_triplets(faces)

    assert redundant_faces.size == 0, "Repeated faces found"

    # edg, fe3  = get_edge_code_triples_of_mesh(faces)
    """cannot refactor: fe3 is used """
    if True:
        # >begin of refactorable region
        # unique edges
        f0 = faces[:, np.newaxis, 0:2]
        f1 = faces[:, np.newaxis, 1:3]
        f2 = faces[:, np.newaxis, [0, 2]]
        f0 = f0.copy()
        f0.sort(axis=2)  # changes the order in faces!
        f1 = f1.copy()
        f1.sort(axis=2)
        f2 = f2.copy()
        f2.sort(axis=2)

        fe3 = np.concatenate((f0, f1, f2), axis=1)  # shape==(:,3,2)
        BB = np.array([[1L, B]], dtype=np.int64).transpose().ravel()  # 2x-
        edg = np.dot(fe3, BB)  # fx3
        assert edg.dtype == np.int64
        assert edg.size == 0 or np.min(edg) >= 0
        assert np.max(faces, axis=None) < B
        # < end of refactorable region

    edg_1, fe3_1 = get_edge_code_triples_of_mesh(faces)
    assert np.allclose(edg_1, edg)
    assert np.allclose(fe3, fe3_1)

    # Sort edges to detect repeated edges. Each edge should appear exactly twice.
    q = edg.ravel().copy()
    sort_idx = q.argsort()
    assert sort_idx.shape == (faces.shape[0]*3,)
    q_unsorted = q.copy()
    q.sort()
    assert np.all(q_unsorted[sort_idx] == q)

    # You can get the sorted inices here:
    if not np.all(np.diff(q)[::2] == 0):

        i1 = (np.diff(q)[::2] != 0)

        if False:
            exit()
            del i1

    if CHECK_PAIRED:
        assert np.all(np.diff(q)[::2] == 0), "Not all edges are paired"


class ImplicitFunctionTests(unittest.TestCase):

    def test_mesh_correctness(self):

        examples_list = example_objects.get_all_examples([2])
        for example_name in examples_list:
            print("example_name = ", example_name)
            iobj = example_objects.make_example_vectorized(example_name)
            from example_objects import make_example_vectorized
            iobj = make_example_vectorized(example_name)
            self.check_mesh(iobj, objname=example_name)

    def check_mesh(self, iobj, objname=None):

        """Do the centroids projection """
        if iobj is not None:
            (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-3, +5, 0.2)

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

            check_mesh(faces)


if __name__ == '__main__':
    unittest.main()
