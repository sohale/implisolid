import numpy as np

VERBOSE = False


def centroids(vertex, faces):
    c = np.mean(vertex[faces[:], :], axis=1)
    return c


def make_neighbour_faces_of_vertex(faces):
    """ neighbour_faces_of_vertex is a list. index=vertex, v1,v2,v3 """
    neighbour_faces_of_vertex = {}  # np.zeros( (vertex_count,3) , dtype=np.type) - 1
    for fi in range(faces.shape[0]):
        for vi in range(3):
            v1 = faces[fi, vi]
            if v1 not in neighbour_faces_of_vertex:
                neighbour_faces_of_vertex[v1] = [fi]
            else:
                assert fi not in neighbour_faces_of_vertex[v1]
                neighbour_faces_of_vertex[v1].append(fi)

    return neighbour_faces_of_vertex


def make_edge_lookup_old(faces):
    """ """

    print "Warning: not tested enough"
    nfaces = faces.shape[0]
    # assert nfaces % 2 == 0
    num_edges = nfaces * 3 / 2
    if True:    # VERBOSE:
        print("nfaces= ", nfaces, "num_edges:", num_edges)

    edges_of_faces = np.zeros((nfaces, 3), dtype=np.int) - 1
    faces_of_edges = np.zeros((num_edges, 2), dtype=np.int) - 1
    vertpairs_of_edges = np.zeros((num_edges,), dtype=np.int) - 1

    modulo = long(num_edges)
    eulookup = -np.ones((modulo * num_edges + num_edges,))

    edge_counter = 0
    for fi in range(len(faces)):
        for vj in range(3):
            if VERBOSE:
                print("------------")
            v2j = (vj + 1) % 3
            e = (faces[fi, vj], faces[fi, v2j])
            eu = (e[0], e[1]) if e[1] > e[0] else (e[1], e[0])
            assert e[0] >= 0
            assert e[1] >= 0
            assert not e[1] == e[0]
            eu_pair_int = int(eu[0] + eu[1] * modulo)
            eu_pair_int_signed = eu_pair_int

            if VERBOSE:
                print("eu_pair_int = eu[0] + eu[1] * modulo", eu[0], eu[1], " -> ", eu_pair_int)
                print(" e=", e, " eu_pair_int=", eu_pair_int)

            if eulookup[eu_pair_int] < 0:
                new_edge = True
            else:
                new_edge = False

            if VERBOSE:
                print("new_edge= ", new_edge, " edge_counter=", edge_counter)

            if new_edge:
                e_id = edge_counter
                edges_of_faces[fi, vj] = e_id

                faces_of_edges[e_id, 0] = fi
                assert not vj == v2j
                assert vertpairs_of_edges[e_id] == -1
                vertpairs_of_edges[e_id] = np.abs(eu_pair_int_signed)  # eu

                assert eulookup[eu_pair_int] == -1, "        %d " % (eulookup[eu_pair_int],)
                eulookup[eu_pair_int] = e_id

                edge_counter += 1

            else:
                e_id = eulookup[eu_pair_int]
                assert e_id >= 0
                edges_of_faces[fi, vj] = e_id
                faces_of_edges[e_id, 1] = fi
                assert not vj == v2j

                other_fi = faces_of_edges[e_id, 0]
                if VERBOSE:
                    print(vj, faces[fi, :], faces[other_fi, :])

                if VERBOSE:
                    print(vertpairs_of_edges[e_id], eu_pair_int_signed)
                assert vertpairs_of_edges[e_id] == -eu_pair_int_signed or \
                    vertpairs_of_edges[e_id] == +eu_pair_int_signed     # eu
                assert eulookup[eu_pair_int] == e_id

            if VERBOSE:
                print("edges_of_faces ", edges_of_faces)
                print("faces_of_edges ", faces_of_edges)
                print("eulookup ", eulookup)
            if True:
                eu_paired_int = vertpairs_of_edges[e_id]
                (v1, v2) = (eu_paired_int % modulo, eu_paired_int / modulo)
                if VERBOSE:
                    print("vertpair:", eu_paired_int, " -> ", v1, v2)

    for fi in range(len(faces)):
        e123 = edges_of_faces[fi, :]
        assert e123.size == 3, "Polygons > 3 not allowed"
        e0, e1, e2 = e123[0], e123[1], e123[2]
        if VERBOSE:
            print("edges_of_faces[%d]=" % (fi,), e123)
            print(faces_of_edges[e0, :], " , ", faces_of_edges[e1, :], " , ", faces_of_edges[e2, :])

        assert fi in faces_of_edges[e0, :]
        assert fi in faces_of_edges[e1, :]
        assert fi in faces_of_edges[e2, :]
        va = faces[fi, :]

        e_id0 = e0
        e_id1 = e1
        e_id2 = e2
        vertpairs_of_edges[e_id0] == va[0]
        vertpairs_of_edges[e_id1] == va[1]
        vertpairs_of_edges[e_id2] == va[2]

    for e_id in range(num_edges):
        eu_paired_int = np.abs(vertpairs_of_edges[e_id])
        (v1, v2) = (eu_paired_int % modulo, eu_paired_int / modulo)
        assert eu_paired_int == v1 + v2 * modulo
        assert eu_paired_int >= 0
        if VERBOSE:
            print(eu_paired_int, eulookup[eu_paired_int], e_id)
        assert np.abs(eulookup[eu_paired_int]) == e_id

    assert np.all(np.ravel(edges_of_faces) > -1)
    return (edges_of_faces, faces_of_edges, vertpairs_of_edges)

import scipy.sparse as sp


def make_edge_lookup_sparse(faces):
    """ does make_edge_lookup_old() using sparse matrices. An alternatve solution would be using adictionary."""
    # Todo: Rewrite using asparse matrix (e.g. DOK) as Adjacency matrix. In that case, `modulo` will not be used.
    print "Sparse verion:"

    nfaces = faces.shape[0]
    assert nfaces % 2 == 0
    num_edges = nfaces * 3 / 2
    if True:    # VERBOSE:
        print("nfaces= ", nfaces, "num_edges:", num_edges)

    edges_of_faces = np.zeros((nfaces, 3), dtype=np.int) - 1
    faces_of_edges = np.zeros((num_edges, 2), dtype=np.int) - 1
    vertpairs_of_edges = np.zeros((num_edges,), dtype=np.int) - 1

    modulo = long(num_edges)    # *2
    lookup_array_size = modulo * num_edges + num_edges
    eulookup = sp.dok_matrix((lookup_array_size, 1), dtype=int)
    edge_counter = 0
    for fi in range(len(faces)):
        for vj in range(3):
            if VERBOSE:
                print("------------")
            v2j = (vj + 1) % 3
            e = (faces[fi, vj], faces[fi, v2j])
            # eu = (vj, v2j) if vj>v2j else (v2j, vj)  #unique edge id
            eu = (e[0], e[1]) if e[1] > e[0] else (e[1], e[0])  # unique edge id
            assert e[0] >= 0
            assert e[1] >= 0
            assert not e[1] == e[0]
            eu_pair_int = int(eu[0] + eu[1] * modulo)

            eu_pair_int_signed = eu_pair_int

            if VERBOSE:
                print("eu_pair_int = eu[0] + eu[1] * modulo", eu[0], eu[1], " -> ", eu_pair_int)
                # dont confuse e_id and eu_pair_int. vertpairs_of_edges[e_id] == eu_pair_int
                print(" e=", e, " eu_pair_int=", eu_pair_int)

            # print("********** ", eu_pair_int)
            if eulookup[eu_pair_int, 0]-1 < 0:
                new_edge = True
            else:
                new_edge = False

            if VERBOSE:
                print("new_edge= ", new_edge, " edge_counter=", edge_counter)

            if new_edge:
                # vertpairs_of_edges = ***

                # add a new edge
                e_id = edge_counter
                edges_of_faces[fi, vj] = e_id

                # print faces_of_edges.shape, e_id, edge_counter, num_edges
                faces_of_edges[e_id, 0] = fi
                assert not vj == v2j
                # assert vertpairs_of_edges[e_id] == eu_pair_int, "vertpairs_of_edges[e_id] == eu_pair_int:  vertpairs_of_edges[%d]=%d == %d"%(e_id, vertpairs_of_edges[e_id], eu_pair_int)
                assert vertpairs_of_edges[e_id] == -1
                vertpairs_of_edges[e_id] = np.abs(eu_pair_int_signed)  # eu
                # something wrong here.

                assert eulookup[eu_pair_int, 0]-1 == -1, "        %d " % (eulookup[eu_pair_int, 0],)
                eulookup[eu_pair_int, 0] = e_id + 1    # edge number?

                edge_counter += 1

            else:
                e_id = eulookup[eu_pair_int, 0]-1
                assert e_id >= 0
                edges_of_faces[fi, vj] = e_id
                faces_of_edges[e_id, 1] = fi
                assert not vj == v2j

                other_fi = faces_of_edges[e_id, 0]
                if VERBOSE:
                    print(vj, faces[fi, :], faces[other_fi, :])

                if VERBOSE:
                    print(vertpairs_of_edges[e_id], eu_pair_int_signed)
                assert vertpairs_of_edges[e_id] == -eu_pair_int_signed or \
                    vertpairs_of_edges[e_id] == +eu_pair_int_signed     # eu
                assert eulookup[eu_pair_int, 0] == e_id + 1

            if VERBOSE:
                print("edges_of_faces ", edges_of_faces)
                print("faces_of_edges ", faces_of_edges)
                # print("vertpairs_of_edges ", vertpairs_of_edges)
                print("eulookup ", eulookup-1)
            if True:
                eu_paired_int = vertpairs_of_edges[e_id]
                (v1, v2) = (eu_paired_int % modulo, eu_paired_int / modulo)
                if VERBOSE:
                    print("vertpair:", eu_paired_int, " -> ", v1, v2)

    for fi in range(len(faces)):
        e123 = edges_of_faces[fi, :]
        assert e123.size == 3, "Polygons > 3 not allowed"
        e0, e1, e2 = e123[0], e123[1], e123[2]
        if VERBOSE:
            print("edges_of_faces[%d]=" % (fi,), e123)
            print(faces_of_edges[e0, :], " , ", faces_of_edges[e1, :], " , ", faces_of_edges[e2, :])

        assert fi in faces_of_edges[e0, :]
        assert fi in faces_of_edges[e1, :]
        assert fi in faces_of_edges[e2, :]
        va = faces[fi, :]

        e_id0 = e0
        e_id1 = e1
        e_id2 = e2
        vertpairs_of_edges[e_id0] == va[0]
        vertpairs_of_edges[e_id1] == va[1]
        vertpairs_of_edges[e_id2] == va[2]

    for e_id in range(num_edges):
        eu_paired_int = np.abs(vertpairs_of_edges[e_id])
        (v1, v2) = (eu_paired_int % modulo, eu_paired_int / modulo)
        assert eu_paired_int == v1 + v2 * modulo
        assert eu_paired_int >= 0
        if VERBOSE:
            print(eu_paired_int, eulookup[eu_paired_int, 0]-1, e_id)
        assert np.abs(eulookup[eu_paired_int, 0]) == e_id+1

    assert np.all(np.ravel(edges_of_faces) > -1)
    # edges_of_faces : index = face number, value = [edge number, edge number, edge number]
    # faces_of_edges : index = edge number, value = face number, face number
    # vertpairs_of_edges : index = edge number, value = eu_paired_int([vertex1, vertex2], here vertex1 < vertex2)
    # eulookup[eu_paired_int, 0] = edge number index + 1
    return (edges_of_faces, faces_of_edges)


if __name__ == '__main__':
    pass
