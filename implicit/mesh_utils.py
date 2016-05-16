import numpy as np
from basic_types import check_vector4_vectorized
from ipdb import set_trace
import implicit_config

VERBOSE = False

# class Mesh(object):
#    pass


def centroids(verts, faces):
    # print("faces: ", faces.shape)
    # print(verts[faces[:],:].shape)
    c = np.mean(verts[faces[:], :], axis=1)
    # print(c.shape)
    return c


def face_area(verts, faces):
    """incomplete"""
    i = 30
    # verts[faces,:].shape   # [numfaces,3,3]
    print(verts[faces[:], :])
    abc = verts[faces[i, :], :]
    print (abc.shape)
    # exit()


def is_face_flipped(faces, verts, face_idx):
    pass


def flip_face(faces, verts, face_idx):
    # inplace
    pass


import time

class Timer:
    def __enter__(self):
        self.start = time.clock()
        return self

    def __exit__(self, *args):
        self.end = time.clock()
        self.interval = self.end - self.start
class TooMemoryIntensive(Exception):
    pass


def make_edge_lookup_old(faces):
    """ """
    # edges_of_faces : index = face number, value = [edge number, edge number, edge number]
    # faces_of_edges : index = edge number, value = face number1, face number2
    # vertpairs_of_edges : index = edge number, value = eu_paired_int([vertex1, vertex2], here vertex1 < vertex2)
    # eulookup[eu_paired_int] = edge number index
    #raise
    print "Warning: not tested enough"

    nfaces = faces.shape[0]
    assert nfaces % 2 == 0
    num_edges = nfaces * 3 / 2
    if True:    # VERBOSE:
        print("nfaces= ", nfaces, "num_edges:", num_edges)

    edges_of_faces = np.zeros((nfaces, 3), dtype=np.int) - 1
    faces_of_edges = np.zeros((num_edges, 2), dtype=np.int) - 1
    vertpairs_of_edges = np.zeros((num_edges,), dtype=np.int) - 1

    #set_trace()
    #MAX_MATRIX_SIZE = 50000000

    modulo = long(num_edges)    # *2
    lookup_array_size = modulo * num_edges + num_edges
    print "lookup_array_size",lookup_array_size
    if lookup_array_size > implicit_config.MAX_MATRIX_SIZE:
        raise TooMemoryIntensive()
    eulookup = -np.ones((lookup_array_size,), dtype=int)

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
            eu_pair_int_revsersed = int(eu[1] + eu[0] * modulo)
            eu_pair_int_signed = eu_pair_int
            if not e[1] > e[0]:
                # eu_pair_int_signed = -eu_pair_int
                eu_pair_int_signed = eu_pair_int
            eu_pair_int_signed = eu_pair_int

            if VERBOSE:
                print("eu_pair_int = eu[0] + eu[1] * modulo", eu[0], eu[1], " -> ", eu_pair_int)
                # dont confuse e_id and eu_pair_int. vertpairs_of_edges[e_id] == eu_pair_int
                print(" e=", e, " eu_pair_int=", eu_pair_int)

            # print("********** ", eu_pair_int)
            if eulookup[eu_pair_int] < 0:
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

                #print faces_of_edges.shape, e_id, edge_counter, num_edges
                faces_of_edges[e_id, 0] = fi
                assert not vj == v2j
                # assert vertpairs_of_edges[e_id] == eu_pair_int, "vertpairs_of_edges[e_id] == eu_pair_int:  vertpairs_of_edges[%d]=%d == %d"%(e_id, vertpairs_of_edges[e_id], eu_pair_int)
                assert vertpairs_of_edges[e_id] == -1
                vertpairs_of_edges[e_id] = np.abs(eu_pair_int_signed)  # eu
                # something wrong here.

                assert eulookup[eu_pair_int] == -1, "        %d " % (eulookup[eu_pair_int],)
                eulookup[eu_pair_int] = e_id    # edge number?

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
                # print("vertpairs_of_edges ", vertpairs_of_edges)
                print("eulookup ", eulookup)
            if True:
                eu_paired_int = vertpairs_of_edges[e_id]
                (v1, v2) = (eu_paired_int % modulo, eu_paired_int / modulo)
                if VERBOSE:
                    print("vertpair:", eu_paired_int, " -> ", v1, v2)
            # if VERBOSE and new_edge:
            #    import os
            #    os.system("pause")

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
        # *********

        # def e_id(e_tuple):
        #    eu = (e_tuple[0], e_tuple[1]) if e_tuple[1]>e_tuple[0] else (e_tuple[1], e_tuple[0])  #unique edge id
        #    eu_pair_int = eu[0] + eu[1] * modulo
        #    return eu_pair_int
        # e_id0 = e_id( e0 )
        # e_id1 = e_id( e1 )
        # e_id2 = e_id( e2 )
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
        # face[*]=*
        assert eu_paired_int >= 0
        if VERBOSE:
            print(eu_paired_int,eulookup[eu_paired_int], e_id)
        assert np.abs(eulookup[eu_paired_int]) == e_id

    assert np.all(np.ravel(edges_of_faces) > -1)
    # edges_of_faces : index = face number, value = [edge number, edge number, edge number]
    # faces_of_edges : index = edge number, value = face number, face number
    # vertpairs_of_edges : index = edge number, value = eu_paired_int([vertex1, vertex2], here vertex1 < vertex2)
    # eulookup[eu_paired_int] = edge number index
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
    #print "lookup_array_size",lookup_array_size
    #if lookup_array_size > implicit_config.MAX_MATRIX_SIZE:
    #    raise TooMemoryIntensive()
    #eulookup = np.zeros((lookup_array_size,), dtype=int)
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
            eu_pair_int_revsersed = int(eu[1] + eu[0] * modulo)
            eu_pair_int_signed = eu_pair_int
            if not e[1] > e[0]:
                # eu_pair_int_signed = -eu_pair_int
                eu_pair_int_signed = eu_pair_int
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

                #print faces_of_edges.shape, e_id, edge_counter, num_edges
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
    return (edges_of_faces, faces_of_edges, vertpairs_of_edges)

"""
def make_edge_lookup(faces):
    with Timer() as t1:
            (edges_of_faces1, faces_of_edges1, vertpairs_of_edges1) = make_edge_lookup_old(faces)
    with Timer() as t2:
            (edges_of_faces2, faces_of_edges2, vertpairs_of_edges2) = make_edge_lookup_sparse(faces)
    # 2.3 - 2.5 times slower
    print("Projection two methods: done within ", t1.interval, t2.interval, "RATIO =", t2.interval/t1.interval)

    assert np.allclose(edges_of_faces1, edges_of_faces2)
    assert np.allclose(faces_of_edges1, faces_of_edges2)
    assert np.allclose(vertpairs_of_edges1, vertpairs_of_edges2)
    set_trace()
    return (edges_of_faces2, faces_of_edges2, vertpairs_of_edges2)
"""

make_edge_lookup = make_edge_lookup_sparse


def invariants():
    pass


def mesh_invariant(faces):
    # gets some mesh metrics. Useful for invariants.
    edge_lookup = {}
    for fi in range(faces.shape[0]):
        for vi in range(3):
            (v1, v2) = (faces[fi, vi], faces[fi, (vi + 1) % 3])
            if v1 > v2:
                (v1, v2) = (v2, v1)
            key = "%d-%d" % (v1, v2,)
            # print(key)
            if key in edge_lookup:
                edge_lookup[key] += 1
            else:
                edge_lookup[key] = 1
    odd_edges = 0
    for k, v in edge_lookup.iteritems():
        if v == 2:
            pass
        elif v == 1:
            odd_edges += 1
        else:
            print("****")
            print(k, v)
    open_edges = odd_edges
    print("open edges: ", open_edges)
    if open_edges > 0:
        print("Warning: shape is not closed.")
    # invariant: odd_eges should be zero


def test_make_edge_lookup():
    faces = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
    (edges_of_faces, faces_of_edges, vertpairs_of_edges) = make_edge_lookup(faces)
    pass


def neighbour_faces(faces):
    # The list of faces neighbour to a face
    e123 = edges_of_faces[fi]
    faces_of_edges[e123[0]]
    pass


def make_dual(verts, faces):
    """ neighbour_faces[fi] = index=face, Fx3 : [f1,f2,f3] , i.e. array of indices of the neighbour faces
        neighbour_faces_of_vertex[] , index=vertex
    """

    (edges_of_faces, faces_of_edges, vertpairs_of_edges) = make_edge_lookup(faces)
    dual_verts = centroids(verts, faces)
    dual_faces = None
    for fi in range(faces.shape[0]):
        e123 = edges_of_faces[fi, :]
        f1pair = faces_of_edges[e123[0], :]
        f2pair = faces_of_edges[e123[1], :]
        f3pair = faces_of_edges[e123[2], :]
        print(f1pair, f2pair, f3pair)
        print("incomplete")

# The only useful method in our method. This needs to be done once and can be reused.


def make_neighbour_faces_of_vertex(faces):
    """ neighbour_faces_of_vertex is a list. index=vertex, v1,v2,v3 """
    vertex_count = np.max(np.max(faces)) + 1
    neighbour_faces_of_vertex = {}  # np.zeros( (vertex_count,3) , dtype=np.type) - 1
    #set_trace()
    for fi in range(faces.shape[0]):
        for vi in range(3):
            v1 = faces[fi, vi]
            if v1 not in neighbour_faces_of_vertex:
                neighbour_faces_of_vertex[v1] = [fi]
            else:
                assert fi not in neighbour_faces_of_vertex[v1]
                neighbour_faces_of_vertex[v1].append(fi)
        #print len(neighbour_faces_of_vertex),
    # print(neighbour_faces_of_vertex)
    # print(map( lambda k: len(neighbour_faces_of_vertex[k]), neighbour_faces_of_vertex   )) #for k,v in neighbour_faces_of_vertex:
    #    print(len(v), end="")
    return neighbour_faces_of_vertex


def test_neighbours():
    """ tests if all those faces have that vertex. """
    pass

import vectorized

# some modification
# def optimize_points1(iobj, xa0, tolerance = 0.00001, lambda_=0.1):
#    """ not inplace anymore"""
#    assert vectorized.is_implicit_type(iobj)
#    if not not_using_op: #default is False in the wrong version.
#        xa = xa0.copy()
# ...

# old working version.
# inplace = optimize_points1_not_using_op = True  #good screw


def optimize_points1_inplace_oldworking(iobj, xa, tolerance=0.00001, lambda_=0.1, inplace=True, maxdist=None):

    if maxdist is not None:
        print "Warning: maxdist is not used in the selected algorithm"

    assert vectorized.is_implicit_type(iobj)
    print("--------optimize_points1_inplace_oldworking-------")
    print(inplace)
    if not inplace:     # optimize_points1_not_using_op:
        xa = xa.copy()
    npoints = xa.shape[0]
    assert xa.shape[1] == 4
    vectorized.check_vector4_vectorized(xa)

    ia = np.arange(0, npoints)
    # xa_ = xa[ia, :] [:,0:3]

    va = np.ones((npoints,))
    ga = np.ones((npoints, 4))

    while True:
        # evaluate
        # va = iobj.implicitFunction(xa[ia, :])
        va[ia] = iobj.implicitFunction(xa[ia, :])
        erra = np.abs(va[ia])

        # decide
        ia2 = ia[erra > tolerance]
        ia2_accepted = ia[erra <= tolerance]

        if ia2.size == 0:
            assert np.all(iobj.implicitFunction(xa) <= tolerance)
            return xa

        # update: accepted (output)
        # xa[ia2_accepted, 0:3] = xa[ia2_accepted, :]

        # update: keep updating (loop)
        g1 = iobj.implicitGradient(xa[ia2, :])[:, 0:3]
        g1 = g1 / np.tile(np.linalg.norm(g1, axis=1, ord=2, keepdims=True), (1, 3))
        ga[ia2, 0:3] = g1
        print((-np.tile(va[ia2, np.newaxis], (1, 3))).shape, ga.shape, "ga")
        direction = -np.tile(va[ia2, np.newaxis], (1, 3)) * ga[ia2, 0:3]
        # lambda_ = 0.01 * 10.0
        # print(xa.shape, "max(ia2)=",np.max(ia2))
        xa[ia2, 0:3] += direction * lambda_
        xa[ia2, 3] = 1  # not needed

        # print(len(ia2), "updated")
        print(len(ia2_accepted), "accepted")
        print("left: ", len(ia2))
        print("mean error:", np.sqrt(np.sum(va[ia]**2)))
        assert len(ia2) + len(ia2_accepted) == len(ia)

        # pass along
        ia = ia2

    return "unreachable"


def optimize_points1_inplace(iobj, xa, tolerance=0.00001, lambda_=0.1):
    """ inplace version not tested"""
    assert vectorized.is_implicit_type(iobj)
    npoints = xa.shape[0]
    assert xa.shape[1] == 4
    vectorized.check_vector4_vectorized(xa)
    ia = np.arange(0, npoints)
    va = np.ones((npoints,))
    ga = np.ones((npoints, 4))
    while True:
        va[ia] = iobj.implicitFunction(xa[ia, :])
        erra = np.abs(va[ia])
        # decide
        ia2 = ia[erra > tolerance]
        ia2_accepted = ia[erra <= tolerance]
        if ia2.size == 0:
            assert np.all(iobj.implicitFunction(xa) <= tolerance)
            return
        # update: those accepted (output): not needed.
        # update: those keep updating (loop)
        g1 = iobj.implicitGradient(xa[ia2, :])[:, 0:3]
        g1 = g1 / np.tile(np.linalg.norm(g1, axis=1, ord=2, keepdims=True), (1, 3))
        ga[ia2, 0:3] = g1
        direction = -np.tile(va[ia2, np.newaxis], (1, 3)) * ga[ia2, 0:3]
        xa[ia2, 0:3] += direction * lambda_
        xa[ia2, 3] = 1  # not needed
        assert len(ia2) + len(ia2_accepted) == len(ia)
        # pass along
        ia = ia2
    return "unreachable"


def project_points2(iobj, xa0, tolerance=0.00001, lambda_=0.1, maxdist=np.infty):
    """ using maxdist, maximum iteration count """
    MAX_ITERATION_COUNT = 30
    assert vectorized.is_implicit_type(iobj)
    xa = xa0.copy()
    # xa=xa0
    npoints = xa.shape[0]
    assert xa.shape[1] == 4
    vectorized.check_vector4_vectorized(xa)

    ia = np.arange(0, npoints)
    # xa_ = xa[ia, :] [:,0:3]
    # va[ia,:] = iobj.implicitFunction(xa[ia, :])
    # ga[ia,:] = iobj.implicitFunction(xa[ia, :])
    va = np.ones((npoints,))
    ga = np.ones((npoints, 4))

    counter = 0
    while True:
        # evaluate
        # va = iobj.implicitFunction(xa[ia, :])
        va[ia] = iobj.implicitFunction(xa[ia, :])
        erra = np.abs(va[ia])

        # decide
        assert maxdist == np.infty  # for now
        if maxdist == np.infty:
            ia2 = ia[erra > tolerance]
            ia2_accepted = ia[erra <= tolerance]
        if maxdist < np.infty:
            assert False
            how_far = np.linalg.norm(xa0[ia, :] - xa[ia, :], ord=2)
            ia2 = ia[np.logical_and(erra > tolerance, how_far <= maxdist)]
            ia2_accepted = ia[np.logical_or(erra <= tolerance, how_far > maxdist)]
            ia2_rejected = ia[np.linalg.norm(xa0[ia, :] - xa[ia, :], ord=2) > maxdist]
            xa[ia2_rejected, :] = xa0[ia2_rejected, :]  # something like that
            print("warning: not tested even once *****")

        if ia2.size == 0 or counter > MAX_ITERATION_COUNT:
            # self.unsuccessfully_projected_centroids = iobj.implicitFunction(xa) <= tolerance
            # assert np.all( iobj.implicitFunction(xa) <= tolerance )
            # todo:
            return xa  # xa*0.9 + xa0*0.1

        # update: accepted (output)
        # xa[ia2_accepted, 0:3] = xa[ia2_accepted, :]

        # update: keep updating (loop)
        #ga = iobj.implicitGradient(xa[ia2, :])  [:,0:3]
        #print(ga.shape, va.shape)
        # ga = ga [:,0:3]
        g1 = iobj.implicitGradient(xa[ia2, :])[:, 0:3]
        g1 = g1 / np.tile(np.linalg.norm(g1, axis=1, ord=2, keepdims=True), (1, 3))
        ga[ia2, 0:3] = g1
        print((-np.tile(va[ia2, np.newaxis], (1, 3))).shape, ga.shape, "ga")
        direction = -np.tile(va[ia2, np.newaxis], (1, 3)) * ga[ia2, 0:3]
        print(xa.shape, "max(ia2)=", np.max(ia2))
        xa[ia2, 0:3] += direction * lambda_
        xa[ia2, 3] = 1  # not needed

        print(len(ia2), "updated")
        print(len(ia2_accepted), "accepted")
        print("left: ", len(ia2))
        print("mean error:", np.sqrt(np.sum(va[ia]**2)))
        assert len(ia2) + len(ia2_accepted) == len(ia)

        # pass along
        ia = ia2
        counter += 1
        print(counter)

    return None     # (xa+xa0)/2.0
    # return xa


def project_points3(iobj, xa0, tolerance= 0.00001, lambda_=0.1, maxdist=np.infty):
     """ using maxdist, maximum iteration count.
     This is rubbish. """
     MAX_ITERATION_COUNT = 30

     assert vectorized.is_implicit_type(iobj)
     vectorized.check_vector4_vectorized(xa0)
     assert xa0.shape[1] == 4
     npoints = xa0.shape[0]

     xa = xa0.copy()

     #ia = indices of active points
     ia = np.arange(0, npoints)
     va = np.ones((npoints,))
     ga = np.ones((npoints, 4))
     ga_norm = np.ones((npoints,))
     ga_n = np.ones((npoints, 3))  # ga_normalized
     unsuccessfully_projected_centroids = np.ones((npoints,)) < 0  # [False]

     counter = 0
     while True:
        #evaluate
        va[ia] = iobj.implicitFunction(xa[ia, :])
        ga[ia] = iobj.implicitGradient(xa[ia, :])  [:,0:3]
        ga_norm[ia] = np.linalg.norm(ga[ia], axis=1, ord=2, keepdims=True)
        ga_n[ia] = ga[ia] / ga_norm[ia]
        #g1[ia2 of ia] = g1 / np.tile(g1_norm, (1,3) )

        erra = np.abs(va[ia] / ga_norm[ia])
        #what if some ga_norm[] is near zero?

        #decide
        if maxdist == np.infty:
            #ia2 is a logical array of size of ia
            ia2 = ia[erra > tolerance]
            ia2_dismissed = ia[erra <= tolerance]
        if maxdist < np.infty:
            how_far = np.linalg.norm(xa0[ia,:]-xa[ia, :], ord=2)
            ia2 = ia[np.logical_and(erra > tolerance, how_far <= maxdist)]
            ia2_dismissed = ia[np.logical_or(erra <= tolerance, how_far > maxdist)]
            ia2_rejected = ia[np.linalg.norm(xa0[ia,:]-xa[ia, :], ord=2) > maxdist]
            xa[ia2_rejected,:] = xa0[ia2_rejected, :]  # something like that
            print("warning: not tested even once *****")

        if ia2.size == 0   or  counter > MAX_ITERATION_COUNT:
            #self.unsuccessfully_projected_centroids = iobj.implicitFunction(xa) <= tolerance
            #assert np.all( iobj.implicitFunction(xa) <= tolerance )
            #todo:
            return xa  # xa*0.9 + xa0*0.1

        #update: accepted (output)
        #xa[ia2_dismissed, 0:3] = xa[ia2_dismissed, :]

        #update: keep updating (loop)
        #g1 = iobj.implicitGradient(xa[ia2, :])  [:,0:3]
        #g1_norm = np.linalg.norm(g1, axis=1, ord=2, keepdims=True)
        g1 = ga[ia2, :]
        g1_norm = ga_norm[ia2]
        #g1 = g1 / np.tile(g1_norm, (1,3) )
        ga_n[ia2, 0:3] = g1 / np.tile(g1_norm, (1,3) )
        ##NOT TESTED
        print((-np.tile(va[ia2, np.newaxis],(1,3))).shape, ga_n.shape, "ga")
        direction = -np.tile(va[ia2, np.newaxis],(1,3)) * ga_n[ia2, 0:3]
        #lambda_ = 0.01 * 10.0
        print(xa.shape, "max(ia2)=",np.max(ia2))
        xa[ia2, 0:3] += direction * lambda_
        xa[ia2, 3] = 1 #not needed

        print(len(ia2), "updated")
        print(len(ia2_dismissed), "accepted")
        print("left: ", len(ia2))
        print("mean error:", np.sqrt( np.sum(va[ia]**2) ))
        assert len(ia2) + len(ia2_dismissed) == len(ia)

        # pass along
        ia = ia2
        counter += 1
        print(counter)

     return None # (xa+xa0)/2.0

def project_points4(iobj, xa0, tolerance = 0.00001, lambda_=0.1, maxdist=np.infty):
    pass

#v5
def project_single_point1_mk(iobj, start_x, tolerance = 0.00001, lambda_=0.1, maxdist=np.infty):
        # initialize P = C
        #x = start_x #centroid.reshape(1,4)
        raise ImplementationError()
        tolerance = what
        counter = 0
        projection = start_x
        R = np.ones((4,1))

        if abs(iobj.function(projection)[0]) < tolerance: # if condition is already met return
            return projection
        else:
            Q = projection
        while True:

            counter += 1

            Df = self.gradient(Q)[0][:-1]
            f = self.function(Q)[0]
            Df_norm = np.linalg.norm(self.gradient(Q)[0][:-1])
            R[0][:-1] = Q[0][:-1] + self.lambda_val * (-2.0*f*Df)/(Df_norm)
            #print self.function(R)[0]
            mult = self.function(Q)[0] * self.function(R)[0]
            if  mult < 0 :
                try:
                    projection = self.find_bisection_root(Q,R)
                    break
                except:
                    print "find_bisection_root failed"
            elif mult < tolerance:
                return Q
            else:
                Q = R
        # log lambda values for every simulation so we tune it later on .
        with open('lambda_logs','a') as logfile:
            print('lambda_val, ' + repr(self.lambda_val) + ', counter, '+ repr(counter) + '\n')

        return projection


#first attempt. Not successful.
#v6
def project_single_point2_ohtake(iobj, start_x, lambda_val, max_dist ):
    #dont use this
    """ lambda_val: step size"""
    #max_iter = 20  # config["max_iter"]
    check_vector4_vectorized(start_x)
    assert start_x.shape[0] == 1

    p1 = search_near_using_gradient_ohtake(iobj, start_x, None, lambda_val, max_dist )
    f1 = iobj.implicitFunction(p1)

    p2 = 2*start_x - p1
    f2 = iobj.implicitFunction(p2)
    if f1*f2 < 0:
        direction = -(start_x - p1) # as in Ohtake
        dn = np.linalg.norm(direction)
        if dn>0.000000001:
            direction = direction/dn

        #broken
        p2_ = search_near_using_vector_ohtake(iobj, start_x, direction, lambdaa, max_dist)
        if np.linalg.norm(start_x - p2_) > np.linalg.norm(start_x - p1):
            p = p2_
        else:
            p = p1
    else:
        p = p1
    if np.linalg.norm(start_x - p) <= max_dist:
        return p
    else:
        return None



if __name__ == '__main__':
    pass
