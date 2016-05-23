#Sohail
#note: 3 points should NOT be projected into ONE point. DEFINITELY. => dont remove the faces. Only two vertices at a time.
from vtk_mc import vtk_mc
import sys
import numpy as np
from basic_types import check_vector4_vectorized
from basic_types import normalize_vector4_vectorized

from utils import Timer

from utils import optimised_used
TEST_ON = not optimised_used()


mesh_quality_settings = {
    "min_edge_len":  0.000001   # 0.01  # 0.001 microns   # 0.001, 0.007
    }

CHECK_PAIRED = True

mesh_correction = False
take_it_easy = False

B = 1000000L


class TroubledMesh(Exception):
    def __init__(self, reason):
        self.reason = reason

    def __str__(self):
        return "TroubledMesh: "+str(self.reason)
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
    #print face_triplet_ids[np.diff(face_triplet_ids)==0]
    # Check there is no repeated faces (with exact same vertex list):
    diff0 = (np.diff(face_triplet_ids) == 0)
    if np.sum(diff0) != 0:
        nonz = np.nonzero(diff0)[0]
        diff01 = diff0.copy()
        diff01[nonz+1] = True

    # diff0 versus diff01: diff01 is for print only: to print both sides (Elements) of each "diff==0"
    #nonz = np.nonzero(diff01)[0]  #both of them
    nonz = np.nonzero(diff0)[0]  # only the redundants
    bad_faces = nonz
    #but some repeated ones may remain
    #howver the original idx (of redundant faces) are = face_order[nonz]
    original_indices = face_order[bad_faces]

    #print "***", original_indices
    assert np.sum(diff0) == original_indices.size, "Number of redundant faces"
    #assert np.sum(diff0) == 0, "Repeated faces found"
    #print original_indices.shape, "original_indices.shape"
    #if original_indices.shape != (0,):
    #    assert original_indices.shape[1] == 1
    return original_indices


def check_mesh(facets):
    """ Checks if the output of the Marching Cubes and subdivision are correct. Correction of the mesh: closedness, each edge appears exactly twice, etc"""
    check_faces(facets)

    from mesh_utils import make_edge_lookup
    check_faces(facets)
    check_faces3(facets)

    # The following does a series of `assert`s which check the correctness of the mesh
    (edges_of_faces, faces_of_edges, vertpairs_of_edges) = \
        make_edge_lookup(facets)




def check_faces(faces):
    #print("------ check_faces(faces)")

    redundant_faces = check_face_triplets(faces)
    if not redundant_faces.size == 0:
        set_trace()

    assert redundant_faces.size == 0, "Repeated faces found"


    # def get_edge_code_triples_of_mesh(facets):
    #     """ Returns an array of (F)x(3), containing the 'edge codes' of sides of the faces of a mesh.
    #     There are 3 sides for each face.
    #     An 'edge code' is a long integer (int64) v1+B*v2 where v1,v2 are the indices of the ends (vertices) of the edge, where v1<v2."""
    #     e0 = facets[:, np.newaxis, [0, 1]]
    #     e1 = facets[:, np.newaxis, [1, 2]]
    #     e2 = facets[:, np.newaxis, [2, 0]]   # np view
    #     e012 = np.concatenate((e0, e1, e2), axis=1)  # n x 3 x 2
    #     assert e012.base is None  # make sure it's not a view of faces
    #     e012.sort(axis=2)
    #     BB = np.array([1L, B], dtype=np.int64)
    #     all_edges_triples = np.dot(e012, BB)  # n x 3
    #     assert all_edges_triples.dtype == np.int64
    #     assert all_edges_triples.size == 0 or np.min(all_edges_triples) >= 0
    #     assert np.max(facets, axis=None) < B
    #     assert all_edges_triples.size == 0 or np.min(all_edges_triples) >= 0
    #     assert all_edges_triples.shape == (facets.shape[0], 3)
    #     return all_edges_triples, e012


    #edg, fe3  = get_edge_code_triples_of_mesh(faces)
    """cannot refactor: fe3 is used """
    if True:
        #>begin of refactorable region
        #unique edges
        f0 = faces[:, np.newaxis, 0:2]
        f1 = faces[:, np.newaxis, 1:3]
        f2 = faces[:, np.newaxis, [0, 2]]
        f0 = f0.copy(); f0.sort(axis=2)  # changes the order in faces!
        f1 = f1.copy(); f1.sort(axis=2)
        f2 = f2.copy(); f2.sort(axis=2)

        fe3 = np.concatenate( (f0, f1, f2), axis=1 )  # shape==(:,3,2)
        BB = np.array([[1L, B]], dtype=np.int64).transpose().ravel()  # 2x-
        edg = np.dot(fe3, BB)  # fx3
        assert edg.dtype == np.int64
        assert edg.size == 0 or np.min(edg) >= 0
        assert np.max(faces, axis=None) < B
        #< end of refactorable region

    edg_1, fe3_1 = get_edge_code_triples_of_mesh(faces)
    assert np.allclose(edg_1, edg)
    assert np.allclose(fe3, fe3_1)

    #Sort edges to detect repeated edges. Each edge should appear exactly twice.
    q = edg.ravel().copy()
    sort_idx = q.argsort()
    assert sort_idx.shape == (faces.shape[0]*3,)
    q_unsorted = q.copy()
    q.sort()
    assert np.all(q_unsorted[sort_idx] == q)
    fe3_ravel = fe3.reshape( (np.prod(fe3.shape[0:2]), 2) )

    idx_xy = np.unravel_index( sort_idx[:], fe3.shape[0:2] )
    face_idx = idx_xy[0]  # tuple
    side_idx = idx_xy[1]

    vert_idx_1 = side_idx
    vert_idx_2 = (side_idx + 1) % 3
    # You can get the sorted inices here:
    v1 = faces[face_idx, vert_idx_1]
    v2 = faces[face_idx, vert_idx_2]
    v12 = np.concatenate((v1[:, np.newaxis], v2[:, np.newaxis]), axis=1)

    #Now they are the same: v12 and fe3_ravel: v12, fe3_ravel[sort_idx[:10], :]

    #Keep this:
    #assert np.all(np.diff(q)[1::2] != 0)  # Note that the very first one can be zero

    if not np.all(np.diff(q)[::2] == 0):
        set_trace()
        raise TroubledMesh("Edges are not paired")

        i1 = (np.diff(q)[::2] != 0)
        #print q[::2][i1]
        #print q[1::2][i1]
        #print np.zip(q[::2][i1], q[1::2][i1])
        i00 = np.nonzero(i1)[0][0]
        #print i00, "i1[0]"
        i0 = i00*2-4
        #for i in range(i0, i0+20): #range(zp.size):
        #    print q[i],
        #print
        #yes, some element is repeated 3 times
        if False:
            exit()


            q0 = q[::2][i1][:, np.newaxis]
            q1 = q[1::2][i1][:, np.newaxis]
            #print i1.shape
            #print q[::2].shape
            zp = np.concatenate( (q0, q1), axis=1)
            #print zp.ravel()
            #print zp.shape, "zp"


            del i1
            #for i in range(20): #range(zp.size):
            #    print zp.ravel()[i],
            #print

    if CHECK_PAIRED:
        assert np.all(np.diff(q)[::2] == 0), "Not all edges are paired"
        #print "pair3"
        #set_trace()


def visualise_edge_distribution(elist):
    ea = np.concatenate( (elist[0], elist[1], elist[2]), axis=0)
    #import matplotlib.pyplot as plt
    #plt.hist( ea , 50)
    #plt.show()
    import math
    global STEPSIZE
    special_lengths = np.array([STEPSIZE, math.sqrt(2)*STEPSIZE])
    import matplotlib.pyplot as plt
    #plt.hist( ea , 150)
    #plt.plot(special_lengths, special_lengths*0+100, "*")
    leps = 0.0000001
    plt.hist( np.log10(ea+leps), 150)
    plt.plot(np.log10(special_lengths+leps), special_lengths*0+100, "*")
    plt.show()


def visualise_distance_histogram(verts1, verts2, faces_used):
    used_verts = np.unique(faces_used.ravel())
    #print used_verts
    #set_trace()
    #import math
    global STEPSIZE
    import matplotlib.pyplot as plt

    dists = np.linalg.norm(verts1[used_verts, :] - verts2[used_verts, :], axis=1)
    plt.hist( dists, 50)
    plt.plot([STEPSIZE, STEPSIZE/2.], [2.5, 2.5], "r*")
    plt.show()


#def fix_faces_3div2(faces):
#    from mesh_utils import make_neighbour_faces_of_vertex
#    faceslist_neighbours_of_vertex = make_neighbour_faces_of_vertex(facets)
#    face_neighbours_of_faces_Fx3 = build_faces_of_faces(facets)

def REMOVE_REPEATED_EDGES(faces):
    #see check_faces()

    if True:
        # >begin of refactorable region
        f0 = faces[:, np.newaxis, 0:2]
        f1 = faces[:, np.newaxis, 1:3]
        f2 = faces[:, np.newaxis, [0, 2]]
        f0 = f0.copy(); f0.sort(axis=2)  # changes the order in faces!
        f1 = f1.copy(); f1.sort(axis=2)
        f2 = f2.copy(); f2.sort(axis=2)
        fe3 = np.concatenate( (f0, f1, f2), axis=1)  # shape==(:,3,2)

        BB = np.array([1L, B], dtype=np.int64)
        edg = np.dot(fe3, BB)   # fx3
        assert edg.dtype == np.int64
        assert edg.size == 0 or np.min(edg) >= 0
        assert np.max(faces, axis=None) < B
        # <end of refactorable region

    edg_1, fe3_1 = get_edge_code_triples_of_mesh(faces)
    assert np.allclose(edg_1, edg)
    assert np.allclose(fe3_1, fe3)

    #Sort edges to detect repeated edges. Each edge should appear exactly twice.
    #was q
    edg_sorted = edg.ravel()
    sort_idx = edg_sorted.argsort()
    assert sort_idx.shape == (faces.shape[0]*3,)
    #print sort_idx
    q_unsorted = edg_sorted.copy()
    edg_sorted.sort()
    assert np.all(q_unsorted[sort_idx] == edg_sorted)
    #print edg_sorted[:10]
    #print "diff=", np.diff(edg_sorted)[:10]

    assert np.all(np.diff(edg_sorted)[0::2] == 0), "some edges are not exactly repeated once: not manifold."
    assert np.all(np.diff(edg_sorted)[1::2] != 0), "some edges are not exactly repeated once: not manifold."
    #print np.diff(edg_sorted)
    #print edg_sorted.shape, "edg_sorted.shape", np.array(edg_sorted.shape) / 3.
    #print faces.shape
    #print (faces.shape[0])/3.
    #assert faces.shape[0] % 3 == 0
    assert edg_sorted.shape[0] % 3 == 0

    #back to pairs
    fe3_ravel = fe3.reshape( (np.prod(fe3.shape[0:2]), 2) )
    #exit()
    return

    #print fe3_ravel.shape, "ss"
    #print fe3_ravel[:5*3, :]
    #print faces[:5,:]
    #print f2
    #exit()

    #print locals().keys()

    #idx_xy = np.unravel_index( sort_idx[:10], fe3.shape[0:2] )
    idx_xy = np.unravel_index( sort_idx[:], fe3.shape[0:2] )
    #idx_xy is a tuple
    face_idx = idx_xy[0]
    side_idx = idx_xy[1]
    #print "sort_idx=", sort_idx[:10]
    #print "face, side=", face_idx, side_idx

    vert_idx_1 = side_idx
    vert_idx_2 = (side_idx + 1) % 3
    v1 = faces[face_idx, vert_idx_1]
    v2 = faces[face_idx, vert_idx_2]

    v12 = np.concatenate((v1[:, np.newaxis], v2[:, np.newaxis]), axis=1)

    #Now they are the same: v12 and fe3_ravel
    #print v12
    #print fe3_ravel[sort_idx[:10], :]

    #You can get the sorted inices here:
    #faces[face_idx, vert_idx_1]
    #faces[face_idx, vert_idx_2]

    if not np.all(np.diff(edg_sorted)[::2] == 0):
        pass






def REMOVE_REPEATED_FACES_NEVER_NECESSARY(faces):
    """ Remove faces that are exactly repeated """
    #print "ok"

    # based on check_face_triplets()
    # unique faces
    f3sides = faces.copy()
    f3sides.sort(axis=1)

    BBB = np.array([1L, B, B*B], dtype=np.int64)

    face_triplet_ids = np.dot(f3sides, BBB)
    assert face_triplet_ids.dtype == np.int64
    assert face_triplet_ids.size == 0 or np.min(face_triplet_ids) >= 0
    assert np.max(faces, axis=None) < B
    #face_triplet_ids = face_triplet_ids
    del f3sides

    #print face_triplet_ids
    #print face_triplet_ids.shape
    #exit()

    face_order = face_triplet_ids.argsort()
    face_triplet_ids=face_triplet_ids.copy()
    face_triplet_ids.sort()
    # Check there is no repeated faces (with exact same vertex list):
    diff0 = (np.diff(face_triplet_ids) == 0)
    assert np.sum(diff0) == 0
    #no repeat of faces

    if np.sum(diff0) != 0: # 74 != 0
        #get the original indices from the sorted
        redunds = np.nonzero(diff0)[0]

        # diff01: shifted version, for print only. to print both sides (Elements) of each "diff==0"
        diff01 = diff0.copy()
        diff01[redunds+1] = True
        #print face_triplet_ids[diff0]
        #print face_triplet_ids[diff01]

        #redundants = redunds + 1
        faces = np.delete( faces, redunds + 1, axis=0 )
        print "Redunadnt faces", redunds + 1
        exit()
    #print "no problem"
    #exit()
    return faces


def unused_vertices_slow(verts, facets):
    killed_whichvertices = np.zeros((verts.shape[0],), dtype=np.bool) + True
    killed_whichvertices[facets.ravel()] = False
    idx = np.nonzero(killed_whichvertices)[0]
    idx.sort()
    return idx


def map_vertices_of_nil_faces(faces, nil_areas_whichfaces):
    """ Combines the vertices for removed faces, instead of removing the face altogether. """

    faces_to_kill = faces[nil_areas_whichfaces].copy()
    number_of_faces_to_kill = faces_to_kill.shape[0]
    #print "going to kill", number_of_faces_to_kill, "degenerate faces"

    faces_to_kill.sort(axis=1)
    face_map10 = np.zeros((2, 0), dtype=int)
    for ei in [1, 2]:
        from_idx = faces_to_kill[:, ei]
        to_idx = faces_to_kill[:, 0]
        #print to_idx, to_idx.shape, "*"
        assert np.ndim(from_idx) == 1
        assert np.ndim(to_idx) == 1
        temp_map10 = np.concatenate( (to_idx[np.newaxis, :], from_idx[np.newaxis, :]), axis=0)
        #print to_idx[np.newaxis, :].shape, "--"
        #print temp_map10.shape
        assert temp_map10.shape[0] == 2
        #print np.ndim(temp_map10)
        assert np.ndim(temp_map10) == 2

        face_map10 = np.concatenate( (face_map10, temp_map10), axis=1)
        assert face_map10.shape[0] == 2
        assert np.ndim(face_map10) == 2

    #print faces.shape

    #print face_map10.shape, "******************"  # too many
    #print "???", faces.shape, faces.shape[0]*2  #WRONG  ****************************************

    fmap_ = np.arange(faces.shape[0])
    fmap_[face_map10[1, :]] = face_map10[0, :]

    for i in range(6):
        fmap_ = fmap_[fmap_]
        dosnt_need_change = np.all(fmap_ == fmap_[fmap_])  # 4 times False!
        #print dosnt_need_change
    assert dosnt_need_change

    #print fmap_
    #print fmap_[fmap_]
    assert np.all(fmap_ == fmap_[fmap_])
    #print "number of mapped vertices", np.sum(fmap_ != np.arange(faces.shape[0]))
    dying_vertices = np.arange(faces.shape[0])[fmap_ != np.arange(faces.shape[0])]
    dying_vertices_bool = fmap_ != np.arange(faces.shape[0])
    affected_faces = np.any(dying_vertices_bool[faces], axis=1 )
    #print "affected faces", np.sum(affected_faces)
    #print "full faces", np.sum( np.all(dying_vertices_bool[faces], axis=1 ) )
    #number of faces that all of their elements are mapped onto something

    #print "2 faces", np.sum( np.sum(dying_vertices_bool[faces], axis=1 )>=2 ), "faces that have 2 vertices"
    for s in range(5):
        #print "faces with %d taken vertices:"%(s,), np.sum( np.sum(dying_vertices_bool[faces], axis=1 )==s )
        pass
    #print np.nonzero(affected_faces)

    face__projected_vertices = fmap_[faces]

    v10 = face__projected_vertices[:, 1] == face__projected_vertices[:, 0]
    v20 = face__projected_vertices[:, 2] == face__projected_vertices[:, 0]
    v21 = face__projected_vertices[:, 1] == face__projected_vertices[:, 2]
    v1020 = np.logical_and(v10, v20)
    #print np.sum(v1020), "faces with three (now) identical vertices"

    # Why is it 81 and 80 = number_of_faces_to_kill. (fixme)
    #print faces[v1020, :]
    #p0 = 131  # to what index is p0 projected?
    #print p0, "->",face_map10[0, :][face_map10[1, :]==p0]
    #exit()

    faces_to_annihilate = np.logical_or(v10, v20, v21)  #129 ?!
    #print "faces_to_annihilate", np.sum(faces_to_annihilate)

    #return fmap_
    #return faces[fmap_]
    #return fmap_[faces]
    return face__projected_vertices, faces_to_annihilate


# BUG IS HERE: verts, new_verts, etc
def delete_unused_vertices(verts, faces):
    assert verts.shape[0] > 0
    assert faces.shape[0] > 0
    return verts, faces
    assert np.isreal(verts[0, 0])
    killed_whichvertices = np.zeros((verts.shape[0],), dtype=np.bool) + True
    killed_whichvertices[faces.ravel()] = False
    #assert ***
    not_killed = np.logical_not(killed_whichvertices)
    i123 = np.arange(verts.shape[0])  # np.arange(faces.shape[0])
    verts_newidx = i123[not_killed]
    i123[np.logical_not(killed_whichvertices)] = np.arange(np.sum(not_killed))
    i123[killed_whichvertices] = -1
    #i123: index=old v index  value=new v index
    new_faces = i123[faces]
    assert not np.any(new_faces.ravel() == -1)
    new_verts = verts[verts_newidx]

    assert np.all(new_faces.ravel() >= 0)
    assert np.all(new_faces.ravel() < new_verts.shape[0] )
    assert new_faces.shape[1] == 3
    assert np.ndim(new_faces) == 2

    assert new_verts.shape[1] == 3
    assert np.ndim(new_verts) == 2

    assert issubclass(faces.dtype.type, np.integer)
    assert new_verts.shape[0] > 0
    assert np.isreal(new_verts[0, 0])

    return new_verts, new_faces


#todo: review this code: , separate it, u-test it.
def remove_vertices_and_faces(verts, faces, nil_areas_whichfaces, map12):
    """ nil_areas_whichfaces: indices of faces that need to be removed.
    map12: the vertices that have to be replaced (remove map12[:,1] -> replace as: map12[:,0])."""

    assert faces.shape[0] > 0
    assert verts.shape[0] > 0

    assert map12.shape[0] == 2
    #print "-"
    u =unused_vertices_slow(verts, faces)
    #print "unused vertices: ", u
    #for x in u:
    #    print np.sum(faces.ravel()==x),
    #print ""

    if False:
        #Note: there are cases where the edges are not equal (0.0000001) but the area is small (0.001).
        #print faces[181, :], faces[219, :], faces[952, :] # 181  219  952  993 1281 1295
        f1 = faces[ np.array([181,219,952]), :]
        #print "areas = ", compute_triangle_areas( verts, f1 , return_normals=False, AREA_DEGENERACY_THRESHOLD=None)
        #print verts[f1,:]

    if False:
        #(1) This step can be done either after or before the step (2)
        new_faces = faces[np.logical_not(nil_areas_whichfaces)]
    else:
        new_faces = faces
    #todo: not all edges will be paired

    assert new_faces.shape[0] > 0


    #(2) Make a map_ to change the bad vertices to good ones, even if a -> b -> c.
    # Since this does not change the index of faces, the actual removing of the degenerate faces can be done either after of before this step.
    map_ = np.arange(verts.shape[0])
    map_[map12[1, :]] = map12[0, :]

    for i in range(2):
        map_ = map_[map_]
    assert np.all(map_ == map_[map_])

    #print "*new_faces", new_faces

    if True:
        #for i in range(3):
        #    new_faces = map_[new_faces]
        new_faces = map_[new_faces]
        assert np.all(new_faces == map_[new_faces])

        #assert no map12[1, :] exist in faces

        #wrong test:
        #assert np.intersect1d(map12[1, :], new_faces.ravel()).size == 0, "Confirm no more vertex clustering necessary."
        from_ = np.nonzero(map_ != np.arange(map_.shape[0]))[0]
        assert type(from_) is np.ndarray
        assert np.ndim(from_) == 1
        #print np.intersect1d(from_, new_faces.ravel())
        assert np.intersect1d(from_, new_faces.ravel()).size == 0, "Confirm no more vertex clustering necessary, hence the map is complete, does not need more repeats, and is not cylic.."


    assert new_faces.shape[0] > 0

    #print "map_", map_
    #print "new_faces", new_faces
    new_verts, new_faces = delete_unused_vertices(verts, new_faces)

    assert new_faces.shape[0] > 0

    #killed_whichvertices = np.zeros((verts.shape[0],), dtype=np.bool) + True
    #killed_whichvertices[new_faces.ravel()] = False
    #print "unused count", np.sum(killed_whichvertices)
    #print "unused vertices: ",  unused_vertices_slow(verts, new_faces)
    #new_verts = verts
    unused_vertices = unused_vertices_slow(verts, new_faces)


    new_verts, new_faces = delete_unused_vertices(verts, new_faces)

    assert new_faces.shape[0] > 0
    assert new_verts.shape[0] > 0

    """
    new_faces = faces[np.logical_not(nil_areas_whichfaces)]
    #faces_newidx = np.arange(faces.shape[0])[np.logical_not(nil_areas_whichfaces)]
    #index = new iddex, value = old index

    #bad_vertices = map12[1, :]
    killed_whichvertices = np.zeros((verts.shape[0],), dtype=np.bool)
    killed_whichvertices[map12[1, :]] = True
    #verts_newidx = np.arange(faces.shape[0])[np.logical_not(killed_whichvertices)]
    #index = new index, value = old index
    #new_verts = verts[np.logical_not(killed_whichvertices), :]
    #new_faces = new_faces[verts_newidx]  # ##????


    #todo: write a test
    ***not tested
    return new_faces, new_verts
    """

    #print new_faces[892, :], "892<<<<<<<<<<<<<<<<<<<<<<<"  # was [396 397 476]  -> now: 290,290,290
    #combine the vertices for removed faces, instead of removing the face altogether.

    #faces, faces_to_annihilate = map_vertices_of_nil_faces(faces)[faces]
    new_faces2, faces_to_annihilate = map_vertices_of_nil_faces(new_faces, nil_areas_whichfaces)
    #print new_faces2  #most of them are zero

    #print "n"*100
    #print new_faces.shape,
    #print new_faces2.shape
    ###########################################################
    if False:
        quick_vis(verts, new_faces, nil_areas_whichfaces)
    ####################################

    assert new_faces2.shape[0] > 0
    assert new_faces2.shape[0] > 0


    #checking if the nil_areas_whichfaces matches the faces that are now projected
    eq1 = new_faces2[:, 0] == new_faces2[:, 1]
    eq2 = new_faces2[:, 0] == new_faces2[:, 2]
    eq3 = new_faces2[:, 1] == new_faces2[:, 2]

    eq1and2and3 = np.logical_and(eq1, eq2)  # The triangles which all vertices have the same index number.
    #eq1or2 = np.logical_or(eq1, eq2)
    #set_trace()

    #eq1or2or3 = np.logical_or(eq1, eq2, eq3) # wrong!
    eq1or2or3 = np.logical_or(eq1, np.logical_or(eq2, eq3))  # make it more efficient

    #iii=402
    #print eq1or2or3[iii], eq1[iii], eq2[iii], eq3[iii]
    #print new_faces2[iii, 1] == new_faces2[iii, 2]
    #print  new_faces2[iii, :]
    #exit()

    #print eq1and2and3, " sum1=",np.sum(eq1and2and3)  # 81
    #print " sum1or2=",np.sum(eq1or2)
    #print nil_areas_whichfaces, " sum2=", np.sum(nil_areas_whichfaces)  # 80

    #81 and 80 consequently
    #print np.nonzero( nil_areas_whichfaces != eq1and2and3)

    #print "XOR", np.sum( eq1and2and3 != nil_areas_whichfaces)
    xor_ = eq1and2and3 != nil_areas_whichfaces
    #print new_faces2[xor_, :], np.nonzero(xor_)   # [290,290,290], 892

    #print np.vstack( (eq1or2or3,nil_areas_whichfaces))

    #print np.sum( np.logical_and(np.logical_not(eq1and2and3), nil_areas_whichfaces))
    #The triangles with zero area should be a subset of eq1and2and3,ones that have a repeated edge.
    #print new_faces2[nil_areas_whichfaces, :]
    verts_123_equal_same_as_zero_area = np.sum( np.logical_and(np.logical_not(eq1and2and3), nil_areas_whichfaces)) == 0, "assert nil_areas_whichfaces is-subset eq1and2and3"
    assert verts_123_equal_same_as_zero_area

    a = np.logical_and(np.logical_not(eq1or2or3), nil_areas_whichfaces)
    #print "-"*10
    iii= np.nonzero(a)[0]
    #print iii
    #print eq1or2or3[iii], eq1[iii], eq2[iii], eq3[iii], new_faces2[iii, 1] == new_faces2[iii, 2]
    #print new_faces2[a,:]

    assert np.sum( np.logical_and(np.logical_not(eq1or2or3), nil_areas_whichfaces)) == 0, "assert nil_areas_whichfaces is-subset eq1or2or3"

    #assert np.sum( np.logical_and(np.logical_not(eq1and2and3), nil_areas_whichfaces)) == 0, "assert nil_areas_whichfaces is-subset eq1and2and3"
    #TODO: CHECK HERE
    #todo: also refactor the other code

    def x():
        n = np.logical_not; a = eq1and2and3; b = nil_areas_whichfaces; aNd = np.logical_and
        #print "a & b", np.sum(aNd(a, b))  # 80
        #print "a & ~b", np.sum(aNd(a, n(b)))  # 1   Faces that are not (area=0), but their vertices are not identical. But how?
        #print "~a & b", np.sum(aNd(n(a), b)) # 0
    x()

    a_n_b = np.nonzero(np.sum(np.logical_and(eq1and2and3, np.logical_not(nil_areas_whichfaces))))[0]
    #print a_n_b
    #quick_vis(verts, faces, a_n_b)


    #print new_faces2[eq1and2and3, :]
    #print new_faces2[nil_areas_whichfaces, :] #nice
    #81 and 80 consequently


    #print "DELETING: %d out of %d"%(-np.sum(eq1and2and3), new_faces2.shape[0]), "="*90
    #import time; time.sleep(2)

    #Now: faces with repeats in their vertices should be cut away.
    #It is also the concern of: map_vertices_of_nil_faces
    #new_faces = new_faces2
    #new_faces = new_faces2[ eq1and2and3 ]  # kill! -> oops.  wrong ones.
    if False:
        new_faces = new_faces2[ np.logical_not(eq1and2and3) ]  # kill!

    new_faces = new_faces2[ np.logical_not(eq1or2or3) ]

    assert new_faces.shape[0] > 0

    #Now there are repeated faces.

    redundant_faces = check_face_triplets(new_faces)

    #fails:
    # assert redundant_faces.size==0, "Repeated faces found***"
    #i1n = np.zeros( (new_faces.shape[0],), dtype=bool) + True
    #i1n[redundant_faces] = False
    #print new_faces
    #print redundant_faces

    assert new_faces.shape[0] > 0

    #print "DELETING: %d out of %d"%(redundant_faces.size, new_faces.shape[0]), "="*90
    #import time; time.sleep(2)

    new_faces = np.delete(new_faces, redundant_faces, axis=0)

    assert new_faces.shape[0] > 0

    #todo: D.R.Yself.
    eq1 = new_faces[:, 0] == new_faces[:, 1]
    eq2 = new_faces[:, 0] == new_faces[:, 2]
    eq1and2and3 = np.logical_and(eq1, eq2)

    #print "DELETING: %d out of %d"%(np.sum(eq1and2and3), new_faces.shape[0]), "-"*290
    #print new_faces[eq1and2and3, :]
    #import time; time.sleep(2)

    new_faces = new_faces[np.logical_not(eq1and2and3), :]

    assert new_faces.shape[0] > 0

    #print new_faces
    #print "-"*200
    check_faces(new_faces)
    #passed!

    return new_verts, new_faces



global still_fixed
still_fixed = False
from ipdb import set_trace
# **********************************
def check_degenerate_faces(verts, facets, fix_mode="dontfix"):
    """ This functions checks for Mesh problems, and fied them. It can run in three modes: quit if error, return False if error (dontfix), fix it if error (fix). In the 'fix' mode, the return is verts, faces. Otherwise, the return is boolean. The assert mode returns nothing."""
    # todo: also check facets in-itself.
    #return verts, facets, False
    #set_trace()

    fix_them = (fix_mode == "fix")
    if_assert = (fix_mode == "assert")
    if fix_mode == "fix":
        pass
    elif fix_mode == "assert":
        pass
    elif fix_mode == "dontfix":
        pass
    else:
        raise InvalidUsage()

    #not
    if fix_mode == "assert":
        #if not take_it_easy:
            check_faces(facets)

    #any_correction = False
    e1 = np.linalg.norm(verts[facets[:, 0], :] - verts[facets[:, 2], :], axis=1)
    e2 = np.linalg.norm(verts[facets[:, 1], :] - verts[facets[:, 0], :], axis=1)
    e3 = np.linalg.norm(verts[facets[:, 2], :] - verts[facets[:, 1], :], axis=1)

    assert not np.any(np.isnan(verts.ravel()))
    #tag: pairs of points to combine
    #triple of points to combine.. .
    #points with nan.  => error
    #
    #triangles_need_deletion_because_of_zero_edges = False

    nillfaces_zeroedge = np.zeros((facets.shape[0],), dtype=np.bool)
    #faces_to_remove

    v12 = np.zeros((2, 0), dtype=np.int64)

    any_zero_edge = False
    el = [e1, e2, e3]
    any_nan = False
    # visualise_edge_distribution(el)  # keep this line
    for i in [0, 1, 2]:
            j1 = i; j2 = (i-1+3) % 3
            tv3 = verts[facets[:, j1], :] - verts[facets[:, j2], :]
            ei = np.linalg.norm(tv3, axis=1)
            ineq = ei < mesh_quality_settings["min_edge_len"]
            if np.any(ineq):
                idx = np.nonzero(ineq)  # type: Tuple[np.ndarray]
                #set_trace()
                assert len(idx) == 1
                idx = idx[0]
                assert np.ndim(idx) == 1
                assert np.ndim(ei) == 1
                v1 = facets[idx, j1]
                v2 = facets[idx, j2]
                v12_ = np.concatenate((v1[np.newaxis, :], v2[np.newaxis, :]), axis=0)
                v12 = np.concatenate( (v12, v12_), axis=1)
                # vertices_to_combine

                nillfaces_zeroedge[idx] = True
                any_zero_edge = True
                del idx

            inan = np.isnan(el[i])
            assert not np.any(inan)
            #if np.any(inan):
            #    assert False
            #    idx = np.nonzero(inan)
            #    assert len(idx) == 1
            #    idx = idx[0]
            #    assert np.ndim(idx) == 1
            #    assert np.ndim(ei) == 1
            #    #print idx.size, "nan edges"
            #    any_zero_edge = True
            #    any_nan = True
            #    del idx
            #del inan

    #if still_fixed:
    #    set_trace()
    #    pass
    assert not any_nan,  "not sure"
    v12.sort(axis=0)  # make edges unique
    sv12 = v12  #sv12 is list of all edges that are small. If length zero => good.
    if sv12.size == 0:
        triangles_need_deletion_because_of_zero_edges = False
        pass
    else:
        triangles_need_deletion_because_of_zero_edges = True
        eaa = np.dot(np.array([1L, B], dtype=np.int64), sv12)
        assert eaa.dtype == np.int64
        assert eaa.size == 0 or np.min(eaa) >= 0
        assert np.max(sv12, axis=None) < B

        sort_order = eaa.argsort()
        #eaa_s = eaa.copy(); eaa_s.sort(); #print eaa_s

        # sort sv12 based on their edge codes (eaa)
        sv12 = sv12[:, sort_order]
        if CHECK_PAIRED:
            assert np.all((np.diff(sv12, axis=1)[:, ::2]).ravel() == 0)  # Make sure each edge is repeated exactly once.
        map12 = sv12[:, 0::2]  # project map12[1,:] into map12[0,:]
        edge_vects = verts[sv12[0, :]] - verts[sv12[1, :]]
        # Assert that it's always almost-zeros.
        assert np.allclose(edge_vects, 0, mesh_quality_settings["min_edge_len"])
        # Vector subtractionis sligthly more tight than norm but we use the same tolerance here (more conservatirve).
        # Since the actual different is often actually zero, this should hold. Correct the tolerance if this assert failed.

    average_edge = (np.mean(e1)+np.mean(e2)+np.mean(e3))/3.
    # note: the average edge length may be slightely less than this after removing the repeated vertices
    del average_edge

    if False:
        # Do the projection using map12
        # project map12[1,:] into map12[0,:]
        #faces[ faces == map12[1,:] ] =
        temp = map12[1, :].copy()
        temp.sort()
        # vertices to be removed: map12.shape[1]
        for ei in range(map12.shape[1]):
            e1 = map12[0, ei]
            e2 = map12[1, ei]
            assert np.isscalar(e1)
            assert np.isscalar(e2)
            assert not e1 == e2
            facets[ facets == e2 ] = e1

        #make sure bad vertices are erased.
        bad_vertices = map12[1, :]
        #assert bad_vertices dont exist in faces.ravel()

    facet_areas = compute_triangle_areas(verts, facets, return_normals=False, AREA_DEGENERACY_THRESHOLD=None)
    assert not np.any(np.isnan(facet_areas))

    assert np.all(facet_areas >= 0)
    #AREA_DEGENERACY_THRESHOLD = 0.00001  #  == 0.003 **2
    #AREA_DEGENERACY_THRESHOLD = -1.  # 0.00001 ** 2
    AREA_DEGENERACY_THRESHOLD = 0.000000000001

    ineq = facet_areas < AREA_DEGENERACY_THRESHOLD
    degenerates_count = len(facet_areas[ineq])
    #print "degenerates_count", degenerates_count

    degenerate_faces = ineq
    assert np.ndim(ineq) == 1
    zeroarea_face_indices = np.nonzero(ineq)[0]
    assert np.ndim(zeroarea_face_indices) == 1
    if len(zeroarea_face_indices) == 0:
        # No degenerate area found
        pass
    any_degenerate_area = len(zeroarea_face_indices) > 0

    #
    #obliterate_because_of_midpoint
    obm_faces = []
    obm_map = {}

    #type3_whichfaces = np.zeros((zeroarea_face_indices.size,), dtype=np.bool)
    type3_whichfaces = np.zeros((facets.shape[0],), dtype=np.bool)

    zeroarea_whichfaces = np.zeros((facets.shape[0],), dtype=np.bool)  # also see nillfaces_zeroedge
    zeroarea_whichfaces[zeroarea_face_indices] = True
    #del zeroarea_face_indices  # todo: rename zeroarea_face_indices
    #nf = facets.shape[0]
    for fi in zeroarea_face_indices:  # last use of zeroarea_face_indices
        assert degenerate_faces[fi]
        if degenerate_faces[fi]:
            #print("face:", fi, facets[fi,:])  # ('face:', 181, array([131,  71, 132]))
            degen_triangle = verts[facets[fi, :], :]  # numverts x 3
            v1 = degen_triangle[1, :] - degen_triangle[0, :]
            v2 = degen_triangle[2, :] - degen_triangle[0, :]
            v3 = degen_triangle[2, :] - degen_triangle[1, :]

            #There are Three possibilities for a mesh with zero area:
            # "OBM" case

            #choose the longest side
            v123 = np.vstack( (degen_triangle[1, :] - degen_triangle[0, :],
                degen_triangle[2, :] - degen_triangle[0, :],
                degen_triangle[2, :] - degen_triangle[1, :])
            )
            #print v123[side_index,:]
            #print v123
            side_lengths = np.linalg.norm(v123, axis=0)
            longest_side = np.argmax(side_lengths)
            #print v123[longest_side, :]
            other_vertex = (longest_side+2) % 3 # The vertex that is betwen the other two. The midpoint.
            #other_sides1 = (longest_side+2+1) % 3
            #other_sides2 = (longest_side+2-1) % 3
            ends_vertices = (longest_side, (longest_side+1) % 3)  # WRONG!    #*** FIXME!
            M = degen_triangle[other_vertex, :]
            v1M = degen_triangle[ends_vertices[0], :] - M
            v2M = degen_triangle[ends_vertices[1], :] - M
            cross = np.linalg.norm(np.cross(v1M, v2M))

            #print side_lengths[longest_side], cross  # ~ 2e-16
            type3_whichfaces[fi] = cross

            a = v1M
            b = v2M

            case = "-"
            if np.linalg.norm(a) < 0.00000001 and np.linalg.norm(b) < 0.00000001:
                case = "000"
            else:
                case = "OBM"

            if case == "OBM":
                assert np.linalg.norm(dist) < 0.00000001
                assert  np.linalg.norm(a-b) > 0.000001  #and: two sides: a-b are far    # *** FIXME!
                assert np.linalg.norm(a-b) > 0.000001
                assert np.linalg.norm(a) > 0.000001  # M is at zero
                assert np.linalg.norm(b) > 0.000001

                #proj = v1M - (v1M)
                proj = a - (a-b)*np.dot(a, a-b)/np.dot(a-b, a-b)
                #proj is exactly on the line
                cross_proj = np.linalg.norm(np.cross(v1M-proj, v2M-proj))
                print a,b
                print cross_proj, cross
                assert cross_proj <= cross
                dist = proj-0  # M is already subtracted from proj
                #print "dist", dist, np.dot(a-b, v1M), np.dot(a-b, v2M)
                #print "-"


                if np.linalg.norm(dist) < 0.00000001 and np.linalg.norm(a-b) > 0.000001:  #and: two sides: a-b are far    # *** FIXME!
                    #_type=3 (midpoint)

                    assert np.linalg.norm(a-b) > 0.000001
                    assert np.linalg.norm(a) > 0.000001  # M is at zero
                    assert np.linalg.norm(b) > 0.000001


                    #facets[fi, ends_vertices]

                    def calculate_edge_code(e):
                        assert np.ndim(e) == 2
                        assert e.shape[1] == 2

                        BB = np.array([1L, B], dtype=np.int64)
                        edge_codes = np.dot(e, BB)  # n x ,
                        assert edge_codes.dtype == np.int64
                        assert edge_codes.size == 0 or np.min(edge_codes) >= 0
                        assert np.max(facets, axis=None) < B
                        assert edge_codes.size == 0 or np.min(edge_codes) >= 0
                        return edge_codes

                    (ev1, ev2) = facets[fi, ends_vertices[0]], facets[fi, ends_vertices[1]]
                    edge_code = calculate_edge_code(np.array([[ev1,ev2]]))[0]
                    #sidesabc = sorted(list((np.linalg.norm(b), np.linalg.norm(a), np.linalg.norm(a-b))))
                    sidesabc = sorted(list((np.linalg.norm(a), np.linalg.norm(b), np.linalg.norm(a-b))))
                    #print fi, sidesabc, sidesabc[2]-(sidesabc[1]+sidesabc[0]),
                    #print "\t", edge_code, "->", facets[fi, other_vertex]

                    obm_map[edge_code] = facets[fi, other_vertex]

                    obm_faces += [fi]
                    #set_trace()
                    #print ""

                else:
                    #if not on that line but very small =>
                    assert False


    #set_trace()
    #print obm_map
    #print obm_faces
    if len(obm_faces) > 0:
        #set_trace()
        pass

    #"both area and edge are zero"
    lboth = np.nonzero( np.logical_and(zeroarea_whichfaces, nillfaces_zeroedge) )[0]
    #"zero-area only". Not empty because sometimes area is empty but no edge is zero.
    lb = np.nonzero( np.logical_and(zeroarea_whichfaces, np.logical_not(nillfaces_zeroedge)) )[0]
    #la is not informative
    la = np.nonzero( np.logical_and(np.logical_not(zeroarea_whichfaces), nillfaces_zeroedge) )[0]
    #DOES NOT HOLD:
    # assert len(la) == 0  # la, "zero-edge only", should be empty. Because 'zero edge' => 'zero area'
    #la = 0_E - 0_A
    assert len(la) == 0

    #la,lc = 0, lb==obm_which_faces
    obm_which_faces = np.zeros((facets.shape[0],), dtype=np.bool)
    obm_which_faces[obm_faces] = True
    lc = np.nonzero( np.logical_and(np.logical_not(zeroarea_whichfaces), obm_which_faces) )[0]
    assert len(lc) == 0

    #area_or_edges = np.logical_or(zeroarea_whichfaces, nillfaces_zeroedge)
    #lor = np.nonzero( area_or_edges )[0]

    #keep the following comments:
    #print "zero-area only", lb
    #print "both area and edge are zero", lboth

    #nillfaces_zeroedge is-subset-of zeroarea_whichfaces

    #vertices_to_combine = lboth
    #faces_to_remove = lb

    if False:
        print any_zero_edge, any_degenerate_area
        print la
        print lb
        print lboth
        print facets[lboth, :]
        v012 = verts[facets[lboth, :], :]
            #8x3x3
        print np.diff(v012, axis=1)
        exit()

    if not take_it_easy:
        check_faces(facets)

    FIX_OBM = True # True # False
    #print "*********"*100
    DID =0
    if fix_them:
        assert not triangles_need_deletion_because_of_zero_edges
        if triangles_need_deletion_because_of_zero_edges:

            #set_trace()
            assert map12.shape[0] == 0
            #***lor
            verts, facets = remove_vertices_and_faces(verts, facets, zeroarea_whichfaces, map12)

            set_trace()
            check_faces(facets)
            #print zeroarea_whichfaces
            #print "ZERO AREA",
            DID = 1

        if not take_it_easy:
            check_faces(facets)

        if FIX_OBM and len(obm_faces) > 0:
            assert not triangles_need_deletion_because_of_zero_edges, "If True, we cannot remove both. We need to OR first."
            #print "going to delete ", facets.shape
            #facets = np.delete(facets, obm_faces, axis=0)
            #print "deleted ", facets.shape
            DID = 2
            #check_faces(facets)

    print
    #print "*"*100
    #print "DID", DID

    #check_faces(facets)

    if fix_them:
        if FIX_OBM and len(obm_map) > 0:
            edge_array = np.array(obm_map.keys(), dtype=np.long)

            #set_trace()
            #print "goting to subdiv sides ", facets.shape, edge_array.shape
            old_facets = facets.copy()

            facets = subdivide_1to2_multiple_facets(facets, edge_array, obm_map, careful_for_twosides=False)
            #print "did subdiv sides ", facets.shape, edge_array.shape

            #print "goting to subdiv sides ", facets.shape, edge_array.shape
            facets = subdivide_1to2_multiple_facets(facets, edge_array, obm_map, careful_for_twosides=False)
            #print "did subdiv sides ", facets.shape, edge_array.shape
            #report removed/subdivided


            #print "goting to subdiv sides ", facets.shape, edge_array.shape
            facets = subdivide_1to2_multiple_facets(facets, edge_array, obm_map, careful_for_twosides=False)
            #print "did subdiv sides ", facets.shape, edge_array.shape


            e3, dummy = get_edge_code_triples_of_mesh(facets)  # f x 3
            e01 = e3[:, 0] == e3[:, 1]
            e02 = e3[:, 0] == e3[:, 2]
            e12 = e3[:, 1] == e3[:, 2]
            tr_idx = np.nonzero(np.logical_or(np.logical_or(e01, e02),e12))[0]
            if tr_idx.size > 0:
                facets = np.delete(facets, tr_idx, axis=0)
                #print facets
                #set_trace()

            totally_redundant_faces = check_face_triplets(facets)
            #d = face_triplet_codes(faces)
            if  totally_redundant_faces.size > 0:
                facets = np.delete(facets,  totally_redundant_faces, axis=0)

            #set_trace()
            #set_trace()
            e3, dummy = get_edge_code_triples_of_mesh(facets)
            er = e3.ravel()
            #np.intersect1d(er, edge_array)
            bl = np.lib.arraysetops.in1d(er, edge_array)
            #bl3=bl.reshape(e3.shape)

            #fails because the edges are not paired
            #not yet
            if False:
                check_faces(facets)

    #fails
    #May not be fine yet
    #if True:
    #    check_faces(facets)
    if fix_mode == "assert":
        check_faces(facets)
    #***not tested
    #todo: write unit-test
    #exit()

    any_correction = any_zero_edge or any_degenerate_area

    #print any_zero_edge , any_degenerate_area
    #assert fix_them != if_assert, (fix_them, if_assert)
    if if_assert:
        if any_correction:
            raise AssertionError
        #print "fine"
        return

    if fix_them:
        return verts, facets, any_correction
    else:
        return any_correction


#def compute_triangle_areas(verts, faces, return_normals=False, AREA_DEGENERACY_THRESHOLD=0.00001):
def compute_triangle_areas(verts, faces, return_normals=False, AREA_DEGENERACY_THRESHOLD=None):
    """ facet_normals: can contain NaN if the area is zero.
    If AREA_DEGENERACY_THRESHOLD is None or negative, the NaN is not assiged in output """
    # see mesh1.py ::     def calculate_face_areas(self)
    nfaces = faces.shape[0]
    expand = verts[faces, :]
    assert expand.shape == (nfaces, 3, 3)
    assert expand[:, 2, :].shape == (nfaces, 3)
    a = np.cross(
        expand[:, 1, :] - expand[:, 0, :],
        expand[:, 2, :] - expand[:, 0, :],
        axis=1)
    facet_areas = np.linalg.norm(a, axis=1, ord=2) / 2.0
    #fixme: It is unnecessary to set them to NaN
    if AREA_DEGENERACY_THRESHOLD is not None:
        degenerates_count = len(facet_areas[facet_areas < AREA_DEGENERACY_THRESHOLD])
        facet_areas[facet_areas < AREA_DEGENERACY_THRESHOLD] = np.nan  # -1
        if degenerates_count > 0:
            pass

    if not return_normals:
        return facet_areas
    else:
        #print facet_areas.shape
        assert facet_areas[:, np.newaxis].shape == (nfaces, 1)
        at_autobroadcast = facet_areas[:, np.newaxis] * 2.0
        facet_normals2 = a / at_autobroadcast
        if TEST_ON:
            at = np.tile(facet_areas[:, np.newaxis], (1, 3)) * 2.0
            facet_normals = a / at
            assert np.allclose(facet_normals, facet_normals2, equal_nan=True)
        #facet_normals2[np.isnan(facet_areas), :] = 0.
        return facet_areas, facet_normals2


def compute_centroids(verts, facets):
    # see Mesh1::build_centroids
    #self.calculate_face_areas()
    expand = verts[facets, :]
    nfacets = facets.shape[0]
    assert expand.shape == (nfacets, 3, 3)
    assert np.allclose(verts[facets[:],:], expand)
    centroids3 = np.mean( verts[facets[:],:], axis=1)  # again
    centroids = np.concatenate( (centroids3, np.ones((nfacets,1))), axis=1)
    return centroids

#evaluate_centroid_gradients
def compute_centroid_gradients(centroids, iobj, normalise=True):
    # see mesh1.py :: evaluate_centroid_gradients
    assert centroids is not None
    check_vector4_vectorized(centroids)
    centroid_gradients = iobj.implicitGradient(centroids)
    assert not np.any(np.isnan(centroid_gradients))
    assert not np.any(np.isinf(centroid_gradients))
    if normalise:
        centroid_normals = normalize_vector4_vectorized(centroid_gradients)
        return centroid_normals
    else:
        return centroid_gradients

def check_faces3(facets):
    #REMOVE_REPEATED_FACES_NEVER_NECESSARY(facets)
    REMOVE_REPEATED_EDGES(facets)
    pass


def build_faces_of_faces(facets):
    """ builds lookup tables. The result if an array of nfaces x 3,
    containing the face index of neighbours of each face.
    Since each face has exactly three neighbours, the size of the result is n x 3."""
    from mesh_utils import make_edge_lookup
    #pudb.set_trace()
    #set_trace()
    #set_trace()
    check_faces(facets)
    check_faces3(facets)
    (edges_of_faces, faces_of_edges, vertpairs_of_edges) = \
        make_edge_lookup(facets)
        # ****

    # need: face_neighbours_of_faces_Fx3
    # e__nf_x_3 = edges_of_faces[facets]
    # print e__nf_x_3.shape
    # assert e__nf_x_3.shape == (nfaces, 3, 3)
    #print edges_of_faces.shape
    nfaces = facets.shape[0]
    assert edges_of_faces.shape == (nfaces, 3)
    f1 = faces_of_edges[edges_of_faces, 0]  # first face of all edges of all faces : nf x 3 -> nf
    f2 = faces_of_edges[edges_of_faces, 1]  # second face of all edges of all faces: nf x 3 -> nf   #3 for 3 sides (edges)
    # one of the two (0:2) is equal to index=face. [index,3, 0:2 ]
    f12 = faces_of_edges[edges_of_faces, :]
    #print f12.shape
    assert f12.shape == (nfaces, 3, 2)
    #print f12[0:5, :, :]
    #strip f12 from repeats of the same face as one of its neighbours (at each side: 3 sides). Hence, the size changes from (nfaces,3,2) to (nfaces,3)
    faceindex_eye = np.tile( np.arange(nfaces)[:, np.newaxis, np.newaxis], (1, 3, 2))
    assert faceindex_eye.shape == (nfaces, 3, 2)
    f12 = f12 + 1
    f12[f12 == faceindex_eye+1] = 0  # np.nan
    # print f12[0:5, :, :]
    assert np.allclose(np.prod(f12, axis=2), 0)  # check of each row has at least a zero
    f_uniq = np.sum(f12, axis=2)  # 0s will not be added
    #print f_uniq[0:5,:]
    # f_uniq now contains the neighbour of each face. one face at each side of each face: nfaces x 3 -> face
    assert np.sum(f_uniq.ravel() == 0) == 0  # all faces have at least one neighbour at each side.
    return f_uniq - 1  # fix back the indices


global q
q=0
def process2_vertex_resampling_relaxation(verts, facets, iobj, c=2.0):
    """ @param: c: resampling based on curvature. c= => uniform, regardless of curvature."""
    global q
    q += 1
    #print "q"*1000, q
    #facets = fix_faces_3div2(facets)
    centroids = compute_centroids(verts, facets)
    centroid_normals_normalized = compute_centroid_gradients(centroids, iobj, normalise=True)
    from mesh_utils import make_neighbour_faces_of_vertex
    faceslist_neighbours_of_vertex = make_neighbour_faces_of_vertex(facets)
    face_neighbours_of_faces_Fx3 = build_faces_of_faces(facets)
    assert not np.any(np.isnan(verts.ravel()))  # fine
    new_verts = vertex_resampling(verts, faceslist_neighbours_of_vertex, face_neighbours_of_faces_Fx3, centroids, centroid_normals_normalized, c=c)

    if VISUALISE_RELAXATION_STEPS:
        display_simple_using_mayavi_2( [(verts, facets), (new_verts, facets)],
           mayavi_wireframe=[False, True,], opacity=[0.2, 1,], gradients_at=None, separate_panels=False, gradients_from_iobj=None,
           #minmax=(RANGE_MIN,RANGE_MAX),
           )

    assert not np.any(np.isnan(new_verts.ravel()))  # fails
    return new_verts, facets, centroids  # why does it return facets?



def is_intarray(x):
    return issubclass(x.dtype.type, np.integer)

global failure_pairs
failure_pairs = []

#from mesh1.py
def vertex_resampling(verts, faceslist_neighbours_of_vertex, face_neighbours_of_faces_Fx3, centroids, centroid_normals, c=2.0):

    """ faceslist_neighbours_of_vertex: *** """
    assert is_intarray(face_neighbours_of_faces_Fx3)

    def kij(i, j):
        """ Returns (1/r * Theta), a measure of curvature.
        Theta is the angle between two normals at centroids (dual vertices) i, j.
        The 1/r is the inverse of the distance between the pair.

        Notes:
        Normals should be already normalised (centroid normals).
        centroids should be already [projected] on the implicit surface (dual mesh is optimised).
        """

        # i,j are centroids
        assert i != j
        pi, pj = (centroids[i, 0:3], centroids[j, 0:3])

        # based on gradients. normalised.
        mi, mj = (centroid_normals[i, 0:3], centroid_normals[j, 0:3])
        assert mi.shape == (3,)
        assert mj.shape == (3,)
        assert np.abs(np.linalg.norm(mi) - 1.0) < 0.0000001
        assert np.abs(np.linalg.norm(mj) - 1.0) < 0.0000001
        mimj = np.dot(np.transpose(mi), mj)
        #mimj[mimj>1.0] = 1.0
        if mimj > 1.0:
            mimj = 1.0
        if mimj < -1.0:
            mimj = -1.0
        pipj = np.linalg.norm(pi - pj)
        #print "pipj ", pipj, "  arccos= ",np.arccos(mimj)/np.pi*180 #why is it zero??
        assert pipj == np.linalg.norm(pi - pj, ord=2)

        if pipj == 0.:
            failure_pairs.append( (i, j) )
            #raise TroubledMesh("repeated triangle")
            pipj = 1
        #CAN BE absolute ZERO! ****************
        #FIXME
        assert pipj > 0  # fails # ****

        #if (i, j) == (358, 469):
        #    print "pipj", pipj, "   == %5.20g"%(pipj)
        #    print np.arccos(mimj)
        #    print mimj
        #    #set_trace()
        kij = np.arccos(mimj) / pipj  # radians?
        assert not np.isnan(kij), "i=%d,j=%d"%(i, j)  # fails
        return kij

    def wi(i_facet, ja_facets, c):
        """
        Returns the weight of a facet i_facet.
        Adds kij of all centroids of the neighbour facets.
        ja_facets = list of centroid indices (face index).
        i_facet is a face index. """
        # todo: make a pipj matrix (fxf). Make an acos matrix (fxf). The latter is base on a matrix of gradients: fx3.
        #
        #print i_facet, ja_facets
        assert i_facet not in ja_facets
        assert len(ja_facets) == 3
        # ja_facets = neighbour facets of facet i_facet????
        ki = 0
        for j_facet in ja_facets:
            ki += kij(i_facet, j_facet)

        assert not np.isnan(ki)  # fails

        #set_trace()
        wi = 1.0 + c*ki
        # i_facet is facet (Centroid) index. j_facet is its neighbour facet (centroid). There are three j_facet for an i_facet.
        return wi
    #
    #set_trace()
    c_ = c  # 2.0  # constant
    vertex_index = 1  # vertex
    #assert vertex_index >= 0
    umbrella_facets = faceslist_neighbours_of_vertex[vertex_index]  # A list of facets: The indices of faces that vertex vertex_index belongs to.
    #print("umbrella_facets: ", umbrella_facets)
    #wa = np.zeros()
    w_list = []
    for i_facet in umbrella_facets:
        # neighbour facet i_facet of Vertex vertex_index
        #three_facets = filter(lambda idx: idx != i_facet, umbrella_facets)
        three_facets = face_neighbours_of_faces_Fx3[i_facet, :]
        w = wi(i_facet, three_facets, c_)  # three_facets should be neighbours of the facet i_facet
        # The weight (based on curvature) of neighbour P_i (facet i.e. centroid),
        w_list.append(w)
        #todo: sparse matrix: w[vi=vertex_index, f2=i_facet] = w
        #todo: store in ...
    #print "w_list ",w_list
    #
    #del w_list
    #w seems tobe calculated fine. next: store w_i and cache them for adaptive resampling, for which we need to normalise it across the neighbours.
    nfaces = centroids.shape[0]
    wi_total_array = np.zeros((nfaces,))
    for i_facet in range(nfaces):
        three_facets = face_neighbours_of_faces_Fx3[i_facet, :]
        w = wi(i_facet, three_facets, c_)
        assert not np.isnan(w)  # fails

        wi_total_array[i_facet] = w
    #print wi_total_array
    # The weights are prepared. Now let's resample vertices
    assert not np.any(np.isnan(wi_total_array.ravel()))  # fails

    vertex_index = 1
    #todo: umbrella_Facets = sparse matrix
    #umbrella_facets = np.array(faceslist_neighbours_of_vertex, dtype=np.int)  #empty

    umbrella_facets = np.array(faceslist_neighbours_of_vertex[vertex_index], dtype=np.int)  # empty

    #print "umbrella_facets", umbrella_facets.shape, "****"
    assert np.allclose( wi_total_array[umbrella_facets] - np.array(w_list), 0)

    def lift_verts(verts, centroids):
        new_verts = verts.copy()
        # assign these to a sparse matrix? and  do:  M = M/normalise(M); verts = M * verts
        for vertex_index in range(verts.shape[0]):
            #print vertex_index
            #print "  len=", len(faceslist_neighbours_of_vertex)
            #print "  max=", max(faceslist_neighbours_of_vertex)

            # Note: the vertex may not exist, hence not in  faceslist_neighbours_of_vertex
            if vertex_index in faceslist_neighbours_of_vertex:
                #print faceslist_neighbours_of_vertex[vertex_index]
                umbrella_facets = np.array(faceslist_neighbours_of_vertex[vertex_index], dtype=int)
                w = wi_total_array[umbrella_facets]
                assert not np.any(np.isnan(umbrella_facets.ravel()))  # pass
                assert not np.any(np.isnan(wi_total_array.ravel()))  # fails

                assert not np.any(np.isnan(w.ravel()))  # fails
                w = w / np.sum(w)
                assert not np.any(np.isnan(w.ravel()))  # fails

                new_verts[vertex_index, :] = \
                    np.dot(w, centroids[umbrella_facets, 0:3])  # (n) * (n x 3)
            else:
                pass  # leave new_verts[vertex_index,:] unmodified
        assert not np.any(np.isnan(new_verts.ravel()))  # fails
        return new_verts

    assert not np.any(np.isnan(verts.ravel()))  # fine
    r = lift_verts(verts, centroids)
    #print r.shape
    assert not np.any(np.isnan(r.ravel()))  # fails
    #exit()
    return r






def two_bricks():
    import vectorized, example_objects
    c2 = vectorized.UnitCube1(1.)
    def rotate_scale_(iobj, scale, center, angle=0.):
        ns = vectorized
        import numpy
        m = numpy.eye(4)
        m[0,0] = 0.1
        iobj = ns.Transformed(iobj, m=m)
        iobj  \
            .resize(scale) \
            .move(center[0], center[1], center[2])
        if angle != 0.:
            iobj.rotate(angle, along=make_vector4(1, 1, 1), units="deg")
        return iobj

    c2 = rotate_scale_(c2, 2., [1,1,1])
    iobj = vectorized.CrispUnion( example_objects.rcube_vec(1.), c2 )
    return iobj

def make_bricks():
    import vectorized
    import example_objects
    c2 = vectorized.UnitCube1(1.)

    def rotate_scale_(iobj, scale, center, angle=0.):
        ns = vectorized
        import numpy
        m = numpy.eye(4)
        m[0,0] = 0.1
        iobj = ns.Transformed(iobj, m=m)
        iobj  \
            .resize(scale) \
            .move(center[0], center[1], center[2])
        if angle != 0.:
            iobj.rotate(angle, along=make_vector4(1, 1, 1), units="deg")
        return iobj

    c2 = rotate_scale_(c2, 2., [1, 1, 1])
    c3 = example_objects.rcube_vec(1.)
    iobj = c2 #
    iobj = vectorized.CrispUnion( c3, c2 )
    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-3, +5, 0.2 * 2.)
    return iobj, RANGE_MIN, RANGE_MAX, STEPSIZE




from visual5 import *

#python -O -m kernprof -v -l ohtake_belyaev_5.py
#@profile

def compute_facets_subdivision_curvatures_old(verts, facets, iobj):
    """ Calculates a measure of deviation of the Triangle from object gradients.
    returns: curvature for all triangles.
    This function does not create the subdivisions.
    It just computes. The function subdivide_multiple_facets() does the actual subdivision. """
    facet_areas, facet_normals = compute_triangle_areas(verts, facets, return_normals=True)

    #simple_histogram(facet_areas, "facet_areas")

    nf = facets.shape[0]
    assert facet_areas.shape == (nf,)
    assert facet_normals.shape == (nf, 3)
    assert np.all(np.logical_not(np.isnan(facet_areas[np.logical_not(np.isnan(np.linalg.norm(facet_normals, axis=1)))])))
    #some edges are repeated
    degenerate_faces = np.isnan(facet_areas)

    nn = np.isnan(facet_normals[np.logical_not(degenerate_faces),:])
    #print nn.shape
    #print np.any(nn, axis=1).shape, "or"
    #print facet_areas[np.any(nn, axis=1)]

    assert np.all(np.logical_not(np.isnan(facet_areas.ravel()))), "facet_areas: never nan. But can be zero."

    assert np.all(np.isnan(facet_areas[degenerate_faces]))
    assert np.all(np.logical_not(np.isnan(facet_areas[np.logical_not(degenerate_faces)])))
    assert np.all(np.isnan(facet_normals[degenerate_faces, :]))
    if mesh_correction:
        assert np.all(np.logical_not(facet_areas == 0.)), "Facet area zero."
        assert np.all(np.logical_not(np.isnan(facet_normals[np.logical_not(degenerate_faces),:])))

    centroidmaker_matrix = np.array([
        [1, 0, 0, 1, 0, 1],  # 035
        [0, 1, 0, 1, 1, 0],  # 314
        [0, 0, 1, 0, 1, 1],  # 542
        [0, 0, 0, 1, 1, 1],  # 345
        ]) / 3.

    subdiv_vert_matrix = np.array([
        [1.,   0.,  0.],  # 0
        [0.,   1.,  0.],  # 1
        [0.,   0.,  1.],  # 2

        [0.5,  0.5,  0],  # 3
        [0,  0.5,  0.5],  # 4
        [0.5,  0,  0.5]   # 5
        ])

    #check_degenerate_faces1(verts, facets, degenerate_faces)

    curvatures_array = np.zeros((nf,))
    for fi in range(nf):
        n = facet_normals[fi, :]  # n: (3,)
        nlen = np.linalg.norm(n)
        if nlen > 0:
            n = n / nlen
        if mesh_correction:
            if degenerate_faces[fi]:
                curvatures_array[fi] = 0
                continue
            if np.isnan(n[0]+n[1]+n[2]):
                curvatures_array[fi] = 0
                continue
        triangle = verts[facets[fi, :], :]  # numverts x 3
        assert not np.any(np.isnan(triangle.ravel()))
        assert triangle.shape == (3, 3)
        VVV = triangle  # (nv=3) x 3
        assert not np.any(np.isnan(VVV.ravel()))
        assert not np.any(np.isnan(subdiv_vert_matrix.ravel()))

        m0123 = np.dot( centroidmaker_matrix, np.dot(subdiv_vert_matrix, VVV) )
        assert m0123.shape == (4, 3)
        assert not np.any(np.isnan(m0123.ravel()))
        subdiv_centroids = m0123
        numsubdiv = 4
        subdiv_centroids4 = np.concatenate( (subdiv_centroids, np.ones((numsubdiv, 1))), axis=1)
        if degenerate_faces[fi]:
            print "WARNING: degenerate triangle", fi, " = ",facets[fi,:]
        else:
            assert not degenerate_faces[fi]
        #print "*",subdiv_centroids4.shape  # too many lines printed
        #Problem: there are NaNs
        assert not np.any(np.isnan(subdiv_centroids4.ravel()))
        #check_vector4_vectorized(subdiv_centroids4)
        mm = - iobj.implicitGradient(subdiv_centroids4)[:, 0:3]
        assert mm.shape == (4, 3)
        nn = np.linalg.norm(mm, axis=1)
        nn_tile = np.tile(nn[:,np.newaxis], (1, 3))  # mm: 4 x 3
        if mesh_correction:
            nn_tile[nn_tile<0.00000001] = 100000.
        mm = mm / nn_tile
        mm = mm.transpose()  # 3x4

        #n can be nan. fixme.
        #QD solution:
        if np.any(np.isnan(n)):
            e = 0.
        else:
            e = facet_areas[fi] * np.sum(1. - np.abs(np.dot(n, mm))) / 4.  # sum(,x4)
            if e < 0:
                set_trace()
            assert e >= 0

        # The only reason for NaN should be the exactly-zero gradients

        #assert np.all(np.dot(n, mm) > -0.0000001 ), "ingrown normal!"
        #e = np.sum(1 - np.abs(np.dot(n, mm)))   # sum(,x4)   #forgot the abs!
        curvatures_array[fi] = e
        #print e,
        if not mesh_correction:  # NAN is allowed (and used) for mesh correction
            assert not np.isnan(e)
        #if np.isnan(e):
        #    print e
        #    #print facet_areas[fi] , np.sum(1. - np.abs(np.dot(n, mm))) / 4., np.abs(np.dot(n, mm))
        #    #print "n,mm", np.dot(n, mm), n, mm
        #    print "n", n
        #if e<0:
        #    set_trace()

        if fi % 100 == 0:
            print fi, "*   \r", ;import sys; sys.stdout.flush()

    # The only reason for NaN should be the exactly-zero gradients
    # assert np.sum(np.isnan(curvatures_array)) == 0
    l = curvatures_array[np.logical_not(np.isnan(curvatures_array))].tolist()
    l.sort()
    print "curvature range: min,max = ", l[0], l[-1]   # 3.80127650325e-08, 0.0240651184551
    bad_facets_count = np.sum(degenerate_faces)
    #assert bad_facets_count == 0

    #assert np.sum(np.isnan(curvatures_array)) == 0, "NaN"

    #simple_histogram(curvatures_array, "curvatures_array")

    return curvatures_array, bad_facets_count


def augment4(x):
    return np.concatenate((x, np.ones((x.shape[0], 1))), axis=1)

def compute_facets_curvatures_vectorized(verts, facets, iobj):
    facet_areas, facet_normals = compute_triangle_areas(verts, facets, return_normals=True)

    #simple_histogram(facet_areas, "facet_areas")

    nf = facets.shape[0]
    assert facet_areas.shape == (nf,)
    assert facet_normals.shape == (nf, 3)
    assert np.all(np.logical_not(np.isnan(facet_areas[np.logical_not(np.isnan(np.linalg.norm(facet_normals, axis=1)))])))
    #some edges are repeated
    degenerate_faces = np.isnan(facet_areas)

    nn = np.isnan(facet_normals[np.logical_not(degenerate_faces),:])
    assert np.all(np.logical_not(np.isnan(facet_areas.ravel()))), "facet_areas: never nan. But can be zero."
    assert np.all(np.isnan(facet_areas[degenerate_faces]))
    assert np.all(np.logical_not(np.isnan(facet_areas[np.logical_not(degenerate_faces)])))
    assert np.all(np.isnan(facet_normals[degenerate_faces, :]))
    if mesh_correction:
        assert np.all(np.logical_not(facet_areas == 0.)), "Facet area zero."
        assert np.all(np.logical_not(np.isnan(facet_normals[np.logical_not(degenerate_faces),:])))

    centroidmaker_matrix = np.array([
        [1, 0, 0, 1, 0, 1],  # 035
        [0, 1, 0, 1, 1, 0],  # 314
        [0, 0, 1, 0, 1, 1],  # 542
        [0, 0, 0, 1, 1, 1],  # 345
        ]) / 3.

    subdiv_vert_matrix = np.array([
        [1.,   0.,  0.],  # 0
        [0.,   1.,  0.],  # 1
        [0.,   0.,  1.],  # 2

        [0.5,  0.5,  0],  # 3
        [0,  0.5,  0.5],  # 4
        [0.5,  0,  0.5]   # 5
        ])

    #check_degenerate_faces1(verts, facets, degenerate_faces)

    triangles = verts[facets[:, :], :]  # nf x 3 x 3
    #degenerate_faces: (0:nf)->Bool

    if True:
        nf = triangles.shape[0]
        #np.swapaxes(triangles, 1, 2)
        subdivmat = np.dot(centroidmaker_matrix, subdiv_vert_matrix)  # 4x3
        #m0123 = np.dot( subdivmat, triangles )  # 4x3 *  (nfx3x3)    #4x1784x3??
        #m0123 = np.swapaxes(m0123, 1, 2)
        ##m0123[degenerate_faces, :, :] = ?
        q = np.tensordot(triangles, subdivmat, axes=(1, 1))  #nf x 3 x 4
        del triangles
        q = np.swapaxes(q, 1, 2)  # nf x 4 x 3
        assert q.shape[0]*q.shape[1] == nf*4
        q = q.reshape(q.shape[0]*q.shape[1], 3)
        #set_trace()
        assert q.shape == (nf*4, 3)
        q4 = augment4(q); del q
        check_vector4_vectorized(q4)
        mm = - iobj.implicitGradient(q4)[:, 0:3] ; del q4
        assert mm.shape == (nf*4, 3)
        mmt = mm.reshape(nf, 4, 3); del mm  # nfx4x3
        #set_trace()
        mmt_norm = np.linalg.norm(mmt, axis=2, keepdims=True)
        mmt_norm[mmt_norm < 0.00001] = 1.
        mmt_hat = mmt / mmt_norm ; del mmt; del mmt_norm

        n_norm = np.linalg.norm(facet_normals, axis=1, keepdims=True)
        n_norm[n_norm < 0.00001] = 1.
        n_hat = facet_normals / n_norm; del n_norm

        assert mmt_hat.shape[2] == 3
        assert n_hat.shape[1] == 3
        #nm = np.tensordot(mmt_hat, n_hat, axes=(2, 1))  # Fx4x3, Fx3 -> ?
        #set_trace()
        #assert nm.shape == ()
        #np.sum(, axis=)

        #mmt_hat_Fx3x4 = np.swapaxes(mmt_hat_Fx3x4, 1, 2)  # Fx3x4
        #nm = np.dot(n_hat, mmt_hat_Fx3x4)  # Fx(3), Fx(3)x4  -> ? FxFx4
        #        # np.dot(): For N dimensions it is a sum product over the last axis of a and the second-to-last of b:

        #n_hat, mmt_hat_Fx3x4
        nm = np.sum(n_hat[:, np.newaxis, :] * mmt_hat, axis=2)  # Fx1x3 * Fx4x3 -> Fx4
        #np.sum(1.-np.abs(nm), axis=1)
        curvatures_array = facet_areas * (1. - np.sum(np.abs(nm), axis=1) / 4.)
        curvatures_array[degenerate_faces] = 0.
        bad_facets_count = np.sum(degenerate_faces)
        #print "hey"
        #set_trace()
        return curvatures_array, bad_facets_count


    #set_trace()
    curvatures_array = np.zeros((nf,))
    for fi in range(nf):
        n = facet_normals[fi, :]  # n: (3,)
        nlen = np.linalg.norm(n)
        if nlen > 0:
            n = n / nlen
        if mesh_correction:
            if degenerate_faces[fi]:
                curvatures_array[fi] = 0
                continue
            if np.isnan(n[0]+n[1]+n[2]):
                curvatures_array[fi] = 0
                continue

        if degenerate_faces[fi]:
            print "WARNING: degenerate triangle", fi, " = ",facets[fi,:]
        else:
            assert not degenerate_faces[fi]

        #triangle = verts[facets[fi, :], :]  # numverts x 3
        triangle = triangles[fi, :, :]  # numverts x 3
        assert not np.any(np.isnan(triangle.ravel()))
        assert triangle.shape == (3, 3)  # nv=3
        assert not np.any(np.isnan(triangle.ravel()))
        assert not np.any(np.isnan(subdiv_vert_matrix.ravel()))

        #m0123 = np.dot( centroidmaker_matrix, np.dot(subdiv_vert_matrix, triangle) )
        #m0123 = np.dot( np.dot(centroidmaker_matrix, subdiv_vert_matrix), triangle )
        subdivmat = np.dot(centroidmaker_matrix, subdiv_vert_matrix)
        m0123 = np.dot( subdivmat, triangle )
        assert m0123.shape == (4, 3)
        assert not np.any(np.isnan(m0123.ravel()))
        subdiv_centroids = m0123
        numsubdiv = 4
        subdiv_centroids4 = np.concatenate( (subdiv_centroids, np.ones((numsubdiv, 1))), axis=1)

        #print "*",subdiv_centroids4.shape  # too many lines printed
        #Problem: there are NaNs
        assert not np.any(np.isnan(subdiv_centroids4.ravel()))
        #check_vector4_vectorized(subdiv_centroids4)
        mm = - iobj.implicitGradient(subdiv_centroids4)[:, 0:3]
        assert mm.shape == (4, 3)
        nn = np.linalg.norm(mm, axis=1)
        nn_tile = np.tile(nn[:,np.newaxis], (1, 3))  # mm: 4 x 3
        if mesh_correction:
            nn_tile[nn_tile<0.00000001] = 100000.
        mm = mm / nn_tile
        mm = mm.transpose()  # 3x4

        #n can be nan. fixme.
        #QD solution:
        if np.any(np.isnan(n)):
            e = 0.
        else:
            e = facet_areas[fi] * np.sum(1. - np.abs(np.dot(n, mm))) / 4.  # sum(,x4)
            if e < 0:
                set_trace()
            assert e >= 0

        # The only reason for NaN should be the exactly-zero gradients

        #assert np.all(np.dot(n, mm) > -0.0000001 ), "ingrown normal!"
        #e = np.sum(1 - np.abs(np.dot(n, mm)))   # sum(,x4)   #forgot the abs!
        curvatures_array[fi] = e
        #print e,
        if not mesh_correction:  # NAN is allowed (and used) for mesh correction
            assert not np.isnan(e)
        #if np.isnan(e):
        #    print e
        #    #print facet_areas[fi] , np.sum(1. - np.abs(np.dot(n, mm))) / 4., np.abs(np.dot(n, mm))
        #    #print "n,mm", np.dot(n, mm), n, mm
        #    print "n", n
        #if e<0:
        #    set_trace()

        if fi % 100 == 0:
            print fi, "*   \r", ;import sys; sys.stdout.flush()

    # The only reason for NaN should be the exactly-zero gradients
    # assert np.sum(np.isnan(curvatures_array)) == 0
    l = curvatures_array[np.logical_not(np.isnan(curvatures_array))].tolist()
    l.sort()
    print "curvature range: min,max = ", l[0], l[-1]   # 3.80127650325e-08, 0.0240651184551
    bad_facets_count = np.sum(degenerate_faces)
    #assert bad_facets_count == 0

    #assert np.sum(np.isnan(curvatures_array)) == 0, "NaN"

    #simple_histogram(curvatures_array, "curvatures_array")

    return curvatures_array, bad_facets_count


def compute_facets_subdivision_curvatures(verts, facets, iobj):
    (curvatures_array2, bad_facets_count2) = compute_facets_curvatures_vectorized(verts, facets, iobj)
    TEST_VECTORIZED_CURVATURE_COMPUTAITON = False
    if TEST_VECTORIZED_CURVATURE_COMPUTAITON:
        (curvatures_array1, bad_facets_count1) = compute_facets_subdivision_curvatures_old(verts, facets, iobj)
        assert np.allclose(curvatures_array1, curvatures_array2)
        assert np.allclose(bad_facets_count1, bad_facets_count2)
        print "curvature: tests passed"
        print "*"*300
    return curvatures_array2, bad_facets_count2

#global jk
#jk = 0


def propagated_subdiv(facets, subdivided_edges):
    """ Reports the indices triangles that are to be subdivided
    as a propagation of a previous subdivision.
    It returns separatesly the (not-yet-subdivided) triangles
    with 1,2 and 3 sibdivided sides, in a dictionary.
    returns: edges_need_subdivision: those edges that still exist in mesh that need further subdivision.
    But the type of subdivision will be determined based on the map propag_dict, which organises them based on the number of edges that need subdivision in the traiangle they belong.
    subdivided_edges are both subdivided and not subdivided: they are subdivided previously but need to be subdivided again,
    because they belong to two triangles. The sceond triangle may not have been subdivided yet. So they remain in the mesh as unsubdivided, although they are subdivided previously.
    Each edge belongs to two triangles, hence this contradiction sometimes exit.
    The function returns the edges that remain to be dealt with."""

    if len(subdivided_edges) == 0:
        subdivided_edges = np.zeros((0, 2), dtype=np.int64)
    else:
        subdivided_edges = np.asarray(subdivided_edges)  # doesnt work if empty
    #print subdivided_edges.shape
    assert subdivided_edges.shape[1] == 2
    subdivided_edges.sort(axis=1)

    BB = np.array([1L, B], dtype=np.int64)
    subdivided_edges_codes = np.dot(subdivided_edges, BB)  # can be large
    assert subdivided_edges_codes.dtype == np.int64
    assert subdivided_edges_codes.size == 0 or np.min(subdivided_edges_codes) >= 0
    assert np.max(facets, axis=None) < B

    fc0 = facets.copy()
    f0 = facets[:, np.newaxis, [0, 1]]
    f1 = facets[:, np.newaxis, [1, 2]]
    f2 = facets[:, np.newaxis, [0, 2]]
    f012 = np.concatenate((f0, f1, f2), axis=1)
    f012.sort(axis=2)
    assert np.all(fc0.ravel() == facets.ravel())  # no copy() needed
    all_edges_codes = np.dot(f012, BB)  # *x3

    # now look for subdivided_edges_codes in all_edges_codes

    # todo: refactor: intersec seems to be not used anymore
    intersec = np.intersect1d(all_edges_codes, subdivided_edges_codes)
    # x_ indicates those edges taht are in the mesh's edges, hence remain there.
    # x_ is the bottleneck (both in terms of data and in terms of performance)
    x_ = np.lib.arraysetops.in1d(all_edges_codes, intersec)  # elements of A, A.ravel[x_], that are in B
    assert np.prod(all_edges_codes.shape) == x_.size
    assert np.sum(x_) == intersec.size  # 417
    #assert each_of all_edges_codes.ravel()[x_] in intersec
    edges_need_subdivision = all_edges_codes.ravel()[x_]  # all edges_need_subdivision are in intersec

    sides_booleans_Fx3 = x_.reshape(all_edges_codes.shape)  # *x3
    numsides = np.sum(sides_booleans_Fx3, axis=1)  # numsides_needsubdivision: number of sides that need subdivision. index=face index
    assert sides_booleans_Fx3.shape == (facets.shape[0], 3)

    #now I need those all_edges_codes (i.e. sides_booleans_Fx3) for which numsides==1
    #sides_1 = sides_booleans_Fx3[numsides==1, :] for c==1
    propag_dict = {}
    for c in range(1, 4):
        # Range starts with 1 because we only propagate triangles with subdivided 1,2,3 sides.
        idx = np.nonzero(numsides == c)[0]
        propag_dict[c] = idx
        del idx

    return propag_dict, edges_need_subdivision




def subdivide_multiple_facets(verts_old, facets_old, tobe_subdivided_face_indices, midpoint_map):
    """
    midpoint_map is modified (is input and output).
    midpoint_map is a dictionary that given an edge's unique_int_id, gives you the vertex in the midpoint. It may contain midpoints that are not used anymore.
    Use compute_facets_subdivision_curvatures() to calculate tobe_subdivided_face_indices.
    Does not remove vertices => will be valid. But vertices will change: new elements will be appended to it.
    Returns: new vertices and faces.
    Returns: presubdivision_edges: The edges that have been removed.
    Theses edges will be invalid after this function.
    Any such edges (those that remain somewhere else) also has to be later removed (and replaced by two subdivided ones) from the rest of the mesh.
    This will be used for propagating the subdivision to triangles that their edges are not valid anymore.

    The result will have T-junctions which should be resolved by further propagation of subdivisions.
    """

    # todo: store subdivided gradients (on top of centroids), to avoid unnecessary calculations. When updating vettices, remove the caches.
    # todo: avoid recomputing


    #TODO: INDEX PROBLEM

    centroidmaker_matrix = np.array([
        [1, 0, 0, 1, 0, 1],  # 035
        [0, 1, 0, 1, 1, 0],  # 314
        [0, 0, 1, 0, 1, 1],  # 542
        [0, 0, 0, 1, 1, 1],  # 345
        ]) / 3.

    DIP = 0.05*0
    subdiv_vert_matrix = np.array([
        [1.,   0.,  0.],  # 0
        [0.,   1.,  0.],  # 1
        [0.,   0.,  1.],  # 2
        #todo: remove unnecessary points

        [0.5/(1.+DIP),  0.5/(1.+DIP),  DIP/(1.+DIP)],  # 3
        [DIP/(1.+DIP),  0.5/(1.+DIP),  0.5/(1.+DIP)],  # 4
        [0.5/(1.+DIP),  DIP/(1.+DIP),  0.5/(1.+DIP)]   # 5
        ])  # .transpose()

    global trace_subdivided_facets
    trace_subdivided_facets = []

    verts_old_old = verts_old.copy()

    # Allocate space for the new faces and vertices
    provisional_new_verts_count = 3*len(tobe_subdivided_face_indices)
    provisional_new_facets_count = 3*len(tobe_subdivided_face_indices)
    nverts_old = verts_old.shape[0]
    nfaces_old = facets_old.shape[0]
    new_verts = np.zeros((nverts_old + provisional_new_verts_count, 3), dtype=float)
    new_facets = np.zeros((nfaces_old + provisional_new_facets_count, 3), dtype=int)
    new_verts[:nverts_old, :] = verts_old
    new_facets[:nfaces_old, :] = facets_old

    #some variable names:
    # mapped_midvertices

    #print "Subdividing:"

    #The edges that need to be divided because of the neighbouring triangle subdivided
    presubdivision_edges = []
    redundancy_counter = 0

    new_vertex_counter = nverts_old
    new_facet_counter = nfaces_old
    for subdiv_i in range(len(tobe_subdivided_face_indices)):

        fi = tobe_subdivided_face_indices[subdiv_i]
        triangle_old = verts_old[facets_old[fi, :], :]  # numverts x 3
        assert triangle_old.shape == (3, 3)

        # output: v345_xyz
        _vxyz_0123 = np.dot(subdiv_vert_matrix, triangle_old)  # not efficient
        assert _vxyz_0123.shape == (6, 3)
        v345_xyz = _vxyz_0123[3:6, :]  # only pick the new ones
        del _vxyz_0123

        original_facet_index = fi
        e0 = new_facets[original_facet_index, [0, 1]]
        e1 = new_facets[original_facet_index, [1, 2]]
        e2 = new_facets[original_facet_index, [2, 0]]
        presubdivision_edges .append(tuple(e0.tolist()))
        presubdivision_edges .append(tuple(e1.tolist()))
        presubdivision_edges .append(tuple(e2.tolist()))

        e012 = np.vstack((e0, e1, e2))
        e012.sort(axis=1)
        BB = np.array([1L, B], dtype=np.int64)
        all_edges_triples = np.dot(e012, BB)
        assert all_edges_triples.dtype == np.int64
        assert all_edges_triples.size == 0 or np.min(all_edges_triples) >= 0
        assert np.max(new_facets, axis=None) < B

        #avoid becasue it is redundant
        avoid_which = np.zeros((3,), dtype=np.bool) + False

        idx_counter = new_vertex_counter
        #actual_mapped_midvertices
        actual_3_vertices = np.zeros((3,), dtype=np.int64)
        for i in range(3):
            if all_edges_triples[i] in midpoint_map:
                avoid_which[i] = True
                actual_3_vertices[i] = midpoint_map[all_edges_triples[i]]
                #mapped_midvertices[i] = -1  # for debug
                redundancy_counter += 1
            else:
                assert avoid_which[i] == False
                #x = mapped_midvertices[i]  # wrong!
                #x = idx_counter
                midpoint_map[all_edges_triples[i]] = idx_counter
                idx_counter += 1
                actual_3_vertices[i] = idx_counter - 1  # the new vertex

        use_which = np.logical_not(avoid_which)
        n1 = new_vertex_counter
        n2 = idx_counter
        assert n2 == n1 + np.sum(use_which)
        new_vertex_counter = n2
        new_verts[n1:n2, :] = v345_xyz[use_which, :]

        #Output: mini_faces
        # adding new verts and facets
        # indices of the original vertices.
        _v345 = actual_3_vertices.tolist()
        _v012 = facets_old[fi, :].tolist()
        # facet's vertex indices
        _v012345 = np.array(_v012 + _v345, dtype=int)
        _mini_faces_l = [[0, 3, 5], [3, 1, 4], [5, 4, 2], [3, 4, 5]]  # 0,3,1,4,2,5
        mini_faces = _v012345[np.array(_mini_faces_l, dtype=int)]
        del _v012
        del _v345
        del _v012345

        # facet's vertex indices
        new_facets[original_facet_index, :] = mini_faces[0, :]
        new_facets[new_facet_counter:(new_facet_counter+3), :] = mini_faces[1:(1+3), :]
        assert mini_faces.shape[0] == (1+3)
        trace_subdivided_facets += range(new_facet_counter, (new_facet_counter+3)) + [fi]  # include the face which reuses the old face's index
        # trace_subdivided_facets will contain indices of faces

        new_facet_counter += 3

        if subdiv_i % 100 == 0:
            print subdiv_i , "       \r", ;import sys; sys.stdout.flush()
    #print " "*10

    assert new_verts.shape[0] - new_vertex_counter == redundancy_counter
    new_verts = new_verts[:new_vertex_counter, :]
    #quick_vis(noisy(new_verts, 0.05), new_facets, range(new_facets.shape[0]))
    assert np.max(new_facets.ravel()) < new_verts.shape[0]
    assert new_verts.shape[0] == new_vertex_counter
    assert new_facets.shape[0] == new_facet_counter
    assert provisional_new_verts_count+nverts_old-redundancy_counter == new_vertex_counter, "vector consistency"
    assert provisional_new_facets_count+nfaces_old == new_facet_counter, "face consistency"
    assert len(trace_subdivided_facets) == 0 or np.max(np.array(trace_subdivided_facets)) < new_facet_counter

    return new_verts, new_facets, presubdivision_edges
    # todo: refactor "presubdivision_edges"
    #presubdivision_edges: the cut edges, the old edges that should not exist anymore


def testcase_square():
    v = np.array([[0, 0, 0], [0, 1, 0], [1, 1, 0], [1, 0, 0] ])*10 / 4.
    f = np.array([[0, 1, 2], [0, 2, 3]])
    return v, f

def testcase_cube():
    faces_xyz = []
    for d in range(3):
        dims = np.array([(d+1) % 3, (d+2) %3, d], dtype=np.int)  # new dimensions that x, y, z are mapped into
        for s in range(2):
            side_verts = []
            for x in range(2):
                for y in range(2):
                    xyz = np.array([x, y, s], dtype=np.int)[dims]
                    side_verts.append(tuple(xyz))
            #reorder = np.array([0,1,3,2,0])
            #reorder = [[0, 1], [1, 3], [3, 2], [2, 0]]
            reorder = [[0, 1, 3], [3, 2, 0]]
            verts_reorders = map( lambda face: [side_verts[face[0]], side_verts[face[1]], side_verts[face[2]] ] , reorder)
            faces_xyz += verts_reorders
    #print faces_xyz
    vertdict = {}
    vert_index_dict = {}
    pure_verts = []
    faces = []
    for i in range(len(faces_xyz)):
        f = faces_xyz[i]
        face1 = []
        for v in f:
            key = str(v)
            if key in vertdict:
                idx = vert_index_dict[key]
            else:
                pure_verts.append(v)
                idx = len(pure_verts) - 1
                vert_index_dict[key] = idx
            vertdict[key] = v
            face1.append(idx)
            del idx
        faces.append(face1)


    #for v in  vertdict:
    #    print v
    #print vertdict
    #print vert_index_dict
    #print pure_verts
    for i in range(len(faces_xyz)):
        f = faces_xyz[i]
        for v in f:
            xyz = vertdict[str(v)]
    verts = pure_verts
    #print
    #print faces
    v =  np.array(verts)
    f =  np.array(faces, dtype=np.int64)
    return v, f


#from ipdb import set_trace

def check_pairs(facets): #yet another!!
    #repeated code

    e0 = facets[:, np.newaxis, [0, 1]]
    e1 = facets[:, np.newaxis, [1, 2]]
    e2 = facets[:, np.newaxis, [2, 0]]
    #e012 = np.vstack((e0, e1, e2))
    e012 = np.concatenate((e0, e1, e2), axis=1)  # n x 3 x 2
    e012.sort(axis=2)
    BB = np.array([1, B], dtype=np.int64)
    all_edges_triples = np.dot(e012, BB)  # n x 3
    assert all_edges_triples.dtype == np.int64
    assert all_edges_triples.size == 0 or np.min(all_edges_triples) >= 0
    assert np.max(facets, axis=None) < B
    #print all_edges_triples
    all_edge_triples_ravel = all_edges_triples.ravel()


    e = all_edge_triples_ravel.copy()
    e.sort()
    d = np.diff(e)
    #print d
    #set_trace()
    assert np.all(d[::2] == 0)
    assert np.all(d[1::2] != 0)
    #print "ok"
    #exit()


def pause():
    import sys
    sys.stdout.write(">")
    sys.stdin.readline()


def isomorphic(a, b):
    if not np.ndim(a) == np.ndim(b):
        return False
    if not a.shape == b.shape:
        return False
    #if not np.all(a == b):
    #    return False
    assert a.dtype.type != np.bool
    assert b.dtype.type != np.bool
    return True


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


def subdivide_1to2_multiple_facets(facets, edges_with_1_side, midpoint_map, careful_for_twosides=True):
    """list_edges_with_1_side contains the edges only. The face should be extracted in this function.
    returns: faces.
    careful_for_twosides: whe ntwo sides are being asked for subdivision. """
    #todo: copy some code from propagated_subdiv()
    #check which of these edges still exist in faces. (Each should be there only once. In this context.)
    #Some edges_with_1_side may not be in facets. They are already subdivided twice.
    #remove them and add more.
    #refactor the code copied from propagated_subdiv() into function
    #need to also get the new points. oops!! damn.
    #
    #There is a guarantee that all faces that the edges_with_1_side belong to, have exactly one edge from this list.
    #Note that these edges should be already subdivided


    #yes of course all of them are there in it
    #All e in edges_with_1_side, =>, e in midpoint_map

    assert type(edges_with_1_side) == np.ndarray
    assert edges_with_1_side.size == 0 or np.min(edges_with_1_side) > 0
    el = filter(lambda e: not e in midpoint_map, edges_with_1_side)
    if not len(el) == 0:
        for e in midpoint_map:
            print midpoint_map[e],
        print len(el)
        print "error"
        exit()
    assert len(el) == 0, "assert edges_with_1_side is subset of midpoint_map"

    all_edges_triples, dummy = get_edge_code_triples_of_mesh(facets)

    all_edge_triples_ravel = all_edges_triples.ravel()  # is a view
    #print all_edges_triples.shape
    assert all_edges_triples.shape[1] == 3

    assert np.all(all_edge_triples_ravel.reshape(-1, 3) == all_edges_triples, axis=None)


    #intersec = np.intersect1d(all_edge_triples_ravel, edges_with_1_side)

    #index_of_edges_that_subdiv2

    x_ = np.lib.arraysetops.in1d(all_edge_triples_ravel, edges_with_1_side)  # elements of A, A.ravel[x_], that are in B

    x_1x3 = x_.reshape(3, -1)
    assert np.all(x_1x3.ravel() == x_, axis=None)

    #x3__a = x_.reshape(3, -1)  # wrong. bug
    x3__b_Fx3 = x_.reshape(-1, 3)

    #how many sides are requested to be subdivided
    facemultiplicity = np.sum(x3__b_Fx3, axis=1)

    #Dont want to subdivide 1->2
    #bad2 = np.all(np.sum(x_.reshape(3, -1), axis=0) > 1)  # bug!
    bad2 = np.nonzero(facemultiplicity > 1)[0]
    if careful_for_twosides:
        if bad2.size > 0:
            print midpoint_map
            print bad2
            print facets[bad2, :]
            print edges_with_1_side
            set_trace()
        assert bad2.size == 0
    del bad2

    if careful_for_twosides:
        assert np.all(facemultiplicity <= 1)
    if careful_for_twosides:
        if not np.all(facemultiplicity <= 1):
            #print facemultiplicity.tolist()
            a = facemultiplicity
            #print np.nonzero(a > 1)
            print "FAILED"
        assert np.all(facemultiplicity <= 1)
    #print "THIS FAILS"
    del x3__b_Fx3

    #indices of all edges
    face3_idx = np.nonzero(x_)[0]
    assert np.ndim(face3_idx) == 1

    #todo(refactor): use np.argwhere()
    idx_xy = np.unravel_index(face3_idx, all_edges_triples.shape)
    #idx_xy is a tuple
    #Triangles subject to be subdivided:
    problem_face_idx = idx_xy[0]
    problem_side_idx = idx_xy[1]
    #assert np.all(problem_face_idx < 3)
    def has_repeats(x):
        y = x.copy()
        y.sort()
        return np.any(np.diff(y) == 0)

    #has repeats, becasue [2]is not resolved before
    #if has_repeats(problem_face_idx):
    #    intersec = np.intersect1d(all_edge_triples_ravel, edges_with_1_side)
    #    print intersec
    #    print facets
    #    exit()

    if careful_for_twosides:
        assert not has_repeats(problem_face_idx), "triangles subject to be subdivided"
    for ii in range(problem_face_idx.size):
        #assert all_edge_triples_ravel[problem_face_idx[ii], problem_side_idx[ii]] in edges_with_1_side
        assert all_edges_triples[problem_face_idx[ii], problem_side_idx[ii]] in edges_with_1_side
        assert problem_side_idx[ii] < 3
    #all problem_face_idx should be removed


    # The subdivided edge is between v1 and v2.
    vert_idx_1 = problem_side_idx  # The problem_side will be between (v1,v2) vertices. vert_idx_1 is not a vertex index but it is a vertex index within a face i.e. in (0,1,2).
    vert_idx_2 = (problem_side_idx + 1) % 3
    vert_idx_3 = (problem_side_idx + 2) % 3
    v1 = facets[problem_face_idx, vert_idx_1]
    v2 = facets[problem_face_idx, vert_idx_2]
    v3 = facets[problem_face_idx, vert_idx_3]

    #The sides (vertex pairs) that need to be subdivided
    subdivedges_vertex_pairs = np.vstack((v1, v2))  # size: 2 x F

    #subdivedges_vertex_pairs never actually used apart from assertion tests.

    #edge_s_codes = A flat array of all the edge codes (For the sides that should be replaced with the sibdivided ones)
    #?????????????
    edge_s_codes = all_edge_triples_ravel[x_]  # Intersection from actual edges in mesh and edges requested to get removed/subdivided.
    #subdivedges_vertex_pairs: those edges that*

    #observation: edge_s_codes is (up to morphism) a subset of, but not equal to, subdivedges_vertex_pairs
    if careful_for_twosides:
        assert np.unique(edge_s_codes).size == subdivedges_vertex_pairs.shape[1]  # before applying unique

    #if can tolerate two sides:
    #if not careful_for_twosides:
    #edge_s_codes = np.unique(edge_s_codes)
    if careful_for_twosides:
        assert np.unique(edge_s_codes).size == edge_s_codes.size
    #edge_s_codes are unique but subdivedges_vertex_pairs are not unique

    #####################################################################################################################
    tesort = subdivedges_vertex_pairs.T.copy()
    tesort.sort(axis=1)
    eid9 = np.dot(tesort, np.array([1, B], dtype=np.int64)).copy()
    assert eid9.dtype == np.int64
    assert eid9.size == 0 or np.min(eid9) >= 0
    assert np.max(facets, axis=None) < B
    eid9.sort()
    eid10 = edge_s_codes.copy()
    eid10.sort()
    assert np.all(eid10 == eid9)
    del subdivedges_vertex_pairs

    #exit()
    #map one-to-one between: edge_s_codes and midpoints_third_verts and (v1, v2, ...)
    midpoints_third_verts = np.array(map(lambda edgecode: midpoint_map[edgecode], edge_s_codes), dtype=np.int64)
    #print midpoints_third_verts
    if midpoints_third_verts.size == 0:
        return facets

    assert isomorphic(midpoints_third_verts, v1)
    #(v1,v2,v3) -> (v1, midpoints_third_verts, v3) + (midpoints_third_verts, v2, v3)
    #(v1,v3, midpoints_third_verts),  (v2,v3, midpoints_third_verts)
    new_faces1 = np.vstack(((v1, v3, midpoints_third_verts))).T  # axis is 0. .T.size = N x 3
    new_faces2 = np.vstack(((v2, v3, midpoints_third_verts))).T
    # numpy's zip()
    new_faces = np.concatenate((new_faces1[:, np.newaxis, :], new_faces2[:, np.newaxis, :]), axis=1).reshape(-1, 3)

    def sorted_copy(x):
        y = x.copy()
        y.sort()
        return y
    if careful_for_twosides:
        assert np.all(np.diff(sorted_copy(problem_face_idx)) != 0), "problem_face_idx has repeated elements"
    if careful_for_twosides:
        if not np.all(np.unique(problem_face_idx) == problem_face_idx):
            set_trace()
        assert np.all(np.unique(problem_face_idx) == problem_face_idx)
        todelete = problem_face_idx
    else:
        todelete = np.unique(problem_face_idx)
    f_rm = np.delete(facets, todelete, axis=0)
    appended_faces = np.concatenate((f_rm, new_faces), axis=0)

    #set_trace()
    return appended_faces


def simple_histogram(c, title=None, special_values=[]):
    import matplotlib.pyplot as plt
    #special_values

    plt.hist(c, 20*10)
    #plt.plot([STEPSIZE, STEPSIZE/2.], [2.5, 2.5], "r*")
    if title is not None:
        plt.title(title)
    plt.show()


global highlight
highlight = []
def do_subdivision(verts, facets, iobj, curvature_epsilon, randomized_probability=1.):
    assert not np.any(np.isnan(facets.ravel()))
    assert not np.any(np.isnan(verts.ravel()))  # fails

    #set_trace()
    print "computing curvatures"; sys.stdout.flush()
    curvatures, bad_facets_count = compute_facets_subdivision_curvatures(verts, facets, iobj)
    print "computing curvatures done."; sys.stdout.flush()
    #set_trace()

    #simple_histogram(curvatures, "curvatures")

    assert np.sum(np.isnan(curvatures)) == 0, "NaN"
    curvatures[np.isnan(curvatures)] = 0  # treat NaN curvatures as zero curvature => no subdivision

    which_facets = np.arange(facets.shape[0])[ curvatures > curvature_epsilon ]
    if randomized_probability < 1.:
        n0 = which_facets.shape[0]
        m0 = int(np.ceil(float(randomized_probability)*float(n0)))
        ridx = np.random.choice(n0, m0, replace=False)
        assert len(ridx) == m0
        which_facets = which_facets[ridx]

    print which_facets.shape
    print "applying subdivision on %d triangles."%(int(which_facets.shape[0]))
    midpoint_map = {}
    verts2, facets2, presubdivision_edges = subdivide_multiple_facets(verts, facets, which_facets, midpoint_map)
    global trace_subdivided_facets  # third implicit output

    list_edges_with_1_side = []
    while True:
        propag_dict, edges_which_in1 = propagated_subdiv(facets2, presubdivision_edges)
        facets_with_2_or_3_sides = np.concatenate((propag_dict[2], propag_dict[3]), axis=0)
        # what if those faces dont exist anymore in the next round?
        list_edges_with_1_side += [edges_which_in1]
        #print facets_with_2_or_3_sides.shape
        if facets_with_2_or_3_sides.size == 0:
            break
        verts2, facets2, old_edges2 = subdivide_multiple_facets(verts2, facets2, facets_with_2_or_3_sides, midpoint_map)
        presubdivision_edges += old_edges2  # bug fixed!

    # Finished with 2 or 3 sides.

    # Now 1 side:
    #Append all the lists in list_edges_with_1_side
    n1 = 0
    for i in range(len(list_edges_with_1_side)):
        farr = list_edges_with_1_side[i]
        assert farr.size == farr.shape[0]
        assert len(farr.shape) == 1
        n1 += farr.size
    edges_with_1_side = np.zeros((n1,), dtype=np.int64)
    n1 = 0
    for i in range(len(list_edges_with_1_side)):
        farr = list_edges_with_1_side[i]
        n2 = n1 + farr.size
        edges_with_1_side[n1:n2] = farr
        n1 = n2
    #todo: if length zero dont do it
    assert edges_with_1_side.size == 0 or np.min(edges_with_1_side) > 0

    facets2 = subdivide_1to2_multiple_facets(facets2, edges_with_1_side, midpoint_map)

    ###################
    #check_degenerate_faces(verts2, facets2, "assert")
    #build_faces_of_faces(facets2)

    print("Subdivision applied.");sys.stdout.flush()
    return verts2, facets2


from vectorized import ImplicitFunctionVectorized

class DummyImplicit(ImplicitFunctionVectorized):
    def __init__(self):
        pass

    def implicitFunction(self, p):
        check_vector4_vectorized(p)
        return p[:, 0] * 0. - 1.

    def implicitGradient(self, p):
        check_vector4_vectorized(p)
        g = p * 0. + 1.
        check_vector4_vectorized(g)
        return g

    def hessianMatrix(self, p):
        check_vector4(p)
        h = np.array([[0, 0, 0],  [0, 0, 0],  [0, 0, 0]], ndmin=2)
        return h


def test_example_meshes():

    #np.random.seed(seed=19)

    v, f = testcase_cube()
    #v, f = testcase_square()
    #print f

    v = noisy(v, 0.05)
    #v = v + (np.random.rand(v.shape[0],v.shape[1])*2.-1.) * 0.05

    def make_Cube():
        import vectorized
        ns = vectorized
        c = ns.UnitCube1()
        m2 = np.eye(4)
        m2[0, 0] = 1
        m2[1, 1] = 1
        iobj = ns.Transformed(c, m2)
        return iobj

    #iob = make_Cube()
    iob = DummyImplicit()


    #exit()
    ampl = 0.05
    v2, f2 = v, f
    for i in range(3+2+2):
        v2, f2 = do_subdivision(v2, f2, iob, -1, randomized_probability=0.3 )  # all
        #check_pairs(f2) # does not help yet.
        v2 = noisy(v2, ampl)
        quick_vis(v2, f2, range(f2.shape[0]))
        ampl = ampl * 0.5
    #print f2
    #print v2


def bigprint(text):
    print "# "+("*"*52)+"\n# "+text+"\n# "+("*"*52)
    sys.stdout.flush()


import mesh_utils

VISUALISE_RELAXATION_STEPS = True

def demo_everything(options):

    default_options = {
        "subdiv/epsilon": 1. / 1000.,
        "subdiv/iters": 0,

        "resample/iters": 1,
        "resample/c": 2.0,
        "total_iters": 15,

        "proj/mingradlen": 0.000001,  #  THRESHOLD_minimum_gradient_len =   # kill gradients smaller than this
        "proj/tol": 0.0001, #THRESHOLD_zero_interval
        "proj/maxiter": 20, # MAX_ITER
        "proj/meshnormals": True, #USE_MESH_NORMALS
        "proj/extreme_alpha": False,   #EXTREME_ALPHA = False
        "proj/maxdist_ratio": 1.0,  # Average Edge length * max
        "proj/sequence": [0, 1, 2, 3, 4,5,6],
        #absolute e

        "qem/tau": 680.,
    }

    #Result:
    #{meansqerror}

    #z = options["proj/mingradlen"]
    config = default_options.copy()
    config.update(options)

    curvature_epsilon = config["subdiv/epsilon"]  # 1. / 1000.  # *10. # a>eps  1/a > 1/eps = 2000
    VERTEX_RELAXATION_ITERATIONS_COUNT = 1 # 3
    SUBDIVISION_ITERATIONS_COUNT = 1  # 1  # 2  # 5+4
    VERTEX_RELAXATION_ADD_NOISE = False



        #
    """
        "proj/mingradlen": 0.000001,  #  THRESHOLD_minimum_gradient_len =   # kill gradients smaller than this
        "proj/zero": 0.0001, #THRESHOLD_zero_interval
        "proj/maxiter": 20, # MAX_ITER
        "proj/meshnormals": True, #USE_MESH_NORMALS
        "proj/extreme_alpha": False,   #EXTREME_ALPHA = False
        "proj/maxdist_ratio": 1.0,
        "proj/sequence": [0, 1, 2, 3, 4,5,6],
        """
        #
    """
        mc/bbox/x: (-3, 5),
        mc/bbox/y: (-3, 5),
        mc/bbox/z: (-3, 5),
        mc/step: 0.2,
        """
        #
    """
        mc/bbox/x: (-3, 5),
        mc/bbox/y: (-3, 5),
        mc/bbox/z: (-3, 5),
        mc/step: 0.2,

        mc/step/x: 0.2,
        mc/step/y: 0.2,
        mc/step/z: 0.2,

        mc/step: (0.2, 0.2, 0.2),
        """
        #
    """ todo: yield
        """
    #}

    global STEPSIZE
    from example_objects import make_example_vectorized
    iobj = make_example_vectorized(
        #"rcube_vec"  #
        "rdice_vec"  #
        #"cube_example" # problem: zero facet areas.  otherwise, it works.
        #"ell_example1"  #+
        #"bowl_15_holes"  # works too. But too many faces => too slow, too much memory. 32K?
        #"french_fries_vectorized"
        #"cyl4"
        )
    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-3, +5, 0.2*1.5/1.5  *2. /2.)
    #STEPSIZE = STEPSIZE / 2. /1.5
    #STEPSIZE = STEPSIZE * 2.

    STEPSIZE = STEPSIZE / 2.
    #STEPSIZE = STEPSIZE / 2.

    #cage "cyl4" only
    #(RANGE_MIN, RANGE_MAX, STEPSIZE) = (-32 / 2, +32 / 2, 1.92 / 4.0)

    #"bowl_15_holes" does not work with STEPSIZE= 0.2

    if True:
        iobj, RANGE_MIN, RANGE_MAX, STEPSIZE = make_bricks()
        STEPSIZE = STEPSIZE / 2.
        iobj, RANGE_MIN, RANGE_MAX, STEPSIZE = cube_with_cylinders(1)
        STEPSIZE = 0.2

        if False:
            from vectorized import Transformed
            iobj = Transformed(iobj)
            iobj.rotate(14, along=np.array([1., 1., 1., 1]), units="deg")

    print "STEPSIZE", STEPSIZE
    #set_trace()

    global giobj
    giobj = iobj

    from debug_point_collector import PointCollector
    point_collector = PointCollector(DummyImplicit)
    #iobj = PointCollector(iobj)
    #point_collector = iobj

    #(RANGE_MIN, RANGE_MAX, STEPSIZE) = (-3, +5, 0.2)
    #iobj = two_bricks()

    # ****************************************************
    # Marching Cubes
    # ****************************************************

    from stl_tests import make_mc_values_grid
    gridvals = make_mc_values_grid(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE, old="3")
    verts, facets = vtk_mc(gridvals, (RANGE_MIN, RANGE_MAX, STEPSIZE))
    print("MC calculated.");sys.stdout.flush()

    check_faces(facets)

    if False:
     display_simple_using_mayavi_2( [(verts, facets), ] * 3,
       pointcloud_list=[],
       mayavi_wireframe=[False, True, True], opacity=[0.2, 1, 0.3], gradients_at=None, separate_panels=False, gradients_from_iobj=None,
       minmax=(RANGE_MIN,RANGE_MAX),
       add_noise=[0.05, 0, 0.05], noise_added_before_broadcast=True  )
    #exit()

    #display_simple_using_mayavi_2( [(verts, facets), ],
    #   pointcloud_list=[ ], pointcloud_opacity=0.2,
    #   mayavi_wireframe=[True], opacity=[1],
    #   gradients_at=None, separate_panels=False, gradients_from_iobj=None,
    #   minmax=(RANGE_MIN,RANGE_MAX)  )
    #exit()

    for rep in range(15):

        #mesh_correction
        #take_it_easy
        #any_mesh_correction

        if mesh_correction:
          if not take_it_easy:
            #check_faces(facets)
            #facets = fix_faces_3div2(facets)
            assert not np.any(np.isnan(verts.ravel()))  # fine
            any_mesh_correction = check_degenerate_faces(verts, facets, "dontfix")
            if any_mesh_correction:
                for qq in range(5):
                    verts, facets, any_mesh_correction = check_degenerate_faces(verts, facets, "fix")
                    if not any_mesh_correction:
                        break
            if not take_it_easy:
                assert not any_mesh_correction
            check_degenerate_faces(verts, facets, "assert")

        # *******************************************************
        # VERTEX RELAXATION
        # *******************************************************

        #old_verts, old_facets = verts, facets
        assert not np.any(np.isnan(verts.ravel()))  # fine
        if False:
         display_simple_using_mayavi_2( [(verts, facets), ] * 2,
           pointcloud_list=[],
           mayavi_wireframe=[False, True], opacity=[0.2, 1, 0.9], gradients_at=None, separate_panels=False, gradients_from_iobj=None,
           minmax=(RANGE_MIN,RANGE_MAX)  )

        pre_relaxation_verts = verts.copy()
        for i in range(VERTEX_RELAXATION_ITERATIONS_COUNT):
            VERTEX_RELAXATION_ADD_NOISE = False
            if VERTEX_RELAXATION_ADD_NOISE:
                verts = verts + (np.random.rand(verts.shape[0], verts.shape[1])*2.-1.)/2.* (STEPSIZE/8.)
            #set_trace()
            verts, facets_not_used, centroids = process2_vertex_resampling_relaxation(verts, facets, iobj, c=0*2.0)
            assert not np.any(np.isnan(verts.ravel()))  # fails
            print("Vertex relaxation applied.");sys.stdout.flush()
            if mesh_correction:
                check_degenerate_faces(verts, facets_not_used, "assert")
            #verts, facets_not_used, any_mesh_correction = check_degenerate_faces(verts, facets_not_used, "fix")
            #assert not np.any(np.isnan(verts.ravel()))  # fails
            #assert not any_mesh_correction
            #if any_mesh_correction:
            #    print("mesh correction needed")
            #    exit()

        #visualise_distance_histogram(pre_relaxation_verts, verts, facets)

        if VERTEX_RELAXATION_ITERATIONS_COUNT == 0:
            centroids = compute_centroids(verts, facets)

        if mesh_correction:
            assert len(failure_pairs) == 0, "weighted resampling did not work for some faces"
        if False:
            # list of pairs that have zero distance in weighted resampling.
            fpna = np.asarray(failure_pairs, dtype=np.int64).ravel()
            coords = verts[facets[fpna, :]]
            #quick_vis(verts, facets, fpna)
            #quick_vis(old_verts, old_facets, fpna)
            print coords

        point_collector.reset()


        # *******************************************************
        # PROJECTION
        # *******************************************************
        bigprint("projection")

        preprojection_vf = (verts, facets)

        #non-vectorized version
        #from ohtake_surface_projection import set_centers_on_surface__ohtake
        from ohtake_surface_projection_v2_5 import set_centers_on_surface__ohtake
        #vectorized version
        from ohtake_surface_projection_v2_5 import set_centers_on_surface__ohtake_v3s

        from project_kdtree import set_centers_on_surface__kdtree_v1

        average_edge = compute_average_edge_length(verts, facets)

        if mesh_correction:
            check_degenerate_faces(verts, facets, "assert")
            verts, facets, any_mesh_correction = check_degenerate_faces(verts, facets, "fix")

        c3 = np.mean(verts[facets[:], :], axis=1)
        old_centroids = np.concatenate((c3, np.ones((c3.shape[0], 1))), axis=1)

        nones_map = old_centroids[:, 0]*0 > 100  # all False
        new_centroids = old_centroids.copy()
        new_centroids1 = new_centroids.copy()
        new_centroids2 = new_centroids.copy()
        pre_proj_centroids = new_centroids.copy()

        with Timer() as t1:
            print
            print "1"*100; sys.stdout.flush()
            #set_centers_on_surface__ohtake(iobj, new_centroids1, average_edge*1., nones_map)
            #pass
            print
            print "---"*100; sys.stdout.flush()
        with Timer() as t2:
            print
            print "2"*100; sys.stdout.flush()
            #global _vs
            #global _fs
            #_vs, _fs = verts, facets

            facet_areas, facet_normals = compute_triangle_areas(verts, facets, return_normals=True)
            n = np.linalg.norm(facet_normals[np.logical_not(np.isnan(facet_areas)), :], axis=1)
            bads = np.logical_or(np.isnan(facet_areas), np.abs(facet_areas - 0.) < 0.00000001)
            assert np.allclose(np.linalg.norm(facet_normals[np.logical_not(bads), :], axis=1), 1.)  #is not zero
            facet_normals[bads, :] = 1./np.sqrt(3.)  # not tested
            del facet_areas
            assert facet_normals.shape[1] == 3
            assert np.allclose(np.linalg.norm(facet_normals, axis=1), 1.)

            print "average_edge", average_edge, STEPSIZE  # 0.089840676527 0.1
            #z12 =
            set_centers_on_surface__ohtake_v3s(iobj, new_centroids2, average_edge*1., nones_map, mesh_normals=facet_normals.copy())
                #debug_vf=(verts, facets))
            #new_centroids is the output
            print
            print "---"*100; sys.stdout.flush()
        print("Projection two methods: done within ", t1.interval, t2.interval, "RATIO =", t1.interval/t2.interval)

        with Timer() as t3:
            #set_centers_on_surface__kdtree_v1()
            pass



        #print (" Ratio= "+str(t1.interval/t2.interval))*30   # Linter's bug

        #9 times faster
        #11 times faster: 7.12189273357336, 0.6463103363521201, 'RATIO =', 11.019308114071752)
        # "-O" => 3.90 times faster
        #7 times faster

        #assert np.all(new_centroids1 == new_centroids2)
        new_centroids = new_centroids2
        #exit()

        if False:
            visualise_distance_histogram(pre_proj_centroids, new_centroids, facets)

        if mesh_correction:
            check_degenerate_faces(verts, facets, "assert")
            verts, facets, any_mesh_correction = check_degenerate_faces(verts, facets, "fix")

        # *********************************
        # QEM
        # *********************************

        #faceslist_neighbours_of_vertex
        #Why is this method called twice?
        vertex_neighbour_facelist_dict = mesh_utils.make_neighbour_faces_of_vertex(facets)
        centroid_gradients = compute_centroid_gradients(new_centroids, iobj)
        #nv1  =
        verts_before_qem = verts
        new_verts_qem = \
            vertices_apply_qem3(verts, facets, new_centroids, vertex_neighbour_facelist_dict, centroid_gradients, old_centroids_debug=old_centroids)
        #verts = nv1
        #new_verts_qem = verts
        print "QEM done"

        #print collector
        #set_trace()

        """ The following codevisualises the points that are not projected correctly. """
        del centroids
        if VISUALISE_RELAXATION_STEPS:
            #THRESHOLD_zero_interval = 0.0001
            #zeros2 = np.abs(f2) <= THRESHOLD_zero_interval
            #zeros1 = np.abs(f1) <= THRESHOLD_zero_interval
            #zeros12 = np.logical_or(zeros1, zeros2)  # output
            #todo:

            #Visualising the centroid points that the peojection has filed on them
            f_c = iobj.implicitFunction(new_centroids)
            #nzeros_c = np.nonzero(np.abs(f_c) <= 0.0001)[0]  # THRESHOLD_zero_interval
            nzeros_c = np.nonzero(np.abs(f_c) > 0.0001)[0]
            print "nonzero indices:",
            print "count:",nzeros_c.shape
            #print nzeros_c

            DONT_CUT = False
            sides_xyz = new_verts_qem[facets, :]
            b1 = sides_xyz[:, :, 0] < 0
            b2 = sides_xyz[:, :, 2] > -(1. - 0.1)  # -(1.15-0.1)
            b = np.logical_and(b1, b2)
            if DONT_CUT:
                b = b * 0 + 1
            cut_through = np.nonzero(np.sum(b, axis=1) >= 1)[0]
            print cut_through.shape, facets.shape
            #c3 = new_centroids[z12, :3]
            c3 = new_centroids[nzeros_c, :3]
            display_simple_using_mayavi_2([(new_verts_qem, facets[cut_through, :]), (preprojection_vf[0], preprojection_vf[1][cut_through, :])],
                       mayavi_wireframe=[False, True], opacity=[0.4/2.*2, 0.3],
                       gradients_at=c3,
                       separate_panels=False,
                       gradients_from_iobj=iobj,
                       minmax=(RANGE_MIN/100., RANGE_MAX/100.),
                       #add_noise=[0, 0], noise_added_before_broadcast=True,
                       #pointcloud_list=[new_centroids[z12, :]], pointsizes=[0.02], #pointcloud_list=[point_collector.get_as_array()], pointsizes=[0.01],
                       pointcloud_list=[new_centroids[nzeros_c, :]], pointsizes=[0.02], #pointcloud_list=[point_collector.get_as_array()], pointsizes=[0.01],
                       #labels=(new_centroids, z12), grad_arrow_len=0.2/2.)
                       labels=(new_centroids, nzeros_c), grad_arrow_len=average_edge*1. ,
            # Here you can easily visualise either the effect of projection of QEM shown as arrows.
                       #fromto=(pre_proj_centroids, new_centroids),
                       #fromto=(verts_before_qem, new_verts_qem),
                       fromto=[(pre_proj_centroids[:, :3], new_centroids[:, :3]), (verts_before_qem, new_verts_qem)],
                       )  #


        #no subdivision for now
        if SUBDIVISION_ITERATIONS_COUNT == 0:
            print "immediately after QEM *(without subdivision)"
            assert pre_relaxation_verts.shape == verts.shape
            if False:
             display_simple_using_mayavi_2( [(verts, facets), (new_verts_qem, facets)],
               pointcloud_list=[point_collector.get_as_array()], pointsizes=[0.01],
               mayavi_wireframe=[False, False], opacity=[0, 1],
               gradients_at=None, separate_panels=False, gradients_from_iobj=None,
               minmax=(RANGE_MIN,RANGE_MAX),
               add_noise=[0,0], noise_added_before_broadcast=True,
               fromto=(pre_relaxation_verts, verts) )#fromto=(verts, new_verts_qem)  )
        else:
            print "immediately after QEM *"
            ui = np.unique(facets.ravel())
            #up = find_unusual_points(verts, ui)
            up_i_10 = find_unusual_points(verts[ui, :])
            up = ui[up_i_10]
            #assert np.all(up == up_)

            print "u_p=", up
            if False:
             display_simple_using_mayavi_2([(verts, facets), (new_verts_qem, facets)],
               pointcloud_list=[new_verts_qem[up, :]], pointsizes=[0.02], #pointcloud_list=[point_collector.get_as_array()], pointsizes=[0.01],
               mayavi_wireframe=[False, True], opacity=[0.4, 1],
               gradients_at=None, separate_panels=False, gradients_from_iobj=None,
               minmax=(RANGE_MIN, RANGE_MAX),
               add_noise=[0, 0], noise_added_before_broadcast=True,
               labels=(new_verts_qem, up) )
               # labels is a very powerful method for visualising.

        #exit()


        if mesh_correction:
            any_mesh_correction = check_degenerate_faces(new_verts_qem, facets, "dontfix")
            if not any_mesh_correction:
                print "Mesh health:", any_mesh_correction

        #print "*"*500
        strict_about_mesh = True

        if mesh_correction:
          maxcount = 5
          if strict_about_mesh:
            while any_mesh_correction:
                new_verts_qem, facets, any_mesh_correction = check_degenerate_faces(new_verts_qem, facets, "fix")
                any_mesh_correction = check_degenerate_faces(new_verts_qem, facets, "dontfix")
                print "fixed"
                global still_fixed
                still_fixed = True
                # infinite loop
                maxcount -= 1
                if maxcount <0:
                    break
            if not take_it_easy:
                check_faces(facets)
            print "GOOd "*10


        alpha = 0.
        new_verts_qem_alpha = (new_verts_qem * alpha + verts * (1-alpha))

        if False:
            chosen_facet_indices = np.array(total_subdivided_facets, dtype=int)
            #centroids2, new_centroids2 = old_centroids[chosen_facet_indices], new_centroids[chosen_facet_indices]
            # move the following code into subdivide_multiple_facets() (?)
            if chosen_facet_indices.size == 0:
                chosen_subset_of_facets = np.zeros((0,), dtype=np.int64)
            else:
                chosen_subset_of_facets = facets[chosen_facet_indices, :]


        if False:  # if not take_it_easy:
            check_degenerate_faces(new_verts_qem_alpha, facets, "assert")
            check_degenerate_faces(new_verts_qem, facets, "assert")  # has degenerate face

        #(verts, facets) = (new_verts_qem, facets)


        if False:
            ifnoisy = 0
            #red balls
            highlighted_vertices = np.array([], dtype=np.int)  # np.arange(100, 200)
            hv = new_verts_qem[highlighted_vertices, :]
            display_simple_using_mayavi_2( [(new_verts_qem_alpha, facets), (new_verts_qem, facets), (new_verts_qem, facets), ],
               pointcloud_list=[],
               mayavi_wireframe=[False, True, True,], opacity=[0.2, 1, 0.3], gradients_at=None, separate_panels=False, gradients_from_iobj=None,
               minmax=(RANGE_MIN,RANGE_MAX),
               add_noise=[0.05*ifnoisy, 0, 0.05*ifnoisy], noise_added_before_broadcast=True  )
            #display_simple_using_mayavi_2( [(new_verts_qem_alpha, facets),(new_verts_qem, facets), ],
            #   pointcloud_list=[ hv ], pointcloud_opacity=0.2,
            #   mayavi_wireframe=[False, True], opacity=[0.2, 1, 0.9], gradients_at=None, separate_panels=False, gradients_from_iobj=None,
            #   minmax=(RANGE_MIN,RANGE_MAX)  )


            display_simple_using_mayavi_2( [(new_verts_qem_alpha, facets), (new_verts_qem, facets), (new_verts_qem, facets), ],
               pointcloud_list=[],
               mayavi_wireframe=[False, False, True,], opacity=[0.2, 1, 0.3], gradients_at=None, separate_panels=False, gradients_from_iobj=None,
               minmax=(RANGE_MIN,RANGE_MAX),
               add_noise=[0.05*ifnoisy, 0, 0.05*ifnoisy], noise_added_before_broadcast=True  )


        #]iobj = make_example_vectorized("ell_example1")

        verts = new_verts_qem
        #facets = facets
        (verts, facets) = (new_verts_qem, facets)

        # END of QEM


        # *******************************************************
        # SUBDIVISION
        # *******************************************************

        #input: verts, facets
        #settings: curvature_epsilon
        #output:verts, facers

        pre_subdiv_vf = (verts, facets)
        total_subdivided_facets = []
        for i in range(SUBDIVISION_ITERATIONS_COUNT):
            #set_trace()

            print "subdivision:"
            verts, facets = do_subdivision(verts, facets, iobj, curvature_epsilon)
            global trace_subdivided_facets  # third implicit output
            verts4_subdivided = verts  # ??
            facets3_subdivided = facets

            if mesh_correction:

                check_degenerate_faces(verts4_subdivided, facets3_subdivided, "assert")
                #verts4_subdivided, facets3_subdivided, any_mesh_correction = check_degenerate_faces(verts4_subdivided, facets3_subdivided, "fix")

                total_subdivided_facets += trace_subdivided_facets  # old face indices remain valid

                #for i in range(VERTEX_RELAXATION_ITERATIONS_COUNT):
                #    #print "i", "="*10, i
                #    verts, facets_not_used, centroids = process2_vertex_resampling_relaxation(verts, facets, iobj)
                #    print("Vertex relaxation2 applied again.");sys.stdout.flush()
                #    check_degenerate_faces(verts, facets_not_used, "assert")
                #    verts, facets_not_used, any_mesh_correction = check_degenerate_faces(verts, facets_not_used, "fix")
        print "subdivision done."

        #verts, facets


        for use_wireframe in [True, False]:
            ifnoisy = 0
            display_simple_using_mayavi_2( [(verts, facets), (verts, facets), (pre_subdiv_vf[0], pre_subdiv_vf[1]), ],
               pointcloud_list=[],
               mayavi_wireframe=[False, use_wireframe, True,], opacity=[0.2, 1, 0.3], gradients_at=None, separate_panels=False, gradients_from_iobj=None,
               minmax=(RANGE_MIN,RANGE_MAX),
               add_noise=[0.05*ifnoisy, 0, 0.05*ifnoisy], noise_added_before_broadcast=True  )



def compute_average_edge_length(verts, faces):
    nfaces = faces.shape[0]
    expand = verts[faces, :]
    assert expand.shape == (nfaces, 3, 3)
    assert expand[:, 2, :].shape == (nfaces, 3)
    ea_sum = 0.
    for i in range(3):
        i1 = i
        i2 = (i+1) % 3
        e1 = np.linalg.norm(expand[:, i1, :] - expand[:, i2, :], axis=1)  # bug fixed! axis=1 was missing.
        assert e1.shape == (nfaces,)
        ea_sum += np.mean(e1)

    #set_trace()
    return ea_sum / 3.


collector = []

#worst_30 = \
#    np.array([  858,   238,   839,   171,   799,   199,   211,  1068,
#        1135,   998,  1276,   812,   243,   794,   643,   296,
#         957,   529,   294,   453,   960,   531,  1279,  1607,
#        1066,   303,   809,   413,   245,   700])
#
#worst_30 = worst_30[-1:]
worst_30 = np.array([], dtype=int)

def vertices_apply_qem3(verts, facets, centroids, vertex_neighbour_facelist_dict, centroid_gradients, old_centroids_debug):
    """
    result_verts_ranks (is not returned) is -1 if the verts has no neighbour (i.e. not included in facets, i.e. not participating in the mesh."""
    #based on quadratic_optimise_vertices(self, alpha=1.0)
    assert not centroids is None
    assert not vertex_neighbour_facelist_dict is None
    assert not centroid_gradients is None

    #global collector
    #alpha = 1.0
    nvert = verts.shape[0]
    #Check if all vertices are included.
    #assert nvert == len(vertex_neighbour_facelist_dict)  # can be different.
    #if not np.subset(vertex_neighbour_facelist_dict, range(verts.shape[0])):
    #        set_trace()
    #assert np.subset(vertex_neighbour_facelist_dict, range(verts.shape[0]))
    assert set(vertex_neighbour_facelist_dict.keys()).issubset(set(range(verts.shape[0])))

    result_verts_ranks = np.zeros((nvert,), dtype=int)
    assert verts.shape == (nvert, 3)
    new_verts = np.zeros((nvert, 3))

    for vertex_id in range(nvert):
        vi = vertex_id
        if vertex_id not in vertex_neighbour_facelist_dict:
            #print vertex_id, "not in the neighbours dict"
            #exit()
            new_verts[vi, 0:3] = verts[vi, 0:3]
            result_verts_ranks[vi] = -1  # never returned
            continue   # leave it alone
        neighbours_facelist = vertex_neighbour_facelist_dict[vertex_id]
        neighbours_faces = np.array(neighbours_facelist, dtype=np.int64)
        #qem_origin = np.zeros((3, 1))
        qem_origin = verts[vertex_id, :].reshape(3, 1)
        A, b = get_A_b(vi, neighbours_faces, centroids, centroid_gradients, qem_origin)
        #print "Ax+b=0: A,b:", A, b

        ###
        #A, b = self.get_A_b(vi)

        u, s, v = np.linalg.svd(A)
        assert np.allclose(A, np.dot(u, np.dot(np.diag(s), v)))
        #print(s)  # [  1.48148148e+01   1.67928330e-15   1.01592270e-50]
        assert s[0] == np.max(s)
        #print( s / s[0] )  # [  1.00000000e+00   1.13351623e-16   6.85747820e-52]

        #tau = 1000.
        tau = 680.
        #tau = (26.077)**2
        #tau = (0.11 * 237.)**2

        s[s / s[0] < 1.0/tau] = 0
        #print(s , s[0] , tau)
        rank = np.sum(s / s[0] > 1.0/tau)
        #if rank==1:
        # Threshold_minimum_sigma
        #      rank = np.sum(s / s[0] > Threshold_minimum_sigma)
        # assert rank <= 1

        #print(s)
        #print("rank = ", rank)

        #rank will never be 0: s[0]/s[0] is always 1, even when s[0] is too small.
        #assert s[0] > 0.000001

        if not  s[0] > 0.000001:
            print("Warning! sigma_1 == 0" )
            print(s)
            print("A", A)

            #not tested
            result_verts_ranks[vi] = 0
            new_verts[vi, 0:3] = new_x[:, 0]

        assert np.all(s[:rank]/s[0] >= 1.0/tau)

        x = verts[vi, 0:3, np.newaxis] - qem_origin
        assert x.shape == (3, 1)
        #print "x_OLD", x
        global x_old
        x_old = x
        y = np.dot(v, x).copy()
        #x,y are default x,y elements (in case of degeneracy)

        utb = np.dot(-np.transpose(u), b)
        #print("rank", rank, "1/tau=", 1./tau)
        #print s
        for i in range(rank):
            #print(np.dot(-np.transpose(u), b), "scalar")
            assert np.dot(-np.transpose(u), b).shape == (3,1)
            #print s[i] , 1.0/tau
            #assert s[i] >= 1.0/tau #fails when s[0] is small
            assert s[i]/s[0] >= 1.0/tau
            y[i] = utb[i] / s[i]
            #print 1./ s[i], "*"
            #collector.append((1./ s[i], vertex_id))
        if vertex_id in worst_30:
            #set_trace()
            #print
            #print
            pass
        assert qem_origin.shape == (3, 1)
        new_x = np.dot(np.transpose(v), y) + qem_origin
        #print(x.ravel(), " -> ", new_x.ravel())
        #print("    delta=", (new_x - x).ravel())

        new_verts[vi, 0:3] = new_x[:, 0]
        #self.new_verts[vi,3] = 1

        displacement = verts[vi,0:3]-new_verts[vi, 0:3]
        if np.linalg.norm(displacement) > STEPSIZE*3:
            #vi==310: norm is > 10.
            #set_trace()
            collector.append((vi, np.linalg.norm(displacement)))
            pass


        assert x.shape == (3, 1)
        # Apply alpha
        #new_verts[vi, 0:3] = new_x[:, 0] * alpha + x[:, 0] * (1.0-alpha)
        new_verts[vi, 0:3] = new_x[:, 0]

        if not np.all(np.abs(utb.ravel()[rank:] ) < 0.0001):
            #print("s", s.ravel()/s[0], "   utb", utb.ravel()/s[0])
            pass
        result_verts_ranks[vi] = rank



        #254,
        if vi==179: #254: #179: #81: #310:
            print "here"
            def q():
                print "x", x_old
                print "new_x", new_x - qem_origin
                set_trace()
                print "neighbours_faces:", neighbours_faces
                for x in neighbours_faces:
                    #873
                    print "face", x, ":", facets[x, :], "centroid ", centroids[x, :]
                fa = facets[neighbours_faces, :]
                pc = centroids[neighbours_faces, :]
                pc_old = old_centroids_debug[neighbours_faces, :]
                # old_centroids_debug
                global giobj

                display_simple_using_mayavi_2(
                    [(verts, fa), (verts, fa), (verts, facets),],
                    pointcloud_list=[ pc, pc_old, new_verts[np.newaxis, vi, 0:3] ], pointcloud_opacity=0.2, pointsizes=[0.2*0.2, 0.2*0.1, 0.07]*2,
                    mayavi_wireframe=[False, True, False], opacity=[0.2, 1, 0.1],
                    gradients_at=pc[:, :3], separate_panels=False, gradients_from_iobj=giobj,
                    add_noise=[0, 0, 0.], noise_added_before_broadcast=True, random_interior_point=False )
                pass
            if False:
                q()
            """
            def
                display_simple_using_mayavi_2( [(verts, facets), (verts, facets[face_idx, :]), (verts, facets[face_idx, :])],
                   pointcloud_list=[ ], pointcloud_opacity=0.2,
                   mayavi_wireframe=[True, True, True], opacity=[0.1, 1, 0.2],
                   gradients_at=None, separate_panels=False, gradients_from_iobj=None, add_noise=[0, 0, 0.01],
                   noise_added_before_broadcast=True
                   )
            """
            pass

        #exit()
    if result_verts_ranks.size > 0:
        print("max rank = ", np.max(result_verts_ranks))
        print("min rank = ", np.min(result_verts_ranks))
        if not np.min(result_verts_ranks) >= 1:
            print("Warning: assertion: np.min(result_verts_ranks) >= 1 failed." )
    else:
        print("max,min rank = []")


    print np.array(collector)
    #set_trace()
    if False:
        assert result_verts_ranks.size == 0 or np.min(result_verts_ranks) >= 1
    return new_verts

def get_A_b(vertex_id, neighbours_faces, centroids, centroid_gradients, qem_origin):
    #refactor: faces_array contains faces.  -> neightbours_facelist or neightbours_faces
    #neighbours_facelist = self.vertex_neighbour_facelist_dict[vertex_id]
    #neighbours_faces = np.array(neighbours_facelist, dtype=int)
    #neighbours_faces = neighbours_faces
    assert qem_origin.shape == (3, 1)
    center_array = centroids[neighbours_faces, :]

    #note some centers may not be projected successfully in the previous step
    not_projected_successfully = np.isnan(center_array[:].ravel())
    if np.any(not_projected_successfully):
        pass

    normals = centroid_gradients[neighbours_faces, :]  #why do we have repeats??
    #note : not normalised. But it works.

    norms = np.linalg.norm(normals, ord=2, axis=1)
    #can be either 0, 1 or Nan
    if np.any(norms < 0.000001):  #can be exactly 0.0
        print("Error: bad normal", normals)

    TH_N = 0.0000001  # 0.000001 = I millions of millimeter = 1 nanometer
    #can be 0,0,0, inf, nonsharp, degenerate, ...
    degenerate_normals = np.logical_or(np.isnan( np.sum(normals, axis=1)), norms < TH_N )
    #simpler: degenerate_normals = np.logical_or(np.isnan(norms), norms < 0.0000001 )
    #todo:



    #print(normals)
    assert not np.any(np.isnan(normals) )
    assert not np.any(np.isinf(normals) )

    #normals = normalize_vector4_vectorized( normals ) #todo: move it to evaluate_centroid_gradients or self.centroid_normals

    #print("normals", normals) # can be all 0,0,0

    #qem_origin = np.zeros((3, 1))

    assert normals.shape[1] == 4
    #normals = normals   # (4)x4
    #grad = Ax+b
    A = np.zeros((3, 3))
    b = np.zeros((3, 1))
    #assert len(center_array) == len(normals)
    assert normals.shape == center_array.shape
    for i in range(normals.shape[0]):
        n_i = normals[i, 0:3, np.newaxis]
        assert n_i.shape == (3, 1)
        nnt = np.dot(n_i, np.transpose(n_i))
        assert nnt.shape == (3, 3)
        A += nnt
        #It is correct if A contains equal rows. In this case, we have faces that are parallel or on the same plane (e.g. on the same side of a cube)
        p_i = center_array[i, 0:3, np.newaxis]
        assert p_i.shape == (3, 1)
        b += -np.dot(nnt, p_i - qem_origin)
        if vertex_id==310:
            #set_trace()
            #print "loop end"
            #print "loop end"
            #print "loop end"
            pass

        # IN PROGRESS

    assert not np.any(np.isnan(A.ravel()))
    assert not np.any(np.isnan(b.ravel()))

    return A, b


def fix_windws_control_C():
    if sys.platform != 'win32':
        print "NON-WINDoWS"
        return  # no need
    print "YES WINDoWS"

    #http://stackoverflow.com/questions/15457786/ctrl-c-crashes-python-after-importing-scipy-stats
    import os
    import imp
    import ctypes
    import thread
    import win32api

    # Load the DLL manually to ensure its handler gets
    # set before our handler.
    basepath = imp.find_module('numpy')[1]
    loc = 'C:\\Anaconda3\\pkgs\\mkl-11.3.1-0\\Library\\bin'
    #ctypes.CDLL(os.path.join(basepath, 'core', 'libmmd.dll'))
    #ctypes.CDLL(os.path.join(basepath, 'core', 'libifcoremd.dll'))
    ctypes.CDLL(os.path.join(loc, 'libmmd.dll'))
    ctypes.CDLL(os.path.join(loc, 'libifcoremd.dll'))

    # Now set our handler for CTRL_C_EVENT. Other control event
    # types will chain to the next handler.
    def handler(dwCtrlType, hook_sigint=thread.interrupt_main):
        if dwCtrlType == 0: # CTRL_C_EVENT
            hook_sigint()
            print "CONTROL +C CAUGHT"
            return 1 # don't chain to the next handler
        return 0 # chain to the next handler

    win32api.SetConsoleCtrlHandler(handler, 1)


#@profile
def eval1(iobj, x):
    return iobj.implicitFunction(x)

#@profile
def eval2(iobj, x):
    y = np.zeros(x.shape[0])
    for i in range(x.shape[0]):
        y[i] = iobj.implicitFunction(x[i, np.newaxis, :])
    return y


#@profile
def compare_vectorised_speed():
    import sys
    print "hi"; sys.stdout.flush()
    iobj, RANGE_MIN, RANGE_MAX, STEPSIZE = make_bricks()

    x = np.random.rand(100000, 4); x[:,3] = 1

    print "A"; sys.stdout.flush()
    f = eval1(iobj, x)  # 244 x times faster!

    print "B"; sys.stdout.flush()
    g = eval2(iobj, x)

    print np.sum(f-g)
    sys.stdout.flush()



def cube_with_cylinders(SCALE):
    #Solene
    from basic_types import make_vector4
    import vectorized

    SCALE = 2.  # mm
    sz1 = 2.5

    radius = 0.5 * SCALE
    c_len = 2 * SCALE

    A = make_vector4(-c_len/2.0, 0, 0)
    #A = make_vector4(0, 0, c_len / 2.0)  # bug: aa is wrong
    w = make_vector4(1,0 , 0)
    w = w / np.linalg.norm(w[0:3]); w[3] = 1
    u = make_vector4(0, 1, 0)

    cyl = vectorized.SimpleCylinder(A, w, u, radius, radius, c_len)


    A2 = make_vector4(0, -c_len/2.0, 0)
        #A = make_vector4(0, 0, c_len / 2.0)  # bug: aa is wrong
    w2 = make_vector4(1,0 , 0)
    w2 = w / np.linalg.norm(w[0:3]); w[3] = 1
    u2 = make_vector4(0, 1, 0)

    cyl_2 = vectorized.SimpleCylinder(A2, u2, w2, radius, radius, c_len)
    cube = vectorized.UnitCube1(size=sz1)
    union = vectorized.CrispSubtract(cube, cyl_2)
    final_object = vectorized.CrispUnion(union,cyl)

    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-3, +5, 0.2)
    return final_object, RANGE_MIN, RANGE_MAX, STEPSIZE


def find_unusual_points(verts, k=3):
    """ k is the order."""
    indices = np.arange(verts.shape[0])
    """ Find the points that their distance from the rest of the points is large.
    The distance is the minimum distance from the rest of the points.
    This test is done only on indices"""
    from sklearn.neighbors import NearestNeighbors
    assert verts.shape == (verts.shape[0], 3)
    nbrs = NearestNeighbors(n_neighbors=k, algorithm='ball_tree').fit(verts)
    distances1, ind1 = nbrs.kneighbors(verts[indices, :])
    #print distances1, ind1
    #distances1[:, 1]  #[:,0] are all zero
    d1 = distances1[:, k-1]
    i1 = ind1[:, k-1]
    ten_most_far1 = np.argsort(d1)[-10:]
    #print ten_most_far1
    print indices[ten_most_far1]

    from sklearn.neighbors import KDTree
    kdt = KDTree(verts, leaf_size=30, metric='euclidean')
    #ind2= kdt.query(verts[indices, :], k=2, return_distance=False)
    dist2, ind2 = kdt.query(verts[indices, :], k=k, return_distance=True)
    dist2, ind2 = dist2[:, k-1], ind2[:, k-1]
    ten_most_far2 = np.argsort(dist2)[-10:]
    #print ten_most_far2
    print indices[ten_most_far2]
    #set_trace()
    #print "returning"
    assert np.all(indices[ten_most_far1] == indices[ten_most_far2])

    from sklearn.neighbors import KDTree
    kdt = KDTree(verts, leaf_size=30, metric='euclidean')
    dist3, ind3 = kdt.query(verts[indices, :], k=k, return_distance=True)
    dist3, ind3 = dist3[:, k-1], ind3[:, k-1]
    ten_most_far3 = np.argsort(dist3)[-10:]

    return indices[ten_most_far3]

if __name__ == '__main__':
    fix_windws_control_C()


    #test_example_meshes()
    #print "good"
    #exit()

    demo_everything(options={})

    #compare_vectorised_speed()


#from stl_tests import display_simple_using_mayavi_vf1
from ohtake_surface_projection import display_simple_using_mayavi_


def fix_degenerate_Faces():
    pass




def make_test_vf_1():
    faces = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
    verts = np.random.randn((4, 3))
    return verts, faces

def make_test_vf_2():
    faces = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
    verts = np.random.randn((8, 3))
    return verts, faces
