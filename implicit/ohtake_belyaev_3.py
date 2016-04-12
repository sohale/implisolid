#note: 3 points should NOT be projected into ONE point. DEFINITELY. => dont remove the faces. Only two vertices at a time.
from vtk_mc import vtk_mc
import sys
import numpy as np
from basic_types import check_vector4_vectorized
from basic_types import normalize_vector4_vectorized

mesh_quality_settings = {
    "min_edge_len":  0.000001   # 0.01  # 0.001 microns   # 0.001, 0.007
    }

CHECK_PAIRED = False


class TroubledMesh(Exception):
    pass


def check_face_triplets(faces):
    # unique faces
    f3sides = faces.copy()
    f3sides.sort(axis=1)

    B = 100000
    BBB = np.array([[1, B, B*B]]).transpose()  # 3x1

    d = np.dot(f3sides, BBB)
    face_triplet_ids = d.ravel()
    del f3sides

    face_order = face_triplet_ids.argsort()
    face_triplet_ids.sort()
    #print face_triplet_ids[np.diff(face_triplet_ids)==0]
    # Check there is no repeated faces (with exact same vertex list):
    diff0 = (np.diff(face_triplet_ids) == 0)
    if np.sum(diff0) != 0:
        print np.sum(diff0)  # 74
        nonz = np.nonzero(diff0)[0]
        #print nonz
        # diff01: for print only: to print both sides (Elements) of each "diff==0"
        diff01 = diff0.copy()
        diff01[nonz+1] = True
        #print face_triplet_ids[diff0]
        #print face_triplet_ids[diff01]
        #

    nonz = np.nonzero(diff0)[0]
    bad_faces = nonz
    #but some repeated ones may remain
    #howver the original idx (of redundant faces) are = face_order[nonz]
    original_indices = face_order[bad_faces]

    #print "***", original_indices
    assert np.sum(diff0) == original_indices.size, "Number of redundant faces"
    #assert np.sum(diff0) == 0, "Repeated faces found"
    return original_indices


def check_faces(faces):
    print("------ check_faces(faces)")

    redundant_faces = check_face_triplets(faces)
    assert redundant_faces.size==0, "Repeated faces found"


    #unique edges
    f0 = faces[:, np.newaxis, 0:2]
    f1 = faces[:, np.newaxis, 1:3]
    f2 = faces[:, np.newaxis, np.array([0, 2])]
    #print faces[:,:]
    #print f2
    f0 = f0.copy(); f0.sort(axis=2)  # changes the order in faces!
    f1 = f1.copy(); f1.sort(axis=2)
    f2 = f2.copy(); f2.sort(axis=2)
    #print faces[:,:]
    #print f2
    #exit()

    fe3 = np.concatenate( (f0, f1, f2), axis=1 )  # shape==(:,3,2)
    #fe3 = np.swapaxes(fe3, 0, 1)
    print fe3.shape
    #print fe3.shape
    B = 100000
    BB = np.array([[1, B]]).transpose().ravel()  # 2x-
    edg = np.dot(fe3, BB)
    #print edg.shape  # fx3
    #print edg

    #Sort edges to detect repeated edges. Each edge should appear exactly twice.
    q = edg.ravel()
    sort_idx = q.argsort()
    assert sort_idx.shape == (faces.shape[0]*3,)
    print sort_idx
    q_unsorted = q.copy()
    #q=-q
    q.sort()
    #q=-q
    assert np.all(q_unsorted[sort_idx] == q)
    #print q.reshape( (-1, 2) ).transpose()
    print q[:10]
    print "diff=", np.diff(q)[:10]

    fe3_ravel = fe3.reshape( (np.prod(fe3.shape[0:2]), 2) )
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
    print "sort_idx=", sort_idx[:10]
    print "face, side=", face_idx, side_idx

    vert_idx_1 = side_idx
    vert_idx_2 = (side_idx + 1) % 3
    v1 = faces[face_idx, vert_idx_1]
    v2 = faces[face_idx, vert_idx_2]
    #print "V1, v2"
    #print v1
    #print v2
    #print
    v12 = np.concatenate((v1[:, np.newaxis], v2[:, np.newaxis]), axis=1)

    #print
    #print fe3[sort_idx[:10], 1, :]
    #print fe3_ravel[sort_idx[:10], :]

    #Now they are the same: v12 and fe3_ravel
    print v12
    print fe3_ravel[sort_idx[:10], :]

    #You can get the sorted inices here:
    #faces[face_idx, vert_idx_1]
    #faces[face_idx, vert_idx_2]
    #exit()

    #assert np.all(np.diff(q)[1::2] != 0)  # what about the very first one

    if not np.all(np.diff(q)[::2] == 0):
        print "q"
        #print q.reshape( (q.size, 1) )

        #for i in range(q.size):
        #    print q[i],
        #print

        #dd = np.diff(q)[::2]
        #for i in range(dd.size):
        #    print dd[i],
        #print

        i1 = (np.diff(q)[::2] != 0)
        #print q[::2][i1]
        #print q[1::2][i1]
        #print np.zip(q[::2][i1], q[1::2][i1])
        i00 = np.nonzero(i1)[0][0]
        print i00, "i1[0]"
        i0 = i00*2-4
        for i in range(i0, i0+20): #range(zp.size):
            print q[i],
        print
        #yes, some element is repeated 3 times
        if False:
            exit()


            q0 = q[::2][i1][:, np.newaxis]
            q1 = q[1::2][i1][:, np.newaxis]
            print i1.shape
            print q[::2].shape
            zp = np.concatenate( (q0, q1), axis=1)
            #print zp.ravel()
            print zp.shape, "zp"


            del i1
            for i in range(20): #range(zp.size):
                print zp.ravel()[i],
            print

    if CHECK_PAIRED:
        assert np.all(np.diff(q)[::2] == 0), "Not all edges are paired"


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


#def fix_faces_3div2(faces):
#    from mesh_utils import make_neighbour_faces_of_vertex
#    neighbour_faces_of_vertex = make_neighbour_faces_of_vertex(facets)
#    faces_of_faces = build_faces_of_faces(facets)

def REMOVE_REPEATED_EDGES(faces):
    #see check_faces()

    f0 = faces[:, np.newaxis, 0:2]
    f1 = faces[:, np.newaxis, 1:3]
    f2 = faces[:, np.newaxis, np.array([0, 2])]
    f0 = f0.copy(); f0.sort(axis=2)  # changes the order in faces!
    f1 = f1.copy(); f1.sort(axis=2)
    f2 = f2.copy(); f2.sort(axis=2)
    fe3 = np.concatenate( (f0, f1, f2), axis=1)  # shape==(:,3,2)

    B = 100000
    BB = np.array([1, B])
    edg = np.dot(fe3, BB)
    print edg.shape  # fx3
    print edg

    #Sort edges to detect repeated edges. Each edge should appear exactly twice.
    #was q
    edg_sorted = edg.ravel()
    sort_idx = edg_sorted.argsort()
    assert sort_idx.shape == (faces.shape[0]*3,)
    print sort_idx
    q_unsorted = edg_sorted.copy()
    edg_sorted.sort()
    assert np.all(q_unsorted[sort_idx] == edg_sorted)
    print edg_sorted[:10]
    print "diff=", np.diff(edg_sorted)[:10]

    assert np.all(np.diff(edg_sorted)[0::2] == 0), "some edges are not exactly repeated once: not manifold."
    assert np.all(np.diff(edg_sorted)[1::2] != 0), "some edges are not exactly repeated once: not manifold."
    print np.diff(edg_sorted)
    print edg_sorted.shape, "edg_sorted.shape", np.array(edg_sorted.shape) / 3.
    print faces.shape
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
    print "sort_idx=", sort_idx[:10]
    print "face, side=", face_idx, side_idx

    vert_idx_1 = side_idx
    vert_idx_2 = (side_idx + 1) % 3
    v1 = faces[face_idx, vert_idx_1]
    v2 = faces[face_idx, vert_idx_2]

    v12 = np.concatenate((v1[:, np.newaxis], v2[:, np.newaxis]), axis=1)

    #Now they are the same: v12 and fe3_ravel
    print v12
    print fe3_ravel[sort_idx[:10], :]

    #You can get the sorted inices here:
    #faces[face_idx, vert_idx_1]
    #faces[face_idx, vert_idx_2]

    if not np.all(np.diff(edg_sorted)[::2] == 0):
        pass






def REMOVE_REPEATED_FACES_NEVER_NECESSARY(faces):
    """ Remove faces that are exactly repeated """
    print "ok"

    # based on check_face_triplets()
    # unique faces
    f3sides = faces.copy()
    f3sides.sort(axis=1)

    B = 100000  # better: B = max of vertex index = 10**(floor(log10(verts.shape[0]+111))+1); assert ...
    BBB = np.array([1, B, B*B])

    face_triplet_ids = np.dot(f3sides, BBB)
    #face_triplet_ids = face_triplet_ids
    del f3sides

    #print face_triplet_ids
    #print face_triplet_ids.shape
    #exit()

    face_order = face_triplet_ids.argsort()
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
        print redunds + 1
        exit()
    print "no problem"
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
    print "going to kill", number_of_faces_to_kill, "degenerate faces"

    faces_to_kill.sort(axis=1)
    face_map10 = np.zeros((2, 0), dtype=int)
    for ei in [1, 2]:
        from_idx = faces_to_kill[:, ei]
        to_idx = faces_to_kill[:, 0]
        print to_idx, to_idx.shape, "*"
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

    print faces.shape

    print face_map10.shape, "******************"  # too many
    print "???", faces.shape, faces.shape[0]*2  #WRONG  ****************************************

    fmap_ = np.arange(faces.shape[0])
    fmap_[face_map10[1, :]] = face_map10[0, :]

    for i in range(6):
        fmap_ = fmap_[fmap_]
        dosnt_need_change = np.all(fmap_ == fmap_[fmap_])  # 4 times False!
        #print dosnt_need_change
    assert dosnt_need_change

    print fmap_
    print fmap_[fmap_]
    assert np.all(fmap_ == fmap_[fmap_])
    print "number of mapped vertices", np.sum(fmap_ != np.arange(faces.shape[0]))
    dying_vertices = np.arange(faces.shape[0])[fmap_ != np.arange(faces.shape[0])]
    dying_vertices_bool = fmap_ != np.arange(faces.shape[0])
    affected_faces = np.any(dying_vertices_bool[faces], axis=1 )
    print "affected faces", np.sum(affected_faces)
    print "full faces", np.sum( np.all(dying_vertices_bool[faces], axis=1 ) )
    #number of faces that all of their elements are mapped onto something

    print "2 faces", np.sum( np.sum(dying_vertices_bool[faces], axis=1 )>=2 ), "faces that have 2 vertices"
    for s in range(5):
        print "faces with %d taken vertices:"%(s,), np.sum( np.sum(dying_vertices_bool[faces], axis=1 )==s )
    #print np.nonzero(affected_faces)

    face__projected_vertices = fmap_[faces]

    v10 = face__projected_vertices[:, 1] == face__projected_vertices[:, 0]
    v20 = face__projected_vertices[:, 2] == face__projected_vertices[:, 0]
    v21 = face__projected_vertices[:, 1] == face__projected_vertices[:, 2]
    v1020 = np.logical_and(v10, v20)
    print np.sum(v1020), "faces with three (now) identical vertices"
    # Why is it 81 and 80 = number_of_faces_to_kill. (fixme)
    #print faces[v1020, :]
    #p0 = 131  # to what index is p0 projected?
    #print p0, "->",face_map10[0, :][face_map10[1, :]==p0]
    #exit()

    faces_to_annihilate = np.logical_or(v10, v20, v21)  #129 ?!
    print "faces_to_annihilate", np.sum(faces_to_annihilate)

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


def quick_vis(verts, facets, face_idx):
    print facets.shape
    print facets[face_idx, :].shape
    display_simple_using_mayavi_2( [(verts, facets), (verts, facets[face_idx, :]), (verts, facets[face_idx, :])],
       pointcloud_list=[ ], pointcloud_opacity=0.2,
       mayavi_wireframe=[False, True, True], opacity=[0.1, 1, 0.2],
       gradients_at=None, separate=False, gradients_from_iobj=None, add_noise=[0, 0, 0.01],
       )

#todo: review this code: , separate it, u-test it.
def remove_vertices_and_faces(verts, faces, nil_areas_whichfaces, map12):
    """ nil_areas_whichfaces: indices of faces that need to be removed.
    map12: the vertices that have to be replaced (remove map12[:,1] -> replace as: map12[:,0])."""

    assert faces.shape[0] > 0
    assert verts.shape[0] > 0

    assert map12.shape[0] == 2
    print "-"
    u =unused_vertices_slow(verts, faces)
    print "unused vertices: ", u
    for x in u:
        print np.sum(faces.ravel()==x),
    print ""

    if False:
        #Note: there are cases where the edges are not equal (0.0000001) but the area is small (0.001).
        print faces[181, :], faces[219, :], faces[952, :] # 181  219  952  993 1281 1295
        f1 = faces[ np.array([181,219,952]), :]
        print "areas = ", compute_triangle_areas( verts, f1 , return_normals=False, AREA_DEGENERACY_THRESHOLD=None)
        print verts[f1,:]

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

    print "*new_faces", new_faces

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
        print np.intersect1d(from_, new_faces.ravel())
        assert np.intersect1d(from_, new_faces.ravel()).size == 0, "Confirm no more vertex clustering necessary, hence the map is complete, does not need more repeats, and is not cylic.."


    assert new_faces.shape[0] > 0

    print "map_", map_
    print "new_faces", new_faces
    new_verts, new_faces = delete_unused_vertices(verts, new_faces)

    assert new_faces.shape[0] > 0

    #killed_whichvertices = np.zeros((verts.shape[0],), dtype=np.bool) + True
    #killed_whichvertices[new_faces.ravel()] = False
    #print "unused count", np.sum(killed_whichvertices)
    print "unused vertices: ",  unused_vertices_slow(verts, new_faces)
    #new_verts = verts


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
    print new_faces2  #most of them are zero
    #exit()
    quick_vis(verts, new_faces, nil_areas_whichfaces)

    assert new_faces2.shape[0] > 0

    assert new_faces2.shape[0] > 0


    #checking if the nil_areas_whichfaces matches the faces that are now projected
    eq1 = new_faces2[:, 0] == new_faces2[:, 1]
    eq2 = new_faces2[:, 0] == new_faces2[:, 2]
    eq12 = np.logical_and(eq1, eq2)  # The triangles which all vertices have the same index number.
    eq1or2 = np.logical_or(eq1, eq2)
    print eq12, " sum1=",np.sum(eq12)  # 81
    print " sum1or2=",np.sum(eq1or2)
    print nil_areas_whichfaces, " sum2=", np.sum(nil_areas_whichfaces)  # 80
    #81 and 80 consequently
    #print np.nonzero( nil_areas_whichfaces != eq12)

    print "XOR", np.sum( eq12 != nil_areas_whichfaces)
    xor_ = eq12 != nil_areas_whichfaces
    #print new_faces2[xor_, :], np.nonzero(xor_)   # [290,290,290], 892

    print np.sum( np.logical_and(np.logical_not(eq12), nil_areas_whichfaces))
    assert np.sum( np.logical_and(np.logical_not(eq12), nil_areas_whichfaces)) == 0, "assert nil_areas_whichfaces is-subset eq12"
    def x():
        n = np.logical_not; a = eq12; b = nil_areas_whichfaces; aNd = np.logical_and
        print "a & b", np.sum(aNd(a, b))  # 80
        print "a & ~b", np.sum(aNd(a, n(b)))  # 1   Faces that are not (area=0), but their vertices are not identical. But how?
        print "~a & b", np.sum(aNd(n(a), b)) # 0
    x()

    a_n_b = np.nonzero(np.sum(np.logical_and(eq12, np.logical_not(nil_areas_whichfaces))))[0]
    print a_n_b
    quick_vis(verts, faces, a_n_b)


    #print new_faces2[eq12, :]
    #print new_faces2[nil_areas_whichfaces, :] #nice
    #81 and 80 consequently


    #print "DELETING: %d out of %d"%(-np.sum(eq12), new_faces2.shape[0]), "="*90
    #import time; time.sleep(2)

    #Now: faces with repeats in their vertices should be cut away.
    #It is also the concern of: map_vertices_of_nil_faces
    #new_faces = new_faces2
    #new_faces = new_faces2[ eq12 ]  # kill! -> oops.  wrong ones.
    if False:
        new_faces = new_faces2[ np.logical_not(eq12) ]  # kill!

    new_faces = new_faces2[ np.logical_not(eq1or2) ]

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
    eq12 = np.logical_and(eq1, eq2)

    #print "DELETING: %d out of %d"%(np.sum(eq12), new_faces.shape[0]), "-"*290
    #print new_faces[eq12, :]
    #import time; time.sleep(2)

    new_faces = new_faces[np.logical_not(eq12), :]

    assert new_faces.shape[0] > 0

    print new_faces
    print "-"*200
    check_faces(new_faces)
    #passed!

    return new_verts, new_faces


# **********************************
def check_degenerate_faces(verts, facets, fix_mode="dontfix"):
    # todo: also check facets in-itself.
    return verts, facets, False

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

    vertices_to_combine = []
    faces_to_remove = []

    nil_edgelen_whichfaces = np.zeros((facets.shape[0],), dtype=np.bool)
    nil_areas_whichfaces = np.zeros((facets.shape[0],), dtype=np.bool)

    v12 = np.zeros((2, 0), dtype=int)

    any_zero_edge = False
    el = [e1, e2, e3]
    # visualise_edge_distribution(el)  # keep this line
    for i in [0, 1, 2]:
            #if if_assert:  # ?????????????????????????
            j1 = i; j2 = (i-1+3) % 3
            tv3 = verts[facets[:, j1], :] - verts[facets[:, j2], :]
            #tv3 = verts[facets[:, j1], :] - verts[facets[:, j2], :]
            ei = np.linalg.norm(tv3, axis=1)
            ineq = ei < mesh_quality_settings["min_edge_len"]
            if np.any(ineq):
                idx = np.nonzero(ineq)
                assert len(idx) == 1
                idx = idx[0]
                assert np.ndim(idx) == 1
                assert np.ndim(ei) == 1
                print idx.size, "zero edges"
                #print idx
                #print tv3 [(ineq), : ] * 1000000.0

                v1 = facets[idx, j1]
                v2 = facets[idx, j2]
                v12_ = np.concatenate((v1[np.newaxis, :], v2[np.newaxis, :]), axis=0)
                v12 = np.concatenate( (v12, v12_), axis=1)

                #fa = idx
                # vertices_to_combine

                nil_edgelen_whichfaces[idx] = True
                any_zero_edge = True

            inan = np.isnan(el[i])
            if np.any(inan):
                idx = np.nonzero(inan)
                assert len(idx) == 1
                idx = idx[0]
                assert np.ndim(idx) == 1
                assert np.ndim(ei) == 1
                print idx.size, "nan edges"

                any_zero_edge = True

                del idx
            del inan

    print v12.shape
    #print v12.transpose()
    v12.sort(axis=0)  # make edges unique
    #sidx = v12.argsort(axis=1)[1, :]  # clump repeated edges
    #sv12 = v12[:, sidx]
    sv12 = v12
    #for di in [1, 0]:
    #    sidx = sv12.argsort(axis=1)[di, :]  # clump repeated edges
    #    sv12 = sv12[:, sidx]
    # Both columns need to be sorted becasue of repeats in each column.
    sidx = np.dot(np.array([1, 10000000]), sv12).argsort()
    sv12 = sv12[:, sidx]
    #print sv12.transpose()
    #print np.diff(sv12, axis=1).transpose()
    #print np.diff(sv12, axis=1)[:, ::2].transpose()
    if CHECK_PAIRED:
        assert np.all((np.diff(sv12, axis=1)[:, ::2]).ravel() == 0)  # Make sure each edge is repeated exactly once.
    map12 = sv12[:, 0::2]  # project map12[1,:] into map12[0,:]

    #print np.diff(sv12, axis=1)[:, :].transpose()
    edge_vects = verts[sv12[0, :]] - verts[sv12[1, :]]


    # Assert that it's always almost-zeros.
    assert np.allclose(edge_vects, 0, mesh_quality_settings["min_edge_len"])
    # Vector subtractionis sligthly more tight than norm but we use the same tolerance here (more conservatirve).
    # Since the actual different is often actually zero, this should hold. Correct the tolerance if this assert failed.


    print( e1[e1 < mesh_quality_settings["min_edge_len"] ])
    print( e1[np.isnan(e1)])

    average_edge = (np.mean(e1)+np.mean(e2)+np.mean(e3))/3.
    # note: the average edge length may be slightely less than this after removing the repeated vertices

    if False:
        # Do the projection using map12
        # project map12[1,:] into map12[0,:]
        #faces[ faces == map12[1,:] ] =
        temp = map12[1, :].copy()
        temp.sort()
        print map12.shape[1], "vertices to be removed. ", temp
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
    degenerate_faces = ineq
    assert np.ndim(ineq) == 1
    idx = np.nonzero(ineq)[0]
    assert np.ndim(idx) == 1
    if len(idx) == 0:
        print "no degenerate area found"
    any_degenerate_area = len(idx) > 0

    nil_areas_whichfaces[idx] = True
    #nf = facets.shape[0]
    for fi in idx:
        assert degenerate_faces[fi]
        if degenerate_faces[fi]:
            #print("face:", fi, facets[fi,:])  # ('face:', 181, array([131,  71, 132]))
            degen_triangle = verts[facets[fi, :], :]  # numverts x 3
            v1 = degen_triangle[1, :] - degen_triangle[0, :]
            v2 = degen_triangle[2, :] - degen_triangle[0, :]
            #print v1,v2
            #print np.cross(v1,v2), np.linalg.norm(np.cross(v1,v2)) * (1000**2) , "(micron^2)"
            #print degen_triangle
            #exit()


            #todo: check all Mesh weirdnesses.
            #then remove those faces and combine those vertices.
            # then remove the faces ith zero area
            # Then fix the T-junctions.
            # Then combine all.
            #Then: optimize the flower.
            #Then: optimise the matrix of arccos & weighted relaxation.
            #Dont do: Decimation: ...


    #print np.nonzero(nil_areas_whichfaces)[0], np.sum(nil_areas_whichfaces)
    #print np.nonzero(nil_edgelen_whichfaces)[0], np.sum(nil_edgelen_whichfaces)
    #"both area and edge are zero"
    lboth = np.nonzero( np.logical_and(nil_areas_whichfaces, nil_edgelen_whichfaces) )[0]
    #"zero-area only"
    lb = np.nonzero( np.logical_and(nil_areas_whichfaces, np.logical_not(nil_edgelen_whichfaces)) )[0]
    #la is not informative
    la = np.nonzero( np.logical_and(np.logical_not(nil_areas_whichfaces), nil_edgelen_whichfaces) )[0]
    #DOES NOT HOLD:
    # assert len(la) == 0  # la, "zero-edge only", should be empty. Because 'zero edge' => 'zero area'
    print "zero-area only", lb
    print "both area and edge are zero", lboth

    #nil_edgelen_whichfaces is-subset-of nil_areas_whichfaces

    #vertices_to_combine = lboth
    #faces_to_remove = lb


    #print "*********"*100
    if fix_them:
        #todo: new_verts
        verts, facets = remove_vertices_and_faces(verts, facets, nil_areas_whichfaces, map12)
        check_faces(facets)
    #***not tested
    #todo: write unit-test

    any_correction = any_zero_edge or any_degenerate_area

    print any_zero_edge , any_degenerate_area
    #assert fix_them != if_assert, (fix_them, if_assert)
    if if_assert:
        if any_correction:
            raise AssertionError
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
    if AREA_DEGENERACY_THRESHOLD is not None:
        degenerates_count = len(facet_areas[facet_areas < AREA_DEGENERACY_THRESHOLD])
        facet_areas[facet_areas < AREA_DEGENERACY_THRESHOLD] = np.nan  # -1
        if degenerates_count > 0:
            print("degenerate triangles", degenerates_count)

    if not return_normals:
        return facet_areas
    else:
        print facet_areas.shape
        assert facet_areas[:, np.newaxis].shape == (nfaces, 1)
        facet_normals = a / np.tile(facet_areas[:, np.newaxis], (1, 3)) / 2.0
        return facet_areas, facet_normals


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
    #from ipdb import set_trace
    from mesh_utils import make_edge_lookup_old
    #pudb.set_trace()
    #set_trace()
    check_faces(facets)
    check_faces3(facets)
    (edges_of_faces, faces_of_edges, vertpairs_of_edges) = \
        make_edge_lookup_old(facets)
        # ****

    # need: faces_of_faces
    # e__nf_x_3 = edges_of_faces[facets]
    # print e__nf_x_3.shape
    # assert e__nf_x_3.shape == (nfaces, 3, 3)
    print edges_of_faces.shape
    nfaces = facets.shape[0]
    assert edges_of_faces.shape == (nfaces, 3)
    f1 = faces_of_edges[edges_of_faces, 0]  # first face of all edges of all faces : nf x 3 -> nf
    f2 = faces_of_edges[edges_of_faces, 1]  # second face of all edges of all faces: nf x 3 -> nf   #3 for 3 sides (edges)
    # one of the two (0:2) is equal to index=face. [index,3, 0:2 ]
    f12 = faces_of_edges[edges_of_faces, :]
    print f12.shape
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
def process2_vertex_resampling_relaxation(verts, facets, iobj):
    global q
    q+=1
    print "q"*1000, q
    #facets = fix_faces_3div2(facets)
    centroids = compute_centroids(verts, facets)
    centroid_normals_normalized = compute_centroid_gradients(centroids, iobj, normalise=True)
    from mesh_utils import make_neighbour_faces_of_vertex
    neighbour_faces_of_vertex = make_neighbour_faces_of_vertex(facets)
    faces_of_faces = build_faces_of_faces(facets)
    assert not np.any(np.isnan(verts.ravel()))  # fine
    new_verts = vertex_resampling(verts, neighbour_faces_of_vertex, faces_of_faces, centroids, centroid_normals_normalized, c=2.0)
    assert not np.any(np.isnan(new_verts.ravel()))  # fails
    return new_verts, facets, centroids  # why does it return facets?



global failure_pairs
failure_pairs = []

#from mesh1.py
def vertex_resampling(verts, neighbour_faces_of_vertex, faces_of_faces, centroids, centroid_normals, c=2.0):
    """ neighbour_faces_of_vertex: *** """

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
        pipj = np.linalg.norm(pi - pj)
        #print "pipj ", pipj, "  arccos= ",np.arccos(mimj)/np.pi*180 #why is it zero??
        assert pipj == np.linalg.norm(pi - pj, ord=2)

        if pipj==0:
            failure_pairs.append( (i, j) )
            #raise TroubledMesh("repeated triangle")
            pipj = 1
        #CAN BE absolute ZERO! ****************
        #FIXME
        assert pipj > 0  # fails # ****

        kij = np.arccos(mimj) / pipj  # radians?
        assert not np.isnan(kij)  # fails
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

        wi = 1.0 + c*ki
        # i_facet is facet (Centroid) index. j_facet is its neighbour facet (centroid). There are three j_facet for an i_facet.
        return wi
    #
    c_ = c  # 2.0  # constant
    vertex_index = 1  # vertex
    #assert vertex_index >= 0
    umbrella_facets = neighbour_faces_of_vertex[vertex_index]  # A list of facets: The indices of faces that vertex vertex_index belongs to.
    print("umbrella_facets: ", umbrella_facets)
    #wa = np.zeros()
    w_list = []
    for i_facet in umbrella_facets:
        # neighbour facet i_facet of Vertex vertex_index
        #three_facets = filter(lambda idx: idx != i_facet, umbrella_facets)
        three_facets = faces_of_faces[i_facet, :]
        w = wi(i_facet, three_facets, c_)  # three_facets should be neighbours of the facet i_facet
        # The weight (based on curvature) of neighbour P_i (facet i.e. centroid),
        w_list.append(w)
        #todo: sparse matrix: w[vi=vertex_index, f2=i_facet] = w
        #todo: store in ...
    #print "w_list ",w_list
    #
    #w seems tobe calculated fine. next: store w_i and cache them for adaptive resampling, for which we need to normalise it across the neighbours.
    nfaces = centroids.shape[0]
    wi_total_array = np.zeros((nfaces,))
    for i_facet in range(nfaces):
        three_facets = faces_of_faces[i_facet, :]
        w = wi(i_facet, three_facets, c_)
        assert not np.isnan(w)  # fails

        wi_total_array[i_facet] = w
    #print wi_total_array
    # The weights are prepared. Now let's resample vertices
    assert not np.any(np.isnan(wi_total_array.ravel()))  # fails

    vertex_index = 1
    #todo: umbrella_Facets = sparse matrix
    #umbrella_facets = np.array(neighbour_faces_of_vertex)  #empty

    umbrella_facets = np.array(neighbour_faces_of_vertex[vertex_index])  # empty
    print "umbrella_facets", umbrella_facets.shape, "****"
    assert np.allclose( wi_total_array[umbrella_facets] - np.array(w_list), 0)

    def lift_verts(verts, centroids):
        new_verts = verts.copy()
        # assign these to a sparse matrix? and  do:  M = M/normalise(M); verts = M * verts
        for vertex_index in range(verts.shape[0]):
            umbrella_facets = np.array(neighbour_faces_of_vertex[vertex_index])
            w = wi_total_array[umbrella_facets]
            assert not np.any(np.isnan(umbrella_facets.ravel()))  # pass
            assert not np.any(np.isnan(wi_total_array.ravel()))  # fails

            assert not np.any(np.isnan(w.ravel()))  # fails
            w = w / np.sum(w)
            assert not np.any(np.isnan(w.ravel()))  # fails

            new_verts[vertex_index, :] = \
                np.dot(w, centroids[umbrella_facets, 0:3])  # (n) * (n x 3)
        assert not np.any(np.isnan(new_verts.ravel()))  # fails
        return new_verts

    assert not np.any(np.isnan(verts.ravel()))  # fine
    r = lift_verts(verts, centroids)
    print r.shape
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


def visualise_gradients(mlab, pos, iobj, arrow_size):
    lm = arrow_size  # 1.  # STEPSIZE
    pos4 = np.concatenate((pos, np.ones((pos.shape[0],1))),axis=1)
    pnormals = - iobj.implicitGradient(pos4)
    pnormals = normalize_vector4_vectorized(pnormals)
    check_vector4_vectorized(pos4)
    xyz = pos
    uvw = pnormals [:,0:3] / 2.
    xx, yy, zz = xyz[:, 0], xyz[:, 1], xyz[:, 2]
    uu, vv, ww = uvw[:, 0], uvw[:, 1], uvw[:, 2]
    #ax.quiver
    #ax.quiver(xx, yy, zz,   uu, vv, ww,  length=np.abs(lm), arrow_length_ratio=0.3, alpha=0.3, pivot="tail")
    #arrow_length_ratio=   length=np.abs(lm)
    #pivot: tail | middle | tip
    #mlab.quiver3d(x_verts,y_verts,z_verts, UVW_normals[:,0],UVW_normals[:,1],UVW_normals[:,2],color=(0,0,0))
    mlab.quiver3d(xx, yy, zz, uu, vv, ww, color=(0, 0, 0), scale_factor=np.abs(lm), line_width=0.5)


def display_simple_using_mayavi_2(vf_list, pointcloud_list, minmax=(-1,1), mayavi_wireframe=False, opacity=1.0, 
        separate=True, gradients_at=None, gradients_from_iobj=None, pointsizes=None, pointcloud_opacity=1.,
        add_noise=[]):
    """Two separate panels"""

    print"Mayavi.", ; sys.stdout.flush()

    from mayavi import mlab
    print "Imported."; sys.stdout.flush()

    if pointsizes is None:
        pointsizes = [0.2]*10

    if type(opacity) is list:
        opacities = opacity  # 1.0
    else:
        opacities = [opacity] + [0.2]*(len(vf_list)-1)  # 1.0, 0.2 #0.1


    for fi in range(len(vf_list)):
        if separate:
            mlab.figure()

        vf = vf_list[fi]
        verts, faces = vf

        if verts is None:
            continue
        if verts.size == 0:
            print("Warning: empty vertices")
            continue
        if faces.size == 0:
            print("Warning: no faces")
            continue

        assert verts.ndim == 2
        assert faces.ndim == 2
        assert verts.shape == (verts.shape[0], 3), str(verts.shape)
        assert faces.shape == (faces.shape[0], 3), str(faces.shape)
        if type(mayavi_wireframe) is list:
            wire_frame1 = mayavi_wireframe[fi]
            assert len(mayavi_wireframe) == len(vf_list)
        else:
            wire_frame1 = mayavi_wireframe
        M = 0.1*0
        if len(add_noise)>0:
            assert len(add_noise) == len(vf_list)
            M = float(add_noise[fi])


        if False:
            _v = verts
            _f = faces
        else:
            _v = verts[faces, :]
            print _v.shape

            _nv = np.prod(faces.shape)
            _v = _v.reshape( (_nv, 3) )
            assert _nv % 3 == 0
            _f = np.arange(_nv).reshape( (_nv/3, 3) )
            #print _nv
            #print _nv % 3, _nv / 3
            #print _f
            #print np.max(_f.ravel())
            #print _v.shape

            #print _f.shape
            #print _v.shape
            #print np.max(_f.ravel()), _v.shape[0]

            #print _f.shape

            qq = _v[_f,:]

            #print np.nonzero( np.isnan(_f.ravel()) )
            #print np.nonzero( np.isinf(_f.ravel()) )
            vv=_f.ravel(); vv.sort()
            #print np.nonzero(vv != np.arange(vv.size))

            #_v = np.concatenate( (_v, np.zeros( (10000, 3) )), axis=0)
            _v = _v + (np.random.rand( _v.shape[0], _v.shape[1] ) -0.5) * M

            #qq = _v[_f, :]
            #assert  np.all( _f.ravel() == np.arange( _f.size ) )
            #print qq
            #print _v
            #assert np.all(qq.ravel() == _v.ravel())

        #v1 = [_v[0]+M*np.random.rand() for vert in verts],
        #v2 = [_v[1]+M*np.random.rand() for vert in verts],
        #v3 = [_v[2]+M*np.random.rand() for vert in verts],
        v1 = [v[0] for v in _v]
        v2 = [v[1] for v in _v]
        v3 = [v[2] for v in _v]
        mlab.triangular_mesh(
                        v1, v2, v3, _f,
                        representation="surface" if not wire_frame1 else "wireframe",
                        opacity=opacities[fi], scale_factor = 100.0)
        #opacity = 0.2 #0.1


        #allpoints are plottedon all panels?
        color_list = [(1, 0, 0), (0, 0, 0), (1, 1, 0), (0, 0, 1), (0,1,0)]
        i = 0
        for c in pointcloud_list:
            #if separate:
            #    if i != fi:
            #        continue
            #print c[:,0:3]
            mlab.points3d(c[:, 0], c[:, 1], c[:, 2], color=color_list[i], scale_factor=pointsizes[i], opacity=pointcloud_opacity)
            i+=1
        del i

        if minmax is not None:
            (RANGE_MIN, RANGE_MAX) = minmax
            x = np.linspace(RANGE_MIN,RANGE_MAX,2).reshape(2,1)
            y = np.zeros((2,1))
            z = np.zeros((2,1))

            mlab.plot3d(x, y, z,line_width=3,name="x-axis")
            mlab.plot3d(y, x, z,line_width=3,name="y-axis")
            mlab.plot3d(z, y, x,line_width=3,name="z-axis")

            mlab.text3d(RANGE_MAX,0,0, "x", scale=0.3)
            mlab.text3d(0,RANGE_MAX,0, "y", scale=0.3)
            mlab.text3d(0,0,RANGE_MAX, "z", scale=0.3)
            mlab.text3d(RANGE_MIN,0,0, "-x", scale=0.3)
            mlab.text3d(0,RANGE_MIN,0, "-y", scale=0.3)
            mlab.text3d(0,0,RANGE_MIN, "-z", scale=0.3)


    def add_random_interior_points(ax, iobj):
        """ Adding random points """
        n=10000
        import basic_types
        # ******
        print avg_edge_len, "WHY USED BEFORE DEFINED?"
        ampl = avg_edge_len
        #ampl = 2
        x = basic_types.make_random_vector_vectorized(n, ampl, 1, type="rand", normalize=False)
        v = iobj.implicitFunction(x)
        x_sel =  x[ v >= 0 , :]
        if x_sel.size ==0:
            print("No points")
            return
        ax.points3d(x_sel[:,0], x_sel[:,1], x_sel[:,2], color=(0,0,0), scale_factor=0.2)

    if gradients_at is not None:
        verts1, faces1 = vf_list[0]
        avg_edge_len = compute_average_edge_length(verts, faces)
        visualise_gradients(mlab, gradients_at, gradients_from_iobj, avg_edge_len / 20.)
    if gradients_from_iobj is not None:
        add_random_interior_points(mlab, gradients_from_iobj)

    mlab.show()
    return


def compute_facets_subdivision_curvatures(verts, facets, iobj):
    """ Calculates a measure of deviation of the Triangle from object gradients.
    returns: curvature for all triangles.
    This function does not create the subdivisions.
    It just computes. The function subdivide_multiple_facets() does the actual subdivision. """
    facet_areas, facet_normals = compute_triangle_areas(verts, facets, return_normals=True)

    nf = facets.shape[0]
    assert facet_areas.shape == (nf,)
    assert facet_normals.shape == (nf, 3)
    assert np.all(np.logical_not(np.isnan(facet_areas[np.logical_not(np.isnan(np.linalg.norm(facet_normals, axis=1)))])))
    #some edges are repeated
    degenerate_faces = np.isnan(facet_areas)

    nn = np.isnan(facet_normals[np.logical_not(degenerate_faces),:])
    print nn.shape
    print np.any(nn, axis=1).shape, "or"
    print facet_areas[np.any(nn, axis=1)]

    assert np.all(np.logical_not(np.isnan(facet_areas.ravel()))), "facet_areas: never nan. But can be zero."

    assert np.all(np.isnan(facet_areas[degenerate_faces]))
    assert np.all(np.logical_not(np.isnan(facet_areas[np.logical_not(degenerate_faces)])))
    assert np.all(np.isnan(facet_normals[degenerate_faces, :]))
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
        mm = mm / np.tile(nn[:,np.newaxis], (1, 3))  # mm: 4 x 3
        mm = mm.transpose()  # 3x4
        e = facet_areas[fi] * np.sum(1. - np.abs(np.dot(n, mm))) / 4.  # sum(,x4)
        #assert np.all(np.dot(n, mm) > -0.0000001 ), "ingrown normal!"
        #e = np.sum(1 - np.abs(np.dot(n, mm)))   # sum(,x4)   #forgot the abs!
        curvatures_array[fi] = e
        #if e<0:
        #    set_trace()

        if fi % 100 == 0:
            print fi, "*   \r", ;import sys; sys.stdout.flush()
    l = curvatures_array[np.logical_not(np.isnan(curvatures_array))].tolist()
    l.sort()
    print "curvature: min,max = ", l[0], l[-1]   # 3.80127650325e-08, 0.0240651184551
    bad_facets_count = np.sum(degenerate_faces)
    #assert bad_facets_count == 0
    return curvatures_array, bad_facets_count


def propagated_subdiv(verts, facets, old_edges):
    """ Reports the indices triangles that are to be subdivided 
    as a propagation of a previous subdivision.
    It returns separatesly the (not-yet-subdivided) triangles 
    with 1,2 and 3 sibdivided sides, in a dictionary."""
    #facets = new_facets
    #verts = new_verts
    #print old_edges
    if len(old_edges) == 0:
        subdivided_edges = np.zeros((0, 2), dtype=np.int)
    else:
        subdivided_edges = np.asarray(old_edges)  # doesnt work if empty
    print subdivided_edges.shape
    assert subdivided_edges.shape[1] == 2
    subdivided_edges.sort(axis=1)
    B = 100000
    BB = np.array([1, B])

    subdiv_edges_codes = np.dot(subdivided_edges, BB)

    fc0 = facets.copy()
    f0 = facets[:, np.newaxis, [0, 1]]
    f1 = facets[:, np.newaxis, [1, 2]]
    f2 = facets[:, np.newaxis, [0, 2]]
    f012 = np.concatenate( (f0, f1, f2), axis=1 )
    f012.sort(axis=2)
    assert np.all(fc0.ravel() == facets.ravel())  # no copy() needed
    all_edges_codes = np.dot(f012, BB)

    #now look for subdiv_edges_codes in all_edges_codes
    #print all_edges_codes
    #print subdiv_edges_codes  #they are so many! 417
    intersec = np.intersect1d(all_edges_codes, subdiv_edges_codes)
    #print intersec
    x_ = np.lib.arraysetops.in1d(all_edges_codes, intersec) # elements of A, A.ravel[x_], that are in B
    print x_.shape, "x_.shape"
    print all_edges_codes.shape, "all_edges_codes.shape"
    print np.prod(all_edges_codes.shape), "prod(all_edges_codes.shape)"
    print intersec.shape, "intersec.shape"
    #assert each_of all_edges_codes.ravel()[x_] in intersec
    edges_which_in1 = all_edges_codes.ravel()[x_] # all edges_which_in1 are in intersec
    #print edges_which_in1
    #exit()
    #print x_.shape, "x_"
    #print all_edges_codes.shape
    #print intersec.shape
    assert np.sum(x_) == intersec.size  # 417
    ticks = x_.reshape(all_edges_codes.shape)
    #print ticks.shape  # -x3
    sides= np.sum(ticks, axis=1)
    #sides.sort()  # there is no point in sorting

    #now I need whose all_edges_codes (i.e. ticks) for which sides==1
    sides_1 = ticks[sides==1, :]
    #print sides_1
    #print edges_which_in1

    #for c in range(4):
    #    print c, ":", np.sum(sides == c)
    """
    0 : 3083
    1 : 232
    2 : 85
    3 : 5
    """
    # "2" will propagate further. But not 3 and 1.
    propag_list = {}
    for c in range(1, 4):
        idx = np.nonzero(sides == c)[0]
        propag_list[c] = idx
        #only propagate triangles with subdivided 1,2,3 sides.
        #exit()

    #edges1 = -1 #which edges are about one side?

    return propag_list, edges_which_in1


def subdivide_multiple_facets(verts_old, facets_old, tobe_subdivided_face_indices):
    """ Use compute_facets_subdivision_curvatures() to calculate tobe_subdivided_face_indices.
    Does not remove vertices => will be valid. But vertices will change: new elements will be appended to it.
    Returns: new vertices and faces.
    Returns: old_edges: The edges that have been removed. This will be used for propagating the subdivision to triangles that their edges are not valid anymore.
    """

    # todo: store subdivided gradients (on top of centroids), to avoid unnecessary calculations. When updating vettices, remove the caches.
    # todo: avoid recomputing

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

        [0.5/(1.+DIP),  0.5/(1.+DIP),  DIP/(1.+DIP)],  # 3
        [DIP/(1.+DIP),  0.5/(1.+DIP),  0.5/(1.+DIP)],  # 4
        [0.5/(1.+DIP),  DIP/(1.+DIP),  0.5/(1.+DIP)]   # 5
        ])  # .transpose()

    global trace_subdivided_facets
    trace_subdivided_facets = []

    verts_old_old = verts_old.copy()

    #raise "not tested yet. re-read/write step by step"
    #allocate space for them

    #new_verts = verts_old + 3*len(tobe_subdivided_face_indices)
    #new_facets = facets_old + 3*len(tobe_subdivided_face_indices)

    provisional_new_verts_count = 3*len(tobe_subdivided_face_indices)
    provisional_new_facets_count = 3*len(tobe_subdivided_face_indices)
    nverts_old = verts_old.shape[0]
    nfaces_old = facets_old.shape[0]
    new_verts = np.zeros((nverts_old+provisional_new_verts_count, 3), dtype=float)
    new_facets = np.zeros((nfaces_old+provisional_new_facets_count, 3), dtype=int)
    #set_trace()
    new_verts[:nverts_old, :] = verts_old
    new_facets[:nfaces_old, :] = facets_old

    #on number of added vertices:
    #problem: there may be repeated (Redundant) vertices. (as well as T-junctions)
    #also later check for faces with repeated edges. (which can be another cause of null normals)

    print "Subdividing:"

    #edges that need to be divided because of the neighbouring triangle subdivided
    #subdivided_edges = []
    old_edges = []

    new_vertex_counter = nverts_old
    new_facet_counter = nfaces_old
    for subdiv_i in range(len(tobe_subdivided_face_indices)):

        fi = tobe_subdivided_face_indices[subdiv_i]
        oldtriangle = verts_old[facets_old[fi, :], :]  # numverts x 3
        assert oldtriangle.shape == (3, 3)
        VVV = oldtriangle  # (nv=3) x 3

        # new verices
        m0123 = np.dot(np.dot(centroidmaker_matrix, subdiv_vert_matrix), VVV)
        assert m0123.shape == (4, 3)
        subdiv_centroids = m0123

        vxyz_0123 = np.dot(subdiv_vert_matrix, VVV)  # not efficient
        assert vxyz_0123.shape == (6, 3)

        #tobeadded_verts = m0123
        #tobeadded_verts = m123
        #subdivision = oldtriangle

        #mini_verts = np.concatenate( (oldtriangle, tobeadded_verts), axis=0)

        # adding new verts and facets

        #*********
        # indices of original points
        #v012 = facets_old[fi, :]  # range(0, 3)  #
        #v345 = np.arange(new_vertex_counter, new_vertex_counter+3, dtype=int)   #range(3, 6)
        v012 = facets_old[fi, :].tolist()  # range(0, 3)  #
        v345 = range(new_vertex_counter, new_vertex_counter+3)

        v345_xyz = vxyz_0123[3:6, :]  # only pick the new ones

        assert len(v345) == 3
        new_verts[(new_vertex_counter):(new_vertex_counter+3), :] = v345_xyz

        new_vertex_counter += 3

        # facet's vertex indices
        v012345 = np.array(v012 + v345, dtype=int)

        mini_faces_l = [[0, 3, 5], [3, 1, 4], [5, 4, 2], [3, 4, 5]]  # 0,3,1,4,2,5

        mini_faces = v012345[np.array(mini_faces_l)]

        original_facet_index = fi
        #print original_facet_index
        #print new_facets.shape
        #print fi
        old_edges .append(tuple(new_facets[original_facet_index, [0, 1]].tolist()))
        old_edges .append(tuple(new_facets[original_facet_index, [1, 2]].tolist()))
        old_edges .append(tuple(new_facets[original_facet_index, [2, 0]].tolist()))

        new_facets[original_facet_index, :] = mini_faces[0, :]
        new_facets[new_facet_counter:(new_facet_counter+3), :] = mini_faces[1:(1+3), :]
        assert mini_faces.shape[0] == (1+3)
        #trace_subdivided_facets += range(new_facet_counter, (new_facet_counter+3))
        trace_subdivided_facets += range(new_facet_counter, (new_facet_counter+3)) + [fi]  # include the face which reuses the old face's index
        # trace_subdivided_facets will contain indices of faces

        new_facet_counter += 3

        #return mini_verts, mini_faces
        #numsubdiv = 4

        if subdiv_i % 10 == 0:
            print subdiv_i , " @                     \r", ;import sys; sys.stdout.flush()
            #," @", new_facet_counter/3, new_vertex_counter/3
    print ""

    print new_verts.shape[0], new_vertex_counter

    assert new_verts.shape[0] == new_vertex_counter
    assert new_facets.shape[0] == new_facet_counter
    assert provisional_new_verts_count+nverts_old == new_vertex_counter, "vector consistency"
    assert provisional_new_facets_count+nfaces_old == new_facet_counter, "face consistency"
    assert len(trace_subdivided_facets) == 0 or np.max(np.array(trace_subdivided_facets)) < new_facet_counter
    #return new_verts, new_facets

    #axes = tuple(range(np.ndim(new_verts)))
    #print axes
    #print new_verts.shape
    #print np.all(verts_old_old == new_verts, axis=None)  # (0,) ) # axes)
    #exit()
    return new_verts, new_facets, old_edges


def subdivide_1to2_multiple_facets(verts2, facets2, edges_with_1_side):
    #todo: copy some code from propagated_subdiv()
    #check which of these edges still exist in faces. (Each should be there only once. In this context.)
    #remove them and add more.
    #refactor the code copied from propagated_subdiv() into function
    print "good"
    exit()
    return verts2, facets2


def do_subdivision(verts, facets, iobj, curvature_epsilon):
    assert not np.any(np.isnan(facets.ravel()))
    assert not np.any(np.isnan(verts.ravel()))  # fails

    curvatures, bad_facets_count = compute_facets_subdivision_curvatures(verts, facets, iobj)

    assert np.sum(np.isnan(curvatures)) == 0, "NaN"
    curvatures[np.isnan(curvatures)] = 0  # treat NaN curvatures as zero curvature => no subdivision

    which_facets = np.arange(facets.shape[0])[ curvatures > curvature_epsilon ]

    verts2, facets2, old_edges = subdivide_multiple_facets(verts, facets, which_facets)
    global trace_subdivided_facets  # third implicit output
    #verts2, facets2 = verts_subdivided, facets_subdivided

    #list_facets_with_1_side = []
    list_edges_with_1_side = []
    while True:
        propag_list, edges_which_in1 = propagated_subdiv(verts2, facets2, old_edges)
        for k in propag_list:
            print "%d:"%(k,), len(propag_list[k]), propag_list[k].shape

        facets_with_2_or_3_sides = np.concatenate( (propag_list[2], propag_list[3]), axis=0 )
        # what if those faces dont exist anymore in the next round?
        # list_facets_with_1_side += [propag_list[1]]
        list_edges_with_1_side += [edges_which_in1]


        print facets_with_2_or_3_sides.shape
        if facets_with_2_or_3_sides.size ==0:
            print facets_with_2_or_3_sides.shape
            print "COOL"
            break
        verts2, facets2, old_edges = subdivide_multiple_facets(verts2, facets2, facets_with_2_or_3_sides)
        #todo: merge with above call

    # Finished with 2 or 3 sides.

    # Now 1 side:
    n1 = 0
    for i in range(len(list_edges_with_1_side)):
        farr = list_edges_with_1_side[i]
        assert farr.size == farr.shape[0]
        assert len(farr.shape) == 1
        n1 += farr.size
    edges_with_1_side = np.zeros((n1,), dtype=np.int)
    n1 = 0
    for i in range(len(list_edges_with_1_side)):
        farr = list_edges_with_1_side[i]
        n2 = n1 + farr.size
        edges_with_1_side[n1:n2] = farr
        n1 = n2

    #print list_edges_with_1_side
    #print edges_with_1_side


    verts2, facets2 = subdivide_1to2_multiple_facets(verts2, facets2, edges_with_1_side)
    exit()

    #???
    #v5, f5, old_edges = subdivide_multiple_facets(verts2, facets2, facets_with_2_or_3_sides)

    print("Subdivision applied.");sys.stdout.flush()
    return verts2, facets2


def demo_everything():
    """ Base on demo_combination_plus_qem """
    curvature_epsilon = 1. / 1000.  # a>eps  1/a > 1/eps = 2000
    VERTEX_RELAXATION_ITERATIONS_COUNT = 1
    SUBDIVISION_ITERATIONS_COUNT = 1  # 2  # 5+4

    global STEPSIZE
    from example_objects import make_example_vectorized
    iobj = make_example_vectorized(
        #"rcube_vec")  #
        #"rdice_vec")  #
        #"cube_example") # problem: zero facet areas
        "ell_example1")  #+    
        # "bowl_15_holes")  # works too. But too many faces => too slow, too much memory. 32K?
    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-3, +5, 0.2)

    iobj = two_bricks()

    from stl_tests import make_mc_values_grid
    gridvals = make_mc_values_grid(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE, old=False)
    verts, facets = vtk_mc(gridvals, (RANGE_MIN, RANGE_MAX, STEPSIZE))
    print("MC calculated.");sys.stdout.flush()



    #display_simple_using_mayavi_2( [(verts, facets), ],
    #   pointcloud_list=[ ], pointcloud_opacity=0.2,
    #   mayavi_wireframe=[True], opacity=[1],
    #   gradients_at=None, separate=False, gradients_from_iobj=None,
    #   minmax=(RANGE_MIN,RANGE_MAX)  )
    #exit()

    #print verts.shape, verts.shape[0]/2.
    print facets.shape, facets.shape[0]/3.
    #exit()
    #check_faces(facets)


    #facets = fix_faces_3div2(facets)

    #print ("sdfsf")
    #exit()

    assert not np.any(np.isnan(verts.ravel()))  # fine

    any_mesh_correction = check_degenerate_faces(verts, facets, "dontfix")
    print "any_mesh_correction1", any_mesh_correction
    for qq in range(5):
        print "qq=", qq
        verts, facets, any_mesh_correction = check_degenerate_faces(verts, facets, "fix")
        print "any_mesh_correction2:", qq, any_mesh_correction
    #if False:
    check_degenerate_faces(verts, facets, "assert")
    # COOL ! NOW WORKS!

    old_verts, old_facets = verts, facets
    assert not np.any(np.isnan(verts.ravel()))  # fine

    print np.diff(verts[ [352, 363], : ], axis=0)  # is zero!
    print np.diff(verts[ [361, 373], : ], axis=0)  # is zero!
    #*********
    #exit()

    for i in range(VERTEX_RELAXATION_ITERATIONS_COUNT):

        #facets = fix_faces_3div2(facets)

        verts, facets_not_used, centroids = process2_vertex_resampling_relaxation(verts, facets, iobj)

        #facets = fix_faces_3div2(facets)

        assert not np.any(np.isnan(verts.ravel()))  # fails
        print("Vertex relaxation applied.");sys.stdout.flush()
        verts, facets_not_used, any_mesh_correction = check_degenerate_faces(verts, facets_not_used, "assert")
        assert not np.any(np.isnan(verts.ravel()))  # fails
        if any_mesh_correction:
            print("mesh correction needed")
            exit()

    print "failure_pairs"  # list of pairs that have zero distance in weighted resampling.
    #print failure_pairs
    fpna = np.asarray(failure_pairs).ravel()
    #print facets[fpna, :]
    #print "unique triangles:", np.unique(fpna)
    #print "unique vertices:", np.unique(facets[fpna, :].ravel())
    #print "***fpna****"
    #print facets[fpna,:]
    coords = verts[facets[fpna, :]]
    #print coords.shape  # 16x3x3
    #print coords.reshape(16, 9)

    #quick_vis(verts, facets, fpna)

    #quick_vis(old_verts, old_facets, fpna)

    #compute_facets_subdivision_curvatures
    #subdivide_multiple_facets
    #process2_vertex_resampling_relaxation
    #apply_new_projection

    total_subdivided_facets = []
    for i in range(SUBDIVISION_ITERATIONS_COUNT):

        verts, facets = do_subdivision(verts, facets, iobj, curvature_epsilon)
        global trace_subdivided_facets  # third implicit output
        verts4_subdivided = verts  # ??
        facets3_subdivided = facets

        verts4_subdivided, facets3_subdivided, any_mesh_correction = check_degenerate_faces(verts4_subdivided, facets3_subdivided, "assert")

        total_subdivided_facets += trace_subdivided_facets  # old face indices remain valid

        for i in range(VERTEX_RELAXATION_ITERATIONS_COUNT):
            print "i", "="*10, i
            verts, facets_not_used, centroids = process2_vertex_resampling_relaxation(verts, facets, iobj)
            print("Vertex relaxation applied.");sys.stdout.flush()
            verts, facets_not_used, any_mesh_correction = check_degenerate_faces(verts, facets_not_used, "assert")

    #centroids, new_centroids = apply_new_projection(verts, facets, iobj)
    from ohtake_surface_projection import set_centers_on_surface__ohtake

    average_edge = compute_average_edge_length(verts, facets)

    #up to here
    verts, facets, any_mesh_correction = check_degenerate_faces(verts, facets, "assert")

    c3 = np.mean(verts[facets[:], :], axis=1)
    old_centroids = np.concatenate((c3, np.ones((c3.shape[0], 1))), axis=1)

    nones_map = old_centroids[:, 0]*0 > 100  # all False
    new_centroids = old_centroids.copy()
    set_centers_on_surface__ohtake(iobj, new_centroids, average_edge, nones_map)
    #new_centroids is the output


    verts, facets, any_mesh_correction = check_degenerate_faces(verts, facets, "assert")

    #neighbour_faces_of_vertex
    vertex_neighbours_list = mesh_utils.make_neighbour_faces_of_vertex(facets)
    centroid_gradients = compute_centroid_gradients(new_centroids, iobj)
    #nv1  =
    new_verts_qem = \
        vertices_apply_qem3(verts, facets, new_centroids, vertex_neighbours_list, centroid_gradients)
    #verts = nv1
    #new_verts_qem = verts

    new_verts_qem, facets, any_mesh_correction = check_degenerate_faces(new_verts_qem, facets, "assert")

    #

    alpha = 0.
    new_verts_qem_alpha = (new_verts_qem * alpha + verts * (1-alpha))

    chosen_facet_indices = np.array(total_subdivided_facets, dtype=int)

    centroids2, new_centroids2 = old_centroids[chosen_facet_indices], new_centroids[chosen_facet_indices]

    # move the following code into subdivide_multiple_facets() (?)
    if chosen_facet_indices.size == 0:
        chosen_subset_of_facets = np.zeros((0,), dtype=int)
    else:
        chosen_subset_of_facets = facets[chosen_facet_indices, :]

    highlighted_vertices = np.arange(100, 200)
    hv = new_verts_qem[highlighted_vertices, :]

    check_degenerate_faces(new_verts_qem_alpha, facets, "assert")
    check_degenerate_faces(new_verts_qem, facets, "assert")

    display_simple_using_mayavi_2( [(new_verts_qem_alpha, facets),(new_verts_qem, facets), ],
       pointcloud_list=[ hv ], pointcloud_opacity=0.2,
       mayavi_wireframe=[False,False], opacity=[0.4*0, 1, 0.9], gradients_at=None, separate=False, gradients_from_iobj=None,
       minmax=(RANGE_MIN,RANGE_MAX)  )


if __name__ == '__main__':

    demo_choise = 8
    if demo_choise == 7:
        demo_combination_plus_qem()  # subdivision + projection + qem
    elif demo_choise == 8:
        demo_everything()  #
    else:
        print "Error"


#from ipdb import set_trace

#from stl_tests import display_simple_using_mayavi_vf1
from ohtake_surface_projection import display_simple_using_mayavi_


#def adaptive_subdivision(self):
#    pass







def compute_average_edge_length(verts, faces):
    nfaces = faces.shape[0]
    expand = verts[faces, :]
    assert expand.shape == (nfaces, 3, 3)
    assert expand[:, 2, :].shape == (nfaces, 3)
    ea_sum = 0.
    for i in range(3):
        i1 = i
        i2 = (i+1) % 3
        e1 = np.linalg.norm(expand[:, i1, :] - expand[:, i2, :])
        ea_sum += np.mean(e1)
    return ea_sum / 3.



def degenerate_facets():
    facet_areas, facet_normals = compute_triangle_areas(verts, facets, return_normals=True)

    nf = facets.shape[0]
    degenerate_faces = np.isnan(facet_areas)
    indices = np.arange(nf)[degenerate_faces]

    assert facet_areas.shape == (nf,)
    assert facet_normals.shape == (nf, 3)
    assert np.all(np.logical_not(np.isnan(facet_areas[np.logical_not(np.isnan(np.linalg.norm(facet_normals, axis=1)))])))
    assert np.all(np.isnan(facet_areas[degenerate_faces]))
    assert np.all(np.logical_not(np.isnan(facet_areas[np.logical_not(degenerate_faces)])))
    assert np.all(np.isnan(facet_normals[degenerate_faces, :]))
    assert np.all(np.logical_not(np.isnan(facet_normals[np.logical_not(degenerate_faces),:])))

    return indices


def fix_degenerate_Faces():
    pass


#def check_degenerate_faces1(verts, facets, degenerate_faces):
#    nf = facets.shape[0]
#    for fi in range(nf):
#        if degenerate_faces[fi]:
#            degen_triangle = verts[facets[fi, :], :]  # numverts x 3
#            print degen_triangle
#
#            #assert not degenerate_faces[fi]
#
#            triangle = verts[facets[fi, :], :]  # numverts x 3
#            assert triangle.shape == (3, 3)






def apply_new_projection(verts, facets, iobj):
    from ohtake_surface_projection import set_centers_on_surface__ohtake

    average_edge = compute_average_edge_length(verts, facets)

    c3 = np.mean(verts[facets[:], :], axis=1)
    # add extra points
    #c3 = np.concatenate((c3, c3+STEPSIZE*0.1, c3+STEPSIZE*(-0.1)), axis=0)
    #c3 = np.concatenate((c3,), axis=0)
    centroids = np.concatenate((c3, np.ones((c3.shape[0], 1))), axis=1)

    nones_map = centroids[:, 0]*0 > 100  # all False
    new_centroids = centroids.copy()
    set_centers_on_surface__ohtake(iobj, new_centroids, average_edge, nones_map)
    #new_centroids is the output

    return centroids, new_centroids

    #display_simple_using_mayavi_2( [ (verts, facets), ], pointcloud_list=[ centroids, new_centroids],
    #   mayavi_wireframe=[False], opacity=[0.2,], gradients_at=None, separate=False, gradients_from_iobj=None,
    #   pointsizes=[0.02, 0.05]) # minmax=(RANGE_MIN,RANGE_MAX))


def demo_combination_actually_do():
    """ Now combination of vertex relaxation + subdivision only. Both are iterative. """

    #1. / 2000/50. # 43K
    #curvature_epsilon = 1. / 2000/2.  # most points
    #curvature_epsilon = 1. / 1000.
    curvature_epsilon = 1. / 2000.
    VERTEX_RELAXATION_ITERATIONS_COUNT = 0
    SUBDIVISION_ITERATIONS_COUNT = 2  # 5+4

    from example_objects import make_example_vectorized
    iobj = make_example_vectorized("ell_example1")  # "bowl_15_holes" works too
    (RANGE_MIN,RANGE_MAX, STEPSIZE) = (-3, +5, 0.2)


    from stl_tests import make_mc_values_grid
    gridvals = make_mc_values_grid(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE, old=False)
    verts, facets = vtk_mc(gridvals, (RANGE_MIN, RANGE_MAX, STEPSIZE) )
    print("MC calculated.");sys.stdout.flush()

    old_verts, old_facets = verts, facets



    #apply_new_projection(old_verts, old_facets, iobj); exit()

    #from mesh_utils import mesh_invariant
    #mesh_invariant(facets)


    #new_verts3, new_facets3, trace_subdivided_facets = \
    #    process4_combine_both(verts, facets, iobj, curvature_epsilon, 0)

    for i in range(VERTEX_RELAXATION_ITERATIONS_COUNT):
        verts, facets_not_used, centroids = process2_vertex_resampling_relaxation(verts, facets, iobj)
        print("Vertex relaxation applied.");sys.stdout.flush()

    total_subdivided_facets = []
    for i in range(SUBDIVISION_ITERATIONS_COUNT):

        #not tested:
        verts, facets = do_subdivision(verts, facets, iobj, curvature_epsilon)
        global trace_subdivided_facets  # third implicit output

        total_subdivided_facets += trace_subdivided_facets  # old face indices remain valid

        for i in range(VERTEX_RELAXATION_ITERATIONS_COUNT):
            verts, facets_not_used, centroids = process2_vertex_resampling_relaxation(verts, facets, iobj)
            print("Vertex relaxation applied.");sys.stdout.flush()

    #centroids, new_centroids = apply_new_projection(verts, facets, iobj)
    #display_simple_using_mayavi_2( [ (verts, facets), ], pointcloud_list=[ centroids, new_centroids],
    #   mayavi_wireframe=[False], opacity=[0.2,], gradients_at=None, separate=False, gradients_from_iobj=None,
    #   pointsizes=[0.02, 0.05]) # minmax=(RANGE_MIN,RANGE_MAX))
    #exit()

    chosen_facet_indices = np.array(total_subdivided_facets)

    #centroids2, new_centroids2 = centroids[chosen_facet_indices], new_centroids[chosen_facet_indices]

    if chosen_facet_indices.size == 0:
        chosen_subset_of_facets = np.zeros((0,), dtype=int)
    else:
        chosen_subset_of_facets = facets[chosen_facet_indices, :]
    display_simple_using_mayavi_2( [ (old_verts, old_facets), (verts, chosen_subset_of_facets), ],
       mayavi_wireframe=[False, True], opacity=[0.2, 1, 0.1], gradients_at=None, separate=False, gradients_from_iobj=None,
       minmax=(RANGE_MIN,RANGE_MAX),
       pointcloud_list=[], pointsizes=[]
       )


def demo_combination_actually_do_plus_centroid_projection():
    """ Now combination of vertex relaxation + subdivision only. Both are iterative. """

    # 1. / 2000/50. # 43K
    # curvature_epsilon = 1. / 2000/2.  # most points
    curvature_epsilon = 1. / 1000.
    #curvature_epsilon = 1. / 2000.
    VERTEX_RELAXATION_ITERATIONS_COUNT = 0
    SUBDIVISION_ITERATIONS_COUNT = 2  # 5+4

    from example_objects import make_example_vectorized
    iobj = make_example_vectorized( "ell_example1")  #
        # "bowl_15_holes")  # works too. But too many faces => too slow, too much memory. 32K?
    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-3, +5, 0.2)

    from stl_tests import make_mc_values_grid
    gridvals = make_mc_values_grid(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE, old=False)
    verts, facets = vtk_mc(gridvals, (RANGE_MIN, RANGE_MAX, STEPSIZE))
    print("MC calculated.");sys.stdout.flush()

    old_verts, old_facets = verts, facets



    #apply_new_projection(old_verts, old_facets, iobj); exit()

    #from mesh_utils import mesh_invariant
    #mesh_invariant(facets)


    #new_verts3, new_facets3, trace_subdivided_facets = \
    #    process4_combine_both(verts, facets, iobj, curvature_epsilon, 0)

    for i in range(VERTEX_RELAXATION_ITERATIONS_COUNT):
        verts, facets_not_used, centroids = process2_vertex_resampling_relaxation(verts, facets, iobj)
        print("Vertex relaxation applied.");sys.stdout.flush()

    total_subdivided_facets = []
    for i in range(SUBDIVISION_ITERATIONS_COUNT):

        #e_array, bad_facets_count = compute_facets_subdivision_curvatures(verts, facets, iobj)
        #assert np.sum(np.isnan(e_array)) == 0, "NaN"
        #e_array[np.isnan(e_array)] = 0  # treat NaN curvatures as zero curvature => no subdivision
        #which_facets = np.arange(facets.shape[0])[ e_array > curvature_epsilon ]
        #verts4_subdivided, facets3_subdivided, oe = subdivide_multiple_facets(verts, facets, which_facets)
        #global trace_subdivided_facets  # third implicit output
        #verts, facets = verts4_subdivided, facets3_subdivided
        #print("Subdivision applied.");sys.stdout.flush()
        verts, facets = do_subdivision(verts, facets, iobj, curvature_epsilon)
        global trace_subdivided_facets  # third implicit output
        #verts4_subdivided = verts
        #facets3_subdivided = facets

        total_subdivided_facets += trace_subdivided_facets  # old face indices remain valid

        for i in range(VERTEX_RELAXATION_ITERATIONS_COUNT):
            verts, facets_not_used, centroids = process2_vertex_resampling_relaxation(verts, facets, iobj)
            print("Vertex relaxation applied.");sys.stdout.flush()

    centroids, new_centroids = apply_new_projection(verts, facets, iobj)
    #display_simple_using_mayavi_2( [ (verts, facets), ], pointcloud_list=[ centroids, new_centroids],
    #   mayavi_wireframe=[False], opacity=[0.2,], gradients_at=None, separate=False, gradients_from_iobj=None,
    #   pointsizes=[0.02, 0.05]) # minmax=(RANGE_MIN,RANGE_MAX))
    #exit()

    chosen_facet_indices = np.array(total_subdivided_facets)

    centroids2, new_centroids2 = centroids[chosen_facet_indices], new_centroids[chosen_facet_indices]

    if chosen_facet_indices.size == 0:
        chosen_subset_of_facets = np.zeros((0,), dtype=int)
    else:
        chosen_subset_of_facets = facets[chosen_facet_indices, :]
    display_simple_using_mayavi_2( [ (old_verts, old_facets), (verts, chosen_subset_of_facets), ],
       mayavi_wireframe=[False, True], opacity=[0.2, 1, 0.1], gradients_at=None, separate=False, gradients_from_iobj=None,
       minmax=(RANGE_MIN,RANGE_MAX),
       pointcloud_list=[ centroids2, new_centroids2], pointsizes=[0.01, 0.02]    # centroids
       )


def get_A_b(vertex_id, nlist_numpy, centroids, centroid_gradients):
    #nlist = self.vertex_neighbours_list[vertex_id]
    #nai = np.array(nlist)
    nai = nlist_numpy
    center_array = centroids[nai, :]

    #note some centers may not be projected successfully in the previous step
    not_projected_successfully = np.isnan(center_array[:].ravel())
    if np.any(not_projected_successfully):
        pass

    normals = centroid_gradients[nai, :]  #why do we have repeats??
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

    x0 = np.zeros((3, 1))

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
        b += -np.dot(nnt, p_i - x0)

        # IN PROGRESS

    assert not np.any(np.isnan(A.ravel()))
    assert not np.any(np.isnan(b.ravel()))

    return A, b


def vertices_apply_qem3(verts, facets, centroids, vertex_neighbours_list, centroid_gradients):
    #based on quadratic_optimise_vertices(self, alpha=1.0)
    assert not centroids is None
    assert not vertex_neighbours_list is None
    assert not centroid_gradients is None

    #alpha = 1.0
    nvert = verts.shape[0]
    assert nvert == len(vertex_neighbours_list)

    result_verts_ranks = np.zeros((nvert,), dtype=int)
    assert verts.shape == (nvert, 3)
    new_verts = np.zeros((nvert, 3))

    for vertex_id in range(nvert):

        vi = vertex_id
        nlist = vertex_neighbours_list[vertex_id]
        nai = np.array(nlist)
        A, b = get_A_b(vi, nai, centroids, centroid_gradients)
        #print A, b

        ###
        #A, b = self.get_A_b(vi)

        u, s, v = np.linalg.svd(A)
        assert np.allclose(A, np.dot(u, np.dot(np.diag(s), v)))
        #print(s)  # [  1.48148148e+01   1.67928330e-15   1.01592270e-50]
        assert s[0] == np.max(s)
        #print( s / s[0] )  # [  1.00000000e+00   1.13351623e-16   6.85747820e-52]

        tau = 10. ** 3.
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

        x = verts[vi, 0:3, np.newaxis]
        assert x.shape == (3, 1)

        y = np.dot(v, x).copy()
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
        new_x = np.dot(np.transpose(v), y)
        #print(x.ravel(), " -> ", new_x.ravel())
        #print("    delta=", (new_x - x).ravel())

        new_verts[vi, 0:3] = new_x[:, 0]
        #self.new_verts[vi,3] = 1

        assert x.shape == (3, 1)
        # Apply alpha
        #new_verts[vi, 0:3] = new_x[:, 0] * alpha + x[:, 0] * (1.0-alpha)
        new_verts[vi, 0:3] = new_x[:, 0]

        if not np.all(np.abs(utb.ravel()[rank:] ) < 0.0001):
            #print("s", s.ravel()/s[0], "   utb", utb.ravel()/s[0])
            pass
        result_verts_ranks[vi] = rank

        #exit()
    print("max rank = ", np.max(result_verts_ranks))
    print("min rank = ", np.min(result_verts_ranks))
    if not np.min(result_verts_ranks) >= 1:
        print("Warning: assertion: np.min(result_verts_ranks) >= 1 failed." )

    if False:
        assert np.min(result_verts_ranks) >= 1
    return new_verts




import mesh_utils

def demo_combination_plus_qem():
    """ Now with QEM """
    curvature_epsilon = 1. / 1000.  # a>eps  1/a > 1/eps = 2000
    VERTEX_RELAXATION_ITERATIONS_COUNT = 0
    SUBDIVISION_ITERATIONS_COUNT = 0  # 2  # 5+4

    from example_objects import make_example_vectorized
    iobj = make_example_vectorized(
        #"rcube_vec")  #
        #"rdice_vec")  #
        #"cube_example");
        "ell_example1")  #
        # "bowl_15_holes")  # works too. But too many faces => too slow, too much memory. 32K?
    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-3, +5, 0.2)

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


    from stl_tests import make_mc_values_grid
    gridvals = make_mc_values_grid(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE, old=False)
    verts, facets = vtk_mc(gridvals, (RANGE_MIN, RANGE_MAX, STEPSIZE))
    print("MC calculated.");sys.stdout.flush()

    old_verts, old_facets = verts, facets

    for i in range(VERTEX_RELAXATION_ITERATIONS_COUNT):
        verts, facets_not_used, centroids = process2_vertex_resampling_relaxation(verts, facets, iobj)
        print("Vertex relaxation applied.");sys.stdout.flush()



    #compute_facets_subdivision_curvatures
    #subdivide_multiple_facets
    #process2_vertex_resampling_relaxation
    #apply_new_projection


    total_subdivided_facets = []
    for i in range(SUBDIVISION_ITERATIONS_COUNT):

        #e_array, bad_facets_count = compute_facets_subdivision_curvatures(verts, facets, iobj)
        #assert np.sum(np.isnan(e_array)) == 0, "NaN"
        #e_array[np.isnan(e_array)] = 0  # treat NaN curvatures as zero curvature => no subdivision
        #which_facets = np.arange(facets.shape[0])[ e_array > curvature_epsilon ]
        #verts4_subdivided, facets3_subdivided, oe = subdivide_multiple_facets(verts, facets, which_facets)
        #global trace_subdivided_facets  # third implicit output
        #verts, facets = verts4_subdivided, facets3_subdivided
        #print("Subdivision applied.");sys.stdout.flush()
        verts, facets = do_subdivision(verts, facets, iobj, curvature_epsilon)
        global trace_subdivided_facets  # third implicit output
        #verts4_subdivided = verts
        #facets3_subdivided = facets

        total_subdivided_facets += trace_subdivided_facets  # old face indices remain valid

        for i in range(VERTEX_RELAXATION_ITERATIONS_COUNT):
            verts, facets_not_used, centroids = process2_vertex_resampling_relaxation(verts, facets, iobj)
            print("Vertex relaxation applied.");sys.stdout.flush()

    #centroids, new_centroids = apply_new_projection(verts, facets, iobj)
    from ohtake_surface_projection import set_centers_on_surface__ohtake

    average_edge = compute_average_edge_length(verts, facets)

    c3 = np.mean(verts[facets[:], :], axis=1)
    old_centroids = np.concatenate((c3, np.ones((c3.shape[0], 1))), axis=1)

    nones_map = old_centroids[:, 0]*0 > 100  # all False
    new_centroids = old_centroids.copy()
    set_centers_on_surface__ohtake(iobj, new_centroids, average_edge, nones_map)
    #new_centroids is the output


    # The two CHOICEs are equaivalent. Two rewrites of the same method.
    CHOICE = 1
    if CHOICE == 1:
        #neighbour_faces_of_vertex
        vertex_neighbours_list = mesh_utils.make_neighbour_faces_of_vertex(facets)
        centroid_gradients = compute_centroid_gradients(new_centroids, iobj)
        #nv1  =
        new_verts_qem = \
            vertices_apply_qem3(verts, facets, new_centroids, vertex_neighbours_list, centroid_gradients)
        #verts = nv1
        #new_verts_qem = verts

    elif CHOICE == 2:
        import mesh1
        m = mesh1.Mesh_1(facets, verts)
        m.build_centroids()
        m.build_neighbours()

        #m.faces = faces
        #m.verts = verts
        m.centroids = new_centroids
        #m.vertex_neighbours_list = None
        #m.centroid_gradients = None
        #m.facet_areas = None

        m.evaluate_centroid_gradients(iobj)
        do_qem = True
        if do_qem:
            #both not necessary here!
            #m.update_centroids_and_gradients(iobj)
            #m.update_centroids_and_gradients(iobj)
            
            # if not qem_breakdown:
            m.quadratic_optimise_vertices(1)
            m.verts = m.new_verts
            new_verts_qem = m.verts
    #

    alpha = 0.
    new_verts_qem_alpha = (new_verts_qem * alpha + verts * (1-alpha))

    chosen_facet_indices = np.array(total_subdivided_facets, dtype=int)

    centroids2, new_centroids2 = old_centroids[chosen_facet_indices], new_centroids[chosen_facet_indices]

    # move the following code into subdivide_multiple_facets() (?)
    if chosen_facet_indices.size == 0:
        chosen_subset_of_facets = np.zeros((0,), dtype=int)
    else:
        chosen_subset_of_facets = facets[chosen_facet_indices, :]


    highlighted_vertices = np.arange(100, 200)
    hv = new_verts_qem[highlighted_vertices, :]

    display_simple_using_mayavi_2( [(new_verts_qem_alpha, facets),(new_verts_qem, facets), ],
       pointcloud_list=[ hv ],
       mayavi_wireframe=[False,False], opacity=[0.4*0, 1, 0.9], gradients_at=None, separate=False, gradients_from_iobj=None,
       minmax=(RANGE_MIN,RANGE_MAX)  )

    """
    display_simple_using_mayavi_2( [(new_verts_qem_alpha, facets),],
       pointcloud_list=[],   
       mayavi_wireframe=[False], opacity=[0.2, 1, 0.9], gradients_at=None, separate=False, gradients_from_iobj=None,
       minmax=(RANGE_MIN,RANGE_MAX)  )
    """

    #display_simple_using_mayavi_2( [    #(old_verts, old_facets), (verts, chosen_subset_of_facets),
    #    (new_verts_qem_alpha, facets)],
    #   mayavi_wireframe=[False, True], opacity=[0.2, 1, 0.9], gradients_at=None, separate=False, gradients_from_iobj=None,
    #   minmax=(RANGE_MIN,RANGE_MAX)  ) #,
    #   #pointcloud_list=[ centroids2, new_centroids2], pointsizes=[0.01, 0.02]    # centroids
    #   #)





"""
from pycallgraph import PyCallGraph
from pycallgraph.output import GraphvizOutput

with PyCallGraph(output=GraphvizOutput()):
    a = 4
    print (a)
"""


def make_test_vf_1():
    faces = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
    verts = np.random.randn((4, 3))
    return verts, faces

def make_test_vf_2():
    faces = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
    verts = np.random.randn((8, 3))
    return verts, faces
