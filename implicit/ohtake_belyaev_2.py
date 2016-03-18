from ipdb import set_trace

from vtk_mc import vtk_mc
#from stl_tests import display_simple_using_mayavi_vf1
from ohtake_surface_projection import display_simple_using_mayavi_
import sys

import numpy as np
from basic_types import check_vector4_vectorized, normalize_vector4_vectorized

#def adaptive_subdivision(self):
#    pass

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

        kij = np.arccos(mimj) / pipj  # radians?
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

        wi = 1.0 + c*ki
        # i_facet is facet (Centroid) index. j_facet is its neighbour facet (centroid). There are three j_facet for an i_facet.
        #print("w_i=", wi)
        return wi
    #
    c_ = c  # 2.0  # constant
    vertex_index = 1  # vertex
    #assert vertex_index >= 0
    #assert vertex_index < 
    umbrella_facets = neighbour_faces_of_vertex[vertex_index]  # A list of facets: The indices of faces that vertex vertex_index belongs to.
    print("umbrella_facets: ", umbrella_facets)
    #wa = np.zeros()
    w_list = []
    for i_facet in umbrella_facets:
        # neighbour facet i_facet of Vertex vertex_index
        print("i_facet", i_facet)
        #three_facets = filter(lambda idx: idx != i_facet, umbrella_facets)
        three_facets = faces_of_faces[i_facet, :]
        print(i_facet, three_facets)
        w = wi(i_facet, three_facets, c_)  # three_facets should be neighbours of the facet i_facet
        # The weight (based on curvature) of neighbour P_i (facet i.e. centroid),
        print("w_i, i=", i_facet, w)
        #w_list[i] = w
        w_list.append(w)
        #todo: sparse matrix: w[vi=vertex_index, f2=i_facet] = w
        #todo: store in ...
    #exit()
    print "w_list ",w_list
    #[1.9250638933714093, 2.0364604083536744, 1.4236331619142932, 3.4392903610759471, 5.4912745754508183, 3.2499307884393014, 5.0003861534703979]
    print("===============")
    #
    #w seems tobe calculated fine. next: store w_i and cache them for adaptive resampling, for which we need to normalise it across the neighbours.
    nfaces = centroids.shape[0]
    wi_total_array = np.zeros((nfaces,))
    for i_facet in range(nfaces):
        three_facets = faces_of_faces[i_facet, :]
        w = wi(i_facet, three_facets, c_)
        wi_total_array[i_facet] = w
    print wi_total_array
    # The weights are prepared. Now let's resample vertices

    vertex_index = 1
    #todo: umbrella_Facets = sparse matrix
    #umbrella_facets = np.array(neighbour_faces_of_vertex)  #empty

    umbrella_facets = np.array(neighbour_faces_of_vertex[vertex_index])  # empty
    print "umbrella_facets", umbrella_facets.shape, "****"
    assert np.allclose( wi_total_array[umbrella_facets] - np.array(w_list), 0)
    #return wi_total_array

    #

    def lift_verts(verts, centroids):
        new_verts = verts.copy()
        # assign these to a sparse matrix? and  do:  M = M/normalise(M); verts = M * verts
        for vertex_index in range(verts.shape[0]):
            umbrella_facets = np.array(neighbour_faces_of_vertex[vertex_index])
            w = wi_total_array[umbrella_facets]
            #w = w * 0 + 1
            w = w / np.sum(w)
            #print w / np.sum(w), w.shape
            new_verts[vertex_index, :] = \
                np.dot(w, centroids[umbrella_facets, 0:3])  # (n) * (n x 3)
        return new_verts

    return lift_verts(verts, centroids)


def compute_centroids(verts, facets):
    # see Mesh1::build_centroids
    #self.calculate_face_areas()
    expand = verts[facets,:]
    nfacets = facets.shape[0]
    assert expand.shape == (nfacets, 3, 3)
    assert np.allclose(verts[facets[:],:], expand)
    centroids3 = np.mean( verts[facets[:],:], axis=1)  # again
    centroids = np.concatenate( (centroids3, np.ones((nfacets,1))), axis=1)
    return centroids


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


def build_faces_of_faces(facets):
    """ builds lookup tables. The result if an array of nfaces x 3,
    containing the face index of neighbours of each face.
    Since each face has exactly three neighbours, the size of the result is n x 3."""
    from mesh_utils import make_edge_lookup_old
    (edges_of_faces, faces_of_edges, vertpairs_of_edges) = \
        make_edge_lookup_old(facets)

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


def process2_vertex_resampling_relaxation(verts, facets, iobj):
    centroids = compute_centroids(verts, facets)
    centroid_normals_normalized = compute_centroid_gradients(centroids, iobj, normalise=True)

    from mesh_utils import make_neighbour_faces_of_vertex
    neighbour_faces_of_vertex = make_neighbour_faces_of_vertex(facets)

    faces_of_faces = build_faces_of_faces(facets)

    new_verts = vertex_resampling(verts, neighbour_faces_of_vertex, faces_of_faces, centroids, centroid_normals_normalized, c=2.0)

    return new_verts, facets, centroids  # why does it return facets?


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
        separate=True, gradients_at=None, gradients_from_iobj=None, pointsizes=None):
    """Two separate panels"""

    print("Mayavi."); sys.stdout.flush()

    from mayavi import mlab

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
        mlab.triangular_mesh([vert[0] for vert in verts],
                         [vert[1] for vert in verts],
                         [vert[2] for vert in verts],
                         faces,
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
            mlab.points3d(c[:, 0], c[:, 1], c[:, 2], color=color_list[i], scale_factor=pointsizes[i] )
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

def compute_triangle_areas(verts, faces, return_normals=False):
    """ facet_normals: can contain NaN if the area is zero"""
    # see mesh1.py ::     def calculate_face_areas(self)
    DEGENERACY_THRESHOLD = 0.00001
    nfaces = faces.shape[0]
    expand = verts[faces, :]
    assert expand.shape == (nfaces, 3, 3)
    assert expand[:, 2, :].shape == (nfaces, 3)
    a = np.cross(
        expand[:, 1, :] - expand[:, 0, :],
        expand[:, 2, :] - expand[:, 0, :],
        axis=1)
    facet_areas = np.linalg.norm(a, axis=1, ord=2) / 2.0
    degenerates_count = len(facet_areas[facet_areas < DEGENERACY_THRESHOLD])
    facet_areas[facet_areas < DEGENERACY_THRESHOLD] = np.nan  # -1
    if degenerates_count > 0:
        print("degenerate triangles", degenerates_count)
    if not return_normals:
        #print "11"
        return facet_areas
    else:
        print facet_areas.shape
        assert facet_areas[:, np.newaxis].shape == (nfaces, 1)
        facet_normals = a / np.tile(facet_areas[:, np.newaxis], (1, 3)) / 2.0
        #print "22"
        return facet_areas, facet_normals


def process3_subdivide_example(fi, verts, facets, iobj):
    # fi = 100  # triangle T
    triangle = verts[facets[fi, :], :]  # numverts x 3
    assert triangle.shape == (3, 3)
    print triangle.shape
    tra = compute_triangle_areas(verts, facets)
    # A = compute_triangle_area(triangle)
    #l = tra.tolist()
    #l.sort()
    #print l  # .0022 ... 0.69  mm^2
    A = tra[fi]
    midp_matrix = np.array([
        [0.5,  0.5,  0],  # 3
        [0,  0.5,  0.5],  # 4
        [0.5,  0,  0.5]   # 5
        ]) .transpose()
    # (3=nv) x (3 num_verts of new nodes)
    VVV = triangle.transpose()  # 3 x (nv=3)
    m123 = np.dot(VVV, midp_matrix).transpose()  # (new verts) x 3
    #print m123
    new_verts = m123
    subdivision = triangle

    mini_verts = np.concatenate( (triangle, new_verts), axis=0)

    v012 = range(0, 3)  # facets[fi, :]
    v345 = range(3, 6)
    v012345 = np.array(v012 + v345)
    #print v012, v345
    mini_faces_l = [[0, 3, 5], [3, 1, 4], [5, 4, 2], [3, 4, 5]]  # 0,3,1,4,2,5
    #print v012345
    mini_faces = v012345[np.array(mini_faces_l)]
    #print mini_faces

    return mini_verts, mini_faces

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


def compute_facets_subdivision_curvatures(verts, facets, iobj):
    """ Deviation of Mesh from object gradients """

    #fi = 100  # triangle T
    #triangle = verts[facets[fi, :], :]  # numverts x 3
    #assert triangle.shape == (3, 3)
    #print triangle.shape
    facet_areas, facet_normals = compute_triangle_areas(verts, facets, return_normals=True)

    nf = facets.shape[0]
    assert facet_areas.shape == (nf,)
    assert facet_normals.shape == (nf, 3)
    #print "000000000s"

    #  #print np.logical_not(np.isnan(np.linalg.norm(facet_normals, axis=1) ))
    #  print verts[facets[np.isnan(np.linalg.norm(facet_normals, axis=1)), :],:]
    #  #print verts[faces[np.isnan(np.linalg.norm(facet_normals, axis=1)),:],:]
    #
    #  print verts[facets[np.isnan(np.linalg.norm(facet_normals, axis=1)), :],:]
    #  print facet_areas[np.isnan(np.linalg.norm(facet_normals, axis=1))]
    assert np.all(np.logical_not(np.isnan(facet_areas[np.logical_not(np.isnan(np.linalg.norm(facet_normals, axis=1)))])))

    #some edges are repeated

    #zero_normals = np.arange(facets.shape[0])[np.linalg.norm(facet_normals, axis=1) < 0.0000001]

    #print verts[facets[zero_normals, :], :]


    degenerate_faces = np.isnan(facet_areas)
    assert np.all(np.isnan(facet_areas[degenerate_faces]))
    assert np.all(np.logical_not(np.isnan(facet_areas[np.logical_not(degenerate_faces)])))
    assert np.all(np.isnan(facet_normals[degenerate_faces, :]))
    assert np.all(np.logical_not(np.isnan(facet_normals[np.logical_not(degenerate_faces),:])))
    #print len(degenerate_faces)


    #print zero_normals
    #print degenerate_faces
    #assert np.allclose(zero_normals - degenerate_faces, 0)
    #print facet_normals[zero_normals, :]

    #print facet_normals[zero_normals, :]
    #print "facet_areas", facet_areas[zero_normals] 
    #assert np.allclose(np.linalg.norm(facet_normals, axis=1)[np.logical_not(zero_normals)], 1.0)


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
        ])  # .transpose()


    e_array = np.zeros((nf,))
    for fi in range(nf):
        n = facet_normals[fi, :]  # n: (3,)

        triangle = verts[facets[fi, :], :]  # numverts x 3
        assert triangle.shape == (3, 3)
        #print triangle.shape

        assert triangle.shape == (3, 3)
        VVV = triangle  # (nv=3) x 3
        #print np.dot( centroidmaker_matrix, subdiv_vert_matrix)
        #exit()
        m0123 = np.dot( centroidmaker_matrix, np.dot(subdiv_vert_matrix, VVV) )
        assert m0123.shape == (4, 3)
        subdiv_centroids = m0123
        #print subdiv_centroids
        numsubdiv = 4

        subdiv_centroids4 = np.concatenate( (subdiv_centroids, np.ones((numsubdiv, 1))), axis=1)
        mm = - iobj.implicitGradient(subdiv_centroids4)[:, 0:3]
        assert mm.shape == (4, 3)
        nn = np.linalg.norm(mm, axis=1)
        mm = mm / np.tile(nn[:,np.newaxis], (1, 3))  # mm: 4 x 3
        mm = mm.transpose()  # 3x4
        e = facet_areas[fi] * np.sum(1. - np.abs(np.dot(n, mm))) / 4.  # sum(,x4)

        #assert np.all(np.dot(n, mm) > -0.0000001 ), "ingrown normal!"


        #e = np.sum(1 - np.abs(np.dot(n, mm)))   # sum(,x4)   #forgot the abs!
        e_array[fi] = e
        #if e<0:
        #    set_trace()


        if fi % 100 == 0:
            print fi, "\r", ;import sys; sys.stdout.flush()
    #print str(nf) + "   "
    #print e_array
    l = e_array[np.logical_not(np.isnan(e_array))].tolist()
    l.sort()
    print "curvature: min,max = ", l[0], l[-1]   # 3.80127650325e-08, 0.0240651184551
    bad_facets_count = np.sum(degenerate_faces)
    #assert bad_facets_count == 0
    return e_array, bad_facets_count

def process3(verts, facets, iobj, epsilon):  # , centroid_normals):
    e_array, bad_facets_count = compute_facets_subdivision_curvatures(verts, facets, iobj)

    nfaces = facets.shape[0]

    #epsilon = 0.001
    #epsilon = 4  # 0.34  # 1/8
    l = e_array.tolist()
    l.sort()
    print "e_array      ", l[:10], ".,.", l[-10:]
    assert len(l) == nfaces
    print "7/8 median=", l[int(nfaces*( 10./11. ))]
    #set_trace()

    #epsilon = 4. # 0.22
    #epsilon = 4.2
    #epsilon = 0.0001
    #epsilon = 1. / 10.
    #epsilon = 1. / 1.

    a = np.arange(nfaces)[ e_array > epsilon ]
    #a = np.arange(nfaces)[ e_array < epsilon ]
    print "a:", a
    print "count need subdivision", len(a)
    return a
"""
def _prepare_grid_old(rng):
    assert rng.size < 200
    if rng.size > 200:
        raise PolygonizationError(("Grid too large ( >200 ): ", rng.size))

    (xx, yy, zz) = np.meshgrid(rng, rng, rng)
    #xyz = np.mgrid( rng, rng, rng )
    #xyza = xyz.reshape((len(rng)**3, 3))
    #xyza = np.concat( ( np.expand_dims( xyza, axis=3 ), ones(len(rng)**3,1)  ), axis=3 )
    #assert xyza.shape == (len(rng), len(rng), len(rng), 4)

    #X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j] #??

    xyza = np.transpose(np.vstack([xx.ravel(), yy.ravel(), zz.ravel(), (xx*0+1).ravel()]))
    assert xyza.shape[1:] == (4,)

    if VERBOSE:
        print(xyza.shape)
        print("done alloc")
        sys.stdout.flush()

    return xyza
def make_grid_m(iobj, rng):
    if old:
        xyza = _prepare_grid_old(rng)
    else:
        xyza = _prepare_grid(rng)
    #slow_grid__dont_use()

    vgrid_v = iobj.implicitFunction(xyza)
    vgrid = np.reshape(vgrid_v, (len(rng), len(rng), len(rng)), order='C')

    if np.sum(vgrid_v > 0) == 0:
        raise PolygonizationError("The shape is empty. No interior points detected")
    if VERBOSE:
        print("interior points:", np.sum(vgrid_v > 0))
    return vgrid
"""
def visualise_normals_test():
    """ Visualised normals on a given example object"""
    from example_objects import make_example_vectorized
    #exname = "bowl_15_holes"  # "blend_example2_discs" "french_fries_vectorized" "cube_example"
    #exname = "blend_example2_discs" # 
    #exname ="ell_example1" #
    #exname = "first_csg"
    #exname = "bowl_15_holes"
    #iobj = make_example_vectorized("blend_example2_discs")
    ##(RANGE_MIN, RANGE_MAX, STEPSIZE) = (-20., 30., 1/1.)
    #(RANGE_MIN,RANGE_MAX, STEPSIZE) = (-1, +2, 0.1)

    iobj = make_example_vectorized("blend_example2_discs")
    (RANGE_MIN,RANGE_MAX, STEPSIZE) = (-3, +4, 0.2)

    #from example_objects import blend_example2_discs
    #iobj = blend_example2_discs(8.)
    #(RANGE_MIN, RANGE_MAX, STEPSIZE) = (-20., 30., 1/1.)




    #from stl_tests import make_mc_values_grid_mayavi

    from stl_tests import make_mc_values_grid
    gridvals = make_mc_values_grid(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE, old=False)
    verts, facets = vtk_mc(gridvals, (RANGE_MIN, RANGE_MAX, STEPSIZE) )
    print("MC calculated."); sys.stdout.flush()

    from mesh_utils import mesh_invariant
    mesh_invariant(facets)

    verts3 = verts.copy()

    ##############################
    #vv = verts[:, [1, 0, 2]]
    ##############################
    vv = verts


    print("Mayavi.");sys.stdout.flush()
    display_simple_using_mayavi_2( [ (vv, facets), ], pointcloud_list=[],
       mayavi_wireframe=False, opacity=[0.1, 1, 0.1], gradients_at=vv, separate=False, gradients_from_iobj=iobj, pointsizes=0.01)



def single_subdivision_demo():

    #set_trace()
    #dicesize = 16.
    #exname = "udice_vec"  # "blend_example2"
    #import example_objects
    #iobj = example_objects.make_example_vectorized(exname, dicesize)
    #(RANGE_MIN, RANGE_MAX, STEPSIZE) = (-22, +20., 0.8)

    #from example_objects import cyl4
    #iobj, (RANGE_MIN, RANGE_MAX, STEPSIZE) = \
    #    cyl4()
    #STEPSIZE = 1.

    #from example_objects import first_csg
    #iobj = \
    #    first_csg(8.)
    #(RANGE_MIN, RANGE_MAX, STEPSIZE) = (-20, 20, 1)

    from example_objects import blend_example2_discs
    #iobj = blend_example1(); (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-20/4., 20/4., 1/4.)

    iobj = blend_example2_discs(8.)
    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-20., 30., 1/1.)
    #curvature_epsilon = 1. / 4.
    #curvature_epsilon = 10000 # 1. / 40.
    #curvature_epsilon = 1. / 100   # larger==> less points
    #curvature_epsilon = 1. / 1000
    curvature_epsilon = 1. / 2000  # most points



    from example_objects import make_example_vectorized
    #exname = "bowl_15_holes"  # "blend_example2_discs" "french_fries_vectorized" "cube_example"
    #exname = "blend_example2_discs" #
    #exname ="ell_example1" #
    #exname = "first_csg"
    #exname = "bowl_15_holes"
    #(RANGE_MIN, RANGE_MAX, STEPSIZE) = (-20., 30., 1/1.)
    #iobj = make_example_vectorized("???")
    #(RANGE_MIN,RANGE_MAX, STEPSIZE) = (-1, +2, 0.2)

    #"rdice_vec"  too slow
    # screw3: terrible outcome
    (RANGE_MIN,RANGE_MAX, STEPSIZE) = (-3, +5, 0.2)

    #iobj = make_example_vectorized("screw3")
    #(RANGE_MIN,RANGE_MAX, STEPSIZE) = (-2, +2, 0.2)

    #iobj = make_example_vectorized("rdice_vec")
    #(RANGE_MIN,RANGE_MAX, STEPSIZE) = (-1.5, +1.5, 0.1)

    iobj = make_example_vectorized("ell_example1")
    (RANGE_MIN,RANGE_MAX, STEPSIZE) = (-3, +5, 0.2)


    from stl_tests import make_mc_values_grid
    gridvals = make_mc_values_grid(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE, old=False)
    verts, facets = vtk_mc(gridvals, (RANGE_MIN, RANGE_MAX, STEPSIZE) )
    print("MC calculated."); sys.stdout.flush()


    #from stl_tests import make_mc_mesh_scikit
    #verts2, faces2 = make_mc_mesh_scikit(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE)

    from mesh_utils import mesh_invariant
    mesh_invariant(facets)

    #display_simple_using_mayavi_vf1(verts2, faces2)
    #display_simple_using_mayavi_vf1(verts, facets)

    verts3 = verts.copy()
    #for i in range(2):
    #    verts3, faces3, centroids = process2_vertex_resampling_relaxation(verts3, facets, iobj)
    #    print("Process finished.");sys.stdout.flush()

    #print("Mayavi.");sys.stdout.flush()
    #display_simple_using_mayavi_vf1(verts3, faces3)
    #display_simple_using_mayavi_( [(verts, facets)], pointcloud_list=[verts3*1.0 , verts], opacity=0.4)
    #display_simple_using_mayavi_( [(verts3, facets)], pointcloud_list=[])

    #two windows
    #display_simple_using_mayavi_2( [(verts, facets), (verts3, facets)], pointcloud_list=[verts3, verts])

    #display_simple_using_mayavi_( [(verts, facets), (miniverts*1.2, minifaces)], pointcloud_list=[])
    if False:
        #pluck triangle # 100
        miniverts, minifaces = process3_subdivide_example(100, verts, facets, iobj)
        # pluck triangle # 0
        mv, mf = (verts[0:3, :], np.array([[0, 1, 2]]) )
        print("Mayavi.");sys.stdout.flush()

        #display_simple_using_mayavi_( [(mv,mf)], pointcloud_list=[])
        display_simple_using_mayavi_( [(mv*1.0, mf), (verts, facets), (miniverts*1.0, minifaces)], pointcloud_list=[],
            mayavi_wireframe=True, opacity=[1, 0.2, 1])

    #slow part
    chosen_faces = process3(verts, facets, iobj, curvature_epsilon)
    #bypass:
    #chosen_faces=np.array([0,1,2])

    print("Mayavi.+");sys.stdout.flush()
    #display_simple_using_mayavi_( [ (verts, facets), (verts, facets[chosen_faces, :]), ], pointcloud_list=[],
    #    mayavi_wireframe=False, opacity=[0.1, 1, 0.1])
    display_simple_using_mayavi_2( [ (verts, facets), (verts, facets[chosen_faces, :]), ], pointcloud_list=[],
       mayavi_wireframe=False, opacity=[0.2, 1, 0.1], gradients_at=None, separate=False, gradients_from_iobj=None,
       minmax=(RANGE_MIN,RANGE_MAX))
       # gradients_at=verts3, gradients_from_iobj=iobj)


def weighted_resampling_demo():

    #from example_objects import blend_example1, blend_example2_discs
    #iobj = blend_example1(); (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-20/4., 20/4., 1/4.)

    #from example_objects import blend_example2_discs
    #iobj = blend_example2_discs(8.)
    #(RANGE_MIN, RANGE_MAX, STEPSIZE) = (-20., 30., 2.)

    from example_objects import first_csg
    iobj = first_csg(8.)
    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-20.*2, 30.*2, 2.)


    if False:

        #curvature_epsilon = 1. / 4.
        #curvature_epsilon = 10000 # 1. / 40.
        #curvature_epsilon = 1. / 100   # larger==> less points
        #curvature_epsilon = 1. / 1000
        curvature_epsilon = 1. / 2000  # most points



        from example_objects import make_example_vectorized
        #exname = "bowl_15_holes"  # "blend_example2_discs" "french_fries_vectorized" "cube_example"
        #exname = "blend_example2_discs" #
        #exname ="ell_example1" #
        #exname = "first_csg"
        #exname = "bowl_15_holes"
        #(RANGE_MIN, RANGE_MAX, STEPSIZE) = (-20., 30., 1/1.)
        #iobj = make_example_vectorized("???")
        #(RANGE_MIN,RANGE_MAX, STEPSIZE) = (-1, +2, 0.2)

        #"rdice_vec"  too slow
        # screw3: terrible outcome
        (RANGE_MIN,RANGE_MAX, STEPSIZE) = (-3, +5, 0.2)

        #iobj = make_example_vectorized("screw3")
        #(RANGE_MIN,RANGE_MAX, STEPSIZE) = (-2, +2, 0.2)

        #iobj = make_example_vectorized("rdice_vec")
        #(RANGE_MIN,RANGE_MAX, STEPSIZE) = (-1.5, +1.5, 0.1)

        iobj = make_example_vectorized("ell_example1")
        (RANGE_MIN,RANGE_MAX, STEPSIZE) = (-3*3, +5*3, 0.2)


    from stl_tests import make_mc_values_grid
    gridvals = make_mc_values_grid(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE, old=False)
    verts, facets = vtk_mc(gridvals, (RANGE_MIN, RANGE_MAX, STEPSIZE) )
    print("MC calculated.")
    sys.stdout.flush()


    #pre-view
    #display_simple_using_mayavi_( [ (verts, facets), ], pointcloud_list=[], mayavi_wireframe=False, opacity=[1, 0.2, 0.7])


    #from stl_tests import make_mc_mesh_scikit
    #verts2, faces2 = make_mc_mesh_scikit(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE)

    from mesh_utils import mesh_invariant
    mesh_invariant(facets)

    #display_simple_using_mayavi_vf1(verts2, faces2)
    #display_simple_using_mayavi_vf1(verts, facets)

    RESAMPLING_ITERATIONS_COUNT = 5  # 2
    verts3_relaxed = verts.copy()
    for i in range(RESAMPLING_ITERATIONS_COUNT):
        verts3_relaxed, faces3, centroids = process2_vertex_resampling_relaxation(verts3_relaxed, facets, iobj)
        print("Process finished."); sys.stdout.flush()

    #print("Mayavi.");sys.stdout.flush()
    #display_simple_using_mayavi_vf1(verts3_relaxed, faces3)
    #display_simple_using_mayavi_( [(verts, facets)], pointcloud_list=[verts3_relaxed*1.0 , verts], opacity=0.4)
    #display_simple_using_mayavi_( [(verts3_relaxed, facets)], pointcloud_list=[])

    #two windows
    #display_simple_using_mayavi_2( [(verts, facets), (verts3_relaxed, facets)], pointcloud_list=[verts3_relaxed, verts])

    #display_simple_using_mayavi_( [(verts, facets), (miniverts*1.2, minifaces)], pointcloud_list=[])

    #pluck triangle # 100
    miniverts, minifaces = process3_subdivide_example(100, verts, facets, iobj)
    # pluck triangle # 0
    mv, mf = (verts[0:3, :], np.array([[0, 1, 2]]) )
    print("Mayavi.");sys.stdout.flush()

    #display_simple_using_mayavi_( [(mv,mf)], pointcloud_list=[])

    #WIREFRAME
    #display_simple_using_mayavi_( [(mv*1.0, mf), (verts, facets), (miniverts*1.0, minifaces)], pointcloud_list=[verts, verts3_relaxed],
    #    mayavi_wireframe=True, opacity=[1, 0.2, 1])

    #GOOD
    #display_simple_using_mayavi_2( [(mv*1.0, mf), (verts, facets), (miniverts*1.0, minifaces)], pointcloud_list=[verts, verts3_relaxed],
    #    mayavi_wireframe=False, opacity=[1, 0.2, 0.7], separate=False, pointsizes=[0.2, 0.6])
    display_simple_using_mayavi_2( [(verts3_relaxed, facets),  (mv*1.0, mf),  (miniverts*1.0, minifaces)], pointcloud_list=[verts, verts3_relaxed],
        mayavi_wireframe=False, opacity=[1, 0.2, 0.7], separate=False, pointsizes=[0.2, 0.6])


def subdivide_multiple_facets(verts_old, facets_old, tobe_subdivided_face_indices):

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

        new_facets[fi, :] = mini_faces[0, :]
        new_facets[new_facet_counter:(new_facet_counter+3), :] = mini_faces[1:(1+3), :]
        assert mini_faces.shape[0] == (1+3)
        #trace_subdivided_facets += range(new_facet_counter, (new_facet_counter+3))
        trace_subdivided_facets += range(new_facet_counter, (new_facet_counter+3)) + [fi]  # include the face which reuses the old face's index
        # trace_subdivided_facets will contain indices of faces

        new_facet_counter += 3


        #return mini_verts, mini_faces
        #numsubdiv = 4

        if fi % 100 == 0:
            print fi, "\r", ;import sys; sys.stdout.flush()

    print new_verts.shape[0], new_vertex_counter

    assert new_verts.shape[0] == new_vertex_counter
    assert new_facets.shape[0] == new_facet_counter
    print "v", provisional_new_verts_count+nverts_old, new_vertex_counter
    print "f", provisional_new_facets_count+nfaces_old, new_facet_counter
    assert provisional_new_verts_count+nverts_old == new_vertex_counter
    assert provisional_new_facets_count+nfaces_old == new_facet_counter
    assert len(trace_subdivided_facets) == 0 or np.max(np.array(trace_subdivided_facets)) < new_facet_counter
    return new_verts, new_facets



def process4_combine_both(verts, facets, iobj, epsilon, RESAMPLING_ITERATIONS_COUNT):
    # first part updates vertices only
    # second part updates faces only


    # from mesh_utils import mesh_invariant
    # mesh_invariant(facets) #isit really needed?

    #RESAMPLING_ITERATIONS_COUNT = 5  # 2
    verts3_relaxed = verts.copy()
    for i in range(RESAMPLING_ITERATIONS_COUNT):
        verts3_relaxed, faces3, centroids = process2_vertex_resampling_relaxation(verts3_relaxed, facets, iobj)
        print("Process finished.");sys.stdout.flush()
    # return verts3_relaxed

    verts = verts3_relaxed

    e_array, bad_facets_count = compute_facets_subdivision_curvatures(verts, facets, iobj)

    a = np.arange(facets.shape[0])[ e_array > epsilon ]
    #faces3_subdivided = subdivide_multiple_facets(faces, a)
    verts4_subdivided, faces3_subdivided = subdivide_multiple_facets(verts, facets, a)
    global trace_subdivided_facets  # third implicit output

    print "v,f=", verts.shape[0], facets.shape[0], "---->  v=", verts4_subdivided.shape[0], ", f=", faces3_subdivided.shape[0]
    print "max=", np.max(np.array(trace_subdivided_facets))

    # return faces3_subdivided
    #return verts3_relaxed, faces3_subdivided, trace_subdivided_facets
    return verts4_subdivided, faces3_subdivided, trace_subdivided_facets



def multiple_subdivisions_demo():
    #was: combination_actually_do

    from example_objects import blend_example2_discs


    iobj = blend_example2_discs(8.)
    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-20., 30., 1/1.)

    #curvature_epsilon = 1. / 2000.
    curvature_epsilon = 1. / 20000.  # still works nicely



    from example_objects import make_example_vectorized


    iobj = make_example_vectorized("ell_example1")
    (RANGE_MIN,RANGE_MAX, STEPSIZE) = (-3, +5, 0.2)


    from stl_tests import make_mc_values_grid
    gridvals = make_mc_values_grid(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE, old=False)
    verts, facets = vtk_mc(gridvals, (RANGE_MIN, RANGE_MAX, STEPSIZE) )
    print("MC calculated.")
    sys.stdout.flush()



    #from mesh_utils import mesh_invariant
    #mesh_invariant(facets)


    verts3 = verts.copy()  # copies twice


    # 0 => no vertex reaxation
    #faces3 =
    #verts3_relaxed, faces3_subdivided, trace_subdivided_facets = \
    new_verts3, new_faces3, trace_subdivided_facets = \
        process4_combine_both(verts, facets, iobj, curvature_epsilon, 0)


    #set_trace()
    chosen_faces = np.array(trace_subdivided_facets)
    assert chosen_faces.ndim == 1
    print "chosen_faces", chosen_faces.shape
    #print np.max(chosen_faces.ravel())
    print new_verts3.shape
    #print "chosen_faces", chosen_faces

    assert np.all(chosen_faces < new_faces3.shape[0])

    #print("*********")
    print new_verts3.shape
    print new_faces3.shape

    display_simple_using_mayavi_2( [ (new_verts3, new_faces3), (new_verts3, new_faces3[chosen_faces, :]), ], pointcloud_list=[],
       mayavi_wireframe=[False, True], opacity=[0.2, 1, 0.1], gradients_at=None, separate=False, gradients_from_iobj=None,
       minmax=(RANGE_MIN,RANGE_MAX))

    #display_simple_using_mayavi_2( [ (verts, facets), (new_verts3, new_faces3[chosen_faces, :]), ], pointcloud_list=[],
    #   mayavi_wireframe=False, opacity=[0.2, 1, 0.1], gradients_at=None, separate=False, gradients_from_iobj=None,
    #   minmax=(RANGE_MIN,RANGE_MAX))

    #display_simple_using_mayavi_2( [ (verts, facets), (new_verts3, new_faces3[chosen_faces, :]), ], pointcloud_list=[],
    #   mayavi_wireframe=False, opacity=[0.2, 1, 0.1], gradients_at=None, separate=False, gradients_from_iobj=None,
    #   minmax=(RANGE_MIN,RANGE_MAX))


def combination_actually_do():

    #1. / 2000/50. # 43K
    #curvature_epsilon = 1. / 2000/2.  # most points
    curvature_epsilon = 1. / 1000.
    VERTEX_RELAXATION_ITERATIONS_COUNT = 0
    SUBDIVISION_ITERATIONS_COUNT = 1 # 5+4

    from example_objects import make_example_vectorized
    iobj = make_example_vectorized("ell_example1")  # "bowl_15_holes" works too
    (RANGE_MIN,RANGE_MAX, STEPSIZE) = (-3, +5, 0.2)


    from stl_tests import make_mc_values_grid
    gridvals = make_mc_values_grid(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE, old=False)
    verts, facets = vtk_mc(gridvals, (RANGE_MIN, RANGE_MAX, STEPSIZE) )
    print("MC calculated.");sys.stdout.flush()

    old_verts, old_facets = verts, facets


    #from mesh_utils import mesh_invariant
    #mesh_invariant(facets)


    #new_verts3, new_facets3, trace_subdivided_facets = \
    #    process4_combine_both(verts, facets, iobj, curvature_epsilon, 0)

    for i in range(VERTEX_RELAXATION_ITERATIONS_COUNT):
        verts, facets_not_used, centroids = process2_vertex_resampling_relaxation(verts, facets, iobj)
        print("Vertex relaxation applied.");sys.stdout.flush()

    total_subdivided_facets = []
    for i in range(SUBDIVISION_ITERATIONS_COUNT):
        e_array, bad_facets_count = compute_facets_subdivision_curvatures(verts, facets, iobj)

        #print e_array

        #ohtake_belyaev_2.py:1122: RuntimeWarning: invalid value encountered in greater
        e_array[np.isnan(e_array)] = 0  # treat NaN curvatures as zero curvature => no subdivision
        #if np.any(np.isnan(e_array)):
        #    print "funny NaN values around."
        which_facets = np.arange(facets.shape[0])[ e_array > curvature_epsilon ]

        verts4_subdivided, facets3_subdivided = subdivide_multiple_facets(verts, facets, which_facets)
        global trace_subdivided_facets  # third implicit output
        #chosen_facet_indices = np.array(trace_subdivided_facets)
        verts, facets = verts4_subdivided, facets3_subdivided
        print("Subdivision applied.");sys.stdout.flush()

        total_subdivided_facets += trace_subdivided_facets  # old face indices remain valid

    chosen_facet_indices = np.array(total_subdivided_facets)

    if chosen_facet_indices.size == 0:
        chosen_subset_of_facets = np.zeros((0,), dtype=int)
    else:
        chosen_subset_of_facets = facets[chosen_facet_indices, :]
    display_simple_using_mayavi_2( [ (old_verts, old_facets), (verts, chosen_subset_of_facets), ], pointcloud_list=[],
       mayavi_wireframe=[False, True], opacity=[0.2, 1, 0.1], gradients_at=None, separate=False, gradients_from_iobj=None,
       minmax=(RANGE_MIN,RANGE_MAX))



if __name__ == '__main__':
    demo_choise = 5
    if demo_choise == 1:
        visualise_normals_test()   # visualise to check the gradients
    elif demo_choise == 2:
        single_subdivision_demo()  # just shows which ones subdivided
    elif demo_choise == 3:
        weighted_resampling_demo()
    elif demo_choise == 4:
        multiple_subdivisions_demo()  # nice demo! keep it
    elif demo_choise == 5:
        combination_actually_do()
    else:
        print "Error"
