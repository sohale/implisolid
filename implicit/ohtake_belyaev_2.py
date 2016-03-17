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


def process2(verts, facets, iobj):
    centroids = compute_centroids(verts, facets)
    centroid_normals_normalized = compute_centroid_gradients(centroids, iobj, normalise=True)

    from mesh_utils import make_neighbour_faces_of_vertex
    neighbour_faces_of_vertex = make_neighbour_faces_of_vertex(facets)

    faces_of_faces = build_faces_of_faces(facets)

    new_verts = vertex_resampling(verts, neighbour_faces_of_vertex, faces_of_faces, centroids, centroid_normals_normalized, c=2.0)

    return new_verts, facets, centroids


def visualise_gradients(mlab, pos, iobj):
    lm = 1  # STEPSIZE
    pos4 = np.concatenate((pos, np.ones((pos.shape[0],1))),axis=1)
    pnormals = - iobj.implicitGradient(pos4)
    pnormals = normalize_vector4_vectorized(pnormals) 
    check_vector4_vectorized(pos4)
    xyz = pos
    uvw = pnormals
    xx, yy, zz = xyz[:, 0], xyz[:, 1], xyz[:, 2]
    uu, vv, ww = uvw[:, 0], uvw[:, 1], uvw[:, 2]
    #ax.quiver
    #ax.quiver(xx, yy, zz,   uu, vv, ww,  length=np.abs(lm), arrow_length_ratio=0.3, alpha=0.3, pivot="tail")
    #arrow_length_ratio=   length=np.abs(lm)
    #pivot: tail | middle | tip
    #mlab.quiver3d(x_verts,y_verts,z_verts, UVW_normals[:,0],UVW_normals[:,1],UVW_normals[:,2],color=(0,0,0))
    mlab.quiver3d(xx, yy, zz, uu, vv, ww, color=(0, 0, 0))


def display_simple_using_mayavi_2(vf_list, pointcloud_list, minmax=(-1,1), mayavi_wireframe=False, opacity=1.0, 
        separate=True, gradients_at=None, gradients_from_iobj=None):
    """Two separate panels"""
    from mayavi import mlab

    if type(opacity) is list:
        opacities = opacity  # 1.0
    else:
        opacities = [opacity] + [0.2]*(len(vf_list)-1)  # 1.0, 0.2 #0.1


    for fi in range(len(vf_list)):
        if separate:
            mlab.figure()

        vf = vf_list[fi]
        verts, faces = vf

        mlab.triangular_mesh([vert[0] for vert in verts],
                         [vert[1] for vert in verts],
                         [vert[2] for vert in verts],faces,representation="surface" if not mayavi_wireframe else "wireframe",
                         opacity=opacities[fi], scale_factor = 100.0)
        #opacity = 0.2 #0.1


        color_list = [(1, 0, 0), (0, 0, 0), (1, 1, 0), (0, 0, 1), (0,1,0)]
        i = 0
        for c in pointcloud_list:
            #print c[:,0:3]
            mlab.points3d(c[:, 0], c[:, 1], c[:, 2], color=color_list[i], scale_factor=0.2)
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
        ampl = 2
        x = basic_types.make_random_vector_vectorized(n, ampl, 1, type="rand", normalize=False)
        v = iobj.implicitFunction(x)
        x_sel =  x[ v >= 0 , :]
        if x_sel.size ==0:
            print("No points")
            return
        ax.points3d(x_sel[:,0], x_sel[:,1], x_sel[:,2], color=(0,0,0), scale_factor=0.2)

    if gradients_at is not None:
        visualise_gradients(mlab, gradients_at, gradients_from_iobj)
    if gradients_from_iobj is not None:
        add_random_interior_points(mlab, gradients_from_iobj)

    mlab.show()
    return


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
        print "11"
        return facet_areas
    else:
        print facet_areas.shape
        assert facet_areas[:, np.newaxis].shape == (nfaces, 1)
        facet_normals = a / np.tile(facet_areas[:, np.newaxis], (1, 3)) / 2.0
        print "22"
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
    print m123
    new_verts = m123
    subdivision = triangle

    mini_verts = np.concatenate( (triangle, new_verts), axis=0)

    f012 = range(0, 3)  # facets[fi, :]
    f345 = range(3, 6)
    f012345 = np.array(f012 + f345)
    #print f012, f345
    mini_faces_l = [[0, 3, 5], [3, 1, 4], [5, 4, 2], [3, 4, 5]]  # 0,3,1,4,2,5
    #print f012345
    mini_faces = f012345[np.array(mini_faces_l)]
    #print mini_faces

    return mini_verts, mini_faces


#def subdivide_facet():

def subdivide_facets(verts, facets, iobj):
    """ Deviation of Mesh from object gradients """

    #fi = 100  # triangle T
    #triangle = verts[facets[fi, :], :]  # numverts x 3
    #assert triangle.shape == (3, 3)
    #print triangle.shape
    facet_areas, facet_normals = compute_triangle_areas(verts, facets, return_normals=True)

    nf = facets.shape[0]
    assert facet_areas.shape == (nf,)
    assert facet_normals.shape == (nf, 3)
    print "000000000s"

    #print np.logical_not(np.isnan(np.linalg.norm(facet_normals, axis=1) ))
    print verts[facets[np.isnan(np.linalg.norm(facet_normals, axis=1)), :],:]
    #print verts[faces[np.isnan(np.linalg.norm(facet_normals, axis=1)),:],:]

    print verts[facets[np.isnan(np.linalg.norm(facet_normals, axis=1)), :],:]
    print facet_areas[np.isnan(np.linalg.norm(facet_normals, axis=1))]
    assert np.all(np.logical_not(np.isnan(facet_areas[np.logical_not(np.isnan(np.linalg.norm(facet_normals, axis=1)))])))

    #some edges are repeated

    #zero_normals = np.arange(facets.shape[0])[np.linalg.norm(facet_normals, axis=1) < 0.0000001]

    #print verts[facets[zero_normals, :], :]


    degenerate_faces = np.isnan(facet_areas)
    assert np.all(np.isnan(facet_areas[degenerate_faces]))
    assert np.all(np.logical_not(np.isnan(facet_areas[np.logical_not(degenerate_faces)])))
    assert np.all(np.isnan(facet_normals[degenerate_faces, :]))
    assert np.all(np.logical_not(np.isnan(facet_normals[np.logical_not(degenerate_faces),:])))
    print len(degenerate_faces)


    #print zero_normals
    #print degenerate_faces
    #assert np.allclose(zero_normals - degenerate_faces, 0)
    #print facet_normals[zero_normals, :]

    #print facet_normals[zero_normals, :]
    #print "facet_areas", facet_areas[zero_normals] 
    #assert np.allclose(np.linalg.norm(facet_normals, axis=1)[np.logical_not(zero_normals)], 1.0)


    mc = np.array([
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
        #print np.dot( mc, subdiv_vert_matrix)
        #exit()
        m0123 = np.dot( mc, np.dot(subdiv_vert_matrix, VVV) )
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
    print str(nf) + "   "
    #print e_array
    l = e_array[np.logical_not(np.isnan(e_array))].tolist()
    l.sort()
    #print "e=", l[:10], " [...] ", l[-10:]
    print "e: min,max = ", l[0], l[-1]   # 3.80127650325e-08, 0.0240651184551
    return e_array

def process3(verts, facets, iobj, epsilon):  # , centroid_normals):
    e_array = subdivide_facets(verts, facets, iobj)

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
    iobj = make_example_vectorized("blend_example2_discs")
    #(RANGE_MIN, RANGE_MAX, STEPSIZE) = (-20., 30., 1/1.)
    (RANGE_MIN,RANGE_MAX, STEPSIZE) = (-1, +2, 0.1)

    #from stl_tests import make_mc_values_grid_mayavi

    from stl_tests import make_mc_values_grid
    gridvals = make_mc_values_grid(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE, old=False)
    verts, facets = vtk_mc(gridvals, (RANGE_MIN, RANGE_MAX, STEPSIZE) )
    print("MC calculated.")
    sys.stdout.flush()

    from mesh_utils import mesh_invariant
    mesh_invariant(facets)

    verts3 = verts.copy()

    print("Mayavi."); sys.stdout.flush()

    ##############################
    vv = verts[:, [1, 0, 2]]
    ##############################


    print("Mayavi.");sys.stdout.flush()
    display_simple_using_mayavi_2( [ (vv, facets), ], pointcloud_list=[],
       mayavi_wireframe=False, opacity=[0.1, 1, 0.1], gradients_at=vv, separate=False, gradients_from_iobj=iobj)



def ob2_test():

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

    iobj = make_example_vectorized("first_csg")
    (RANGE_MIN,RANGE_MAX, STEPSIZE) = (-3, +5, 0.2)


    from stl_tests import make_mc_values_grid
    gridvals = make_mc_values_grid(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE, old=False)
    verts, facets = vtk_mc(gridvals, (RANGE_MIN, RANGE_MAX, STEPSIZE) )
    print("MC calculated.")
    sys.stdout.flush()


    #from stl_tests import make_mc_mesh_scikit
    #verts2, faces2 = make_mc_mesh_scikit(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE)

    from mesh_utils import mesh_invariant
    mesh_invariant(facets)

    #display_simple_using_mayavi_vf1(verts2, faces2)
    #display_simple_using_mayavi_vf1(verts, facets)

    verts3 = verts.copy()
    #for i in range(2):
    #    verts3, faces3, centroids = process2(verts3, facets, iobj)
    #    print("Process finished.");sys.stdout.flush()

    print("Mayavi.");sys.stdout.flush()
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

    chosen = process3(verts, facets, iobj, curvature_epsilon)
    print("Mayavi.+");sys.stdout.flush()
    #display_simple_using_mayavi_( [ (verts, facets), (verts, facets[chosen, :]), ], pointcloud_list=[],
    #    mayavi_wireframe=False, opacity=[0.1, 1, 0.1])
    display_simple_using_mayavi_2( [ (verts, facets), (verts, facets[chosen, :]), ], pointcloud_list=[],
       mayavi_wireframe=False, opacity=[0.1, 1, 0.1], gradients_at=None, separate=False, gradients_from_iobj=None,
       minmax=None) 
       # gradients_at=verts3, gradients_from_iobj=iobj)


if __name__ == '__main__':
    ob2_test()
    # visualise_normals_test()   # visualise to check the gradients
