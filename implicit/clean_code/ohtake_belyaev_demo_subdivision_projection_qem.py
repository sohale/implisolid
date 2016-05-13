from ipdb import set_trace
import profile_support
from vtk_mc import vtk_mc

import sys
import math

import numpy as np
from basic_functions import check_vector3_vectorized, normalize_vector3_vectorized, normalize_vector4_vectorized

mesh_correction = False

def bisection_3_standard(iobj, p1, p2, f1, f2, MAX_ITER):
    TH1 = 0.001
    TH3 = 0.001

    assert p1.shape[0] == 1
    assert p2.shape[0] == 1

    assert f1 < 0
    assert f1*f2 < 0, "Opposite signs required"
    for j in range(MAX_ITER):
        if np.linalg.norm(p1-p2) < TH1:
            if j == 0:
                p3 = 0.5 * (p1 + p2)
            return p3, j
        p3 = 0.5 * (p1 + p2)
        f3 = iobj.implicitFunction(p3)
        if np.abs(f3) < TH3:
            return p3, j
        elif f1*f3 >= 0:
            p1 = p3
            f1 = f3
        else:
            p2 = p3
            f2 = f3
    # print "Convergence of the bisection did not happen"
    return None, MAX_ITER


# @profile
def bisection_prop_2(iobj, p1, p2, f1, f2, MAX_ITER):

    """ The proportional bisection method. See searchNearPoint1D() """
    # The following is used with MAX_ITER==10 in Ohtake.
    TH2 = 0.0001  # distance between p1,p2
    TH3 = 0.001   # for abs(f3)

    assert p1.shape[0] == 1
    assert p2.shape[0] == 1

    dt = p2 - p1

    assert f1 < 0
    assert f1*f2 < 0, "Opposite signs required"

    for j in range(MAX_ITER):

        assert f1 < 0
        assert f1*f2 < 0
        fp = f1/(f1-f2)
        p3 = p1 + dt * fp
        f3 = iobj.implicitFunction(p3)

        if np.abs(f3) < TH3:
            return True, p3, None, j

        # print f1, f2, " -> fp:", fp  # Sometimes (0.96) it searches too close to f2, and fp converges to 1
        if f3 < 0.:
            if VERBOSE:
                print "A",
            (p1, f1) = (p3, f3)
        else:
            if VERBOSE:
                print "B",
            (p2, f2) = (p3, f3)

        dt = p2 - p1
        # should be moved above

        if np.linalg.norm(dt) < TH2:
            return True, p3, None, j

    if VERBOSE:
        print "Convergence based on the proportional method did not happen"
    # return (False, p1, p2)
    p, ni = bisection_3_standard(iobj, p1, p2, f1, f2, MAX_ITER*4)
    if p is not None:
        if VERBOSE:
            print "Bisection converged "
        return True, p, None, 0  # very unlikely
    else:
        return False, p1, p2, -ni


VERBOSE = False

@profile
def search_near_ohtake(iobj, start_x, direction, lambda_val, MAX_ITER, max_dist):  # max_dist
    """Returns either the point, or None, if not found. lambda_val is the expected distance from the surface.
    The resommended value is half of the average edge length, but a better way is half of the MC'step size (because the expected error is half of the MC grid voxel size).
    Remeber: we may be on a totally irrelevant direction here.
    'direction' should be normalised. lambda*direction is used. Note that lambda is negated by default.
    Note: along_1d mode is in fact the same. It's just initialises direction=gradient(start_x).
    Does both searchNearPoint1D and searchNearPoint()
    :param direction: description
    @param np.array direction
    """
    # lambda should be ~ expected distance?  (that's why it should be half of the average edge size, i.e. half of the MC step size)

    TH1 = 0.001
    # MAX_ITER = 20
    # TH2_L = 0.00000001  # Used by Ohtake. But too small
    # TH2_L = 0.00001  #only used in the along_1d mode
    TH2_L = 0.1/2.
    if direction is not None:
        along_1d = True
    else:
        along_1d = False

    eval_count = 0

    p1 = start_x
    f1 = iobj.implicitFunction(p1)
    eval_count += 1
    f0 = f1

    if math.fabs(f1) < TH1:
        return p1

    p2 = p1  # no need actually

    negative_f1 = -1 if f1 < 0. else +1
    lambda_ = lambda_val * (-negative_f1)  # negation of input is bad practice
#    k = 0
    exit_A = False
    while True:
        # k += 1
        # print k
        max_iter_ = max(min(MAX_ITER, int(math.ceil(max_dist/math.fabs(lambda_)))+2), 2)
        assert max_iter_ >= 2
        # (C) jumps back here.
        for j in range(max_iter_):

            if not along_1d:
                direction = iobj.implicitGradient(start_x)
                dn = np.linalg.norm(direction)
                if dn > 0.0:  # 00000001:
                    direction = direction/dn
                else:
                    pass  # Finding is not going to happen. But it's fine.
            p2 = p2 + lambda_ * direction
            # p2[:, 3] = 1
            f2 = iobj.implicitFunction(p2)
            eval_count += 1

            if f1*f2 < 0.0:
                # (A)
                exit_A = True
                break
            p1 = p2

        else:
            # for loop ended becasue of MAX_ITER

            # (C): next iteration with adaptively decreased lambda_. Revert and start over again using a smaller lambda
            lambda_ = lambda_ / 2.
            # either quit:
            if np.abs(lambda_) < TH2_L:
                # print i, "(B)", eval_count
                print "(B)", eval_count
                return None   #
            # or back to start
            p1 = start_x
            assert f0 * f1 >= 0
            f1 = f0  # This was missing in Ohtake, because the sign of f1 is not expected to change. So I added the assert above.
            p2 = p1

            # restart the loop
            # (C).

        # (A)
        if exit_A:  # f1*f2 < 0.0:
            break
    # ( A)
    assert f1*f2 < 0.0
    if f1 > 0:
        (p1, p2) = (p2, p1)
        (f1, f2) = (f2, f1)

    assert f1 < 0
    assert f2 > 0
    assert math.fabs(f1) < 1000, str(p1)
    assert math.fabs(f2) < 1000, str(p2)

    converged, p1, p2, iter_log = bisection_prop_2(iobj, p1, p2, f1, f2, MAX_ITER/2)
    assert f1 < 0
    assert f2 > 0

    if converged:
        assert p2 is None
        return p1
    else:
        return None


@profile
def set_centers_on_surface_ohtake(iobj, centroids, average_edge):
    # here we consider that the max_dist is the average_edge and lambda = average_edge/2
    # new function who is a combination of sers_on_surface_ohtake and project_point_bidir_ohtake
    lambda_val = average_edge/2
    check_vector3_vectorized(centroids)
    # definition of the matrix that are gonna be used in the rest of the programm
    p1 = np.ndarray(centroids.shape)
    p2 = np.ndarray(centroids.shape)
    p3 = np.ndarray(centroids.shape)
    p = np.ndarray(centroids.shape)
    f1 = np.ndarray(centroids.shape[0])
    f2 = np.ndarray(centroids.shape[0])
    direction_3 = np.ndarray(centroids.shape)
    dn_3 = np.ndarray(centroids.shape[0])

    centroids_new = centroids.copy()
    max_iter = 20

    for i in range(centroids.shape[0]):

        # print i, "Trying to find in the first direction"
        p1[i, :] = search_near_ohtake(iobj, centroids_new[i, :].reshape((1, 3)), None, lambda_val, max_iter, lambda_val*2)
        if p1[i, :] is not None:  # make sure that p are found by the program and they respect the condition enforce by check_vector4_vectorized
            p1[i, :].reshape(1, 3)
            f1[i] = iobj.implicitFunction(p1[i, :].reshape(1, 3))

            # Mirror image: search the opposite way and accept only if it is closer than the best found.
            p2[i, :] = 2*centroids_new[i, :] - p1[i, :]  # p2 correspond to S in the paper
            f2[i] = iobj.implicitFunction(p2[i, :].reshape(1, 3))
            p[i, :] = p1[i, :]

            if f1[i]*f2[i] < 0:
                direction_3[i, :] = (centroids_new[i, :] - p1[i, :].copy())  # as in Ohtake
                dn_3[i] = np.linalg.norm(direction_3[i, :])
                if dn_3[i] > 0:  # dn>0.000000001:
                    direction_3[i, :] = direction_3[i, :]/dn_3[i]
                    p3 = search_near_ohtake(iobj, centroids_new[i, :].reshape(1, 3), direction_3[i, :].reshape(1, 3), lambda_val, max_iter, lambda_val*2)
            #        print i, "Trying to find in the opposite direction"

                if p3 is not None:
                    if np.linalg.norm(centroids_new[i, :] - p3) > np.linalg.norm(centroids_new[i, :] - p1[i, :]):
                        p[i, :] = p3
                        # else:
                        # p = p1

            if np.linalg.norm(centroids[i, :] - p[i, :]) <= average_edge:
                centroids[i, :] = p[i, :]


# def compute_triangle_areas(verts, faces, return_normals=False):
#     """ facet_normals: can contain NaN if the area is zero"""
#     # see mesh1.py ::     def calculate_face_areas(self)
#     DEGENERACY_THRESHOLD = 0.00001
#     nfaces = faces.shape[0]
#     expand = verts[faces, :]
#     assert expand.shape == (nfaces, 3, 3)
#     assert expand[:, 2, :].shape == (nfaces, 3)
#     a = np.cross(
#         expand[:, 1, :] - expand[:, 0, :],
#         expand[:, 2, :] - expand[:, 0, :],
#         axis=1)
#
#     facet_areas = np.linalg.norm(a, axis=1, ord=2) / 2.0
#     degenerates_count = len(facet_areas[facet_areas < DEGENERACY_THRESHOLD])
#
#     non_nul_indices = nfaces - degenerates_count  # the number of indices who have a facets with a non nul area
#     non_nul_facet_area = np.ndarray(non_nul_indices)
#     non_degenerate_faces = np.ndarray((non_nul_indices, 3), dtype=int)
#     non_degenerated_indice = []
#     k = 0
#     for i in range(nfaces):
#         if facet_areas[i] > DEGENERACY_THRESHOLD:
#             non_nul_facet_area[k] = facet_areas[i]
#             non_degenerate_faces[k] = faces[i]
#             non_degenerated_indice.append(i)
#
#             k += 1
#
#     new_expand = verts[faces[non_degenerated_indice], :]
#     new_a = np.cross(
#         new_expand[:, 1, :] - new_expand[:, 0, :],
#         new_expand[:, 2, :] - new_expand[:, 0, :],
#         axis=1)
#
#     if not return_normals:
#         return non_nul_facet_area
#     else:
#         print non_nul_facet_area.shape
#         assert non_nul_facet_area[:, np.newaxis].shape == (non_nul_indices, 1)
#         facet_normals = new_a / np.tile(non_nul_facet_area[:, np.newaxis], (1, 3)) / 2.0
#
#         return non_nul_facet_area, facet_normals, non_degenerate_faces

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
        facet_normals = a / np.tile(facet_areas[:, np.newaxis], (1, 3)) / 2.0
        return facet_areas, facet_normals

def compute_triangle_areas_old(verts, faces, return_normals=False):
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
        # print "11"
        return facet_areas
    else:
        print facet_areas.shape
        assert facet_areas[:, np.newaxis].shape == (nfaces, 1)
        facet_normals = a / np.tile(facet_areas[:, np.newaxis], (1, 3)) / 2.0
        # print "22"
        return facet_areas, facet_normals

def build_faces_of_faces(facets):
    """ builds lookup tables. The result if an array of nfaces x 3,
    containing the face index of neighbours of each face.
    Since each face has exactly three neighbours, the size of the result is n x 3."""
    from mesh_utils import make_edge_lookup_old
    (edges_of_faces, faces_of_edges, vertpairs_of_edges) = \
        make_edge_lookup_old(facets)

    nfaces = facets.shape[0]
    assert edges_of_faces.shape == (nfaces, 3)
    # f1 = faces_of_edges[edges_of_faces, 0]  # first face of all edges of all faces : nf x 3 -> nf
    # f2 = faces_of_edges[edges_of_faces, 1]  # second face of all edges of all faces: nf x 3 -> nf   #3 for 3 sides (edges)
    # one of the two (0:2) is equal to index=face. [index,3, 0:2 ]
    f12 = faces_of_edges[edges_of_faces, :]
    print f12.shape
    assert f12.shape == (nfaces, 3, 2)
    # print f12[0:5, :, :]
    # strip f12 from repeats of the same face as one of its neighbours (at each side: 3 sides). Hence, the size changes from (nfaces,3,2) to (nfaces,3)
    faceindex_eye = np.tile(np.arange(nfaces)[:, np.newaxis, np.newaxis], (1, 3, 2))
    assert faceindex_eye.shape == (nfaces, 3, 2)
    f12 = f12 + 1
    f12[f12 == faceindex_eye+1] = 0  # np.nan
    # print f12[0:5, :, :]
    assert np.allclose(np.prod(f12, axis=2), 0)  # check of each row has at least a zero
    f_uniq = np.sum(f12, axis=2)  # 0s will not be added
    # f_uniq[0:5,:]
    # f_uniq now contains the neighbour of each face. one face at each side of each face: nfaces x 3 -> face
    assert np.sum(f_uniq.ravel() == 0) == 0  # all faces have at least one neighbour at each side.
    return f_uniq - 1  # fix back the indices


def vertex_resampling(verts, faceslist_neighbours_of_vertex, faces_of_faces, centroids, centroid_normals, c=2.0):
    """ faceslist_neighbours_of_vertex: *** """

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
        # mimj[mimj>1.0] = 1.0
        if mimj > 1.0:
            mimj = 1.0
        pipj = np.linalg.norm(pi - pj)
        if pipj == 0:
            return 0
        # print "pipj ", pipj, "  arccos= ",np.arccos(mimj)/np.pi*180 #why is it zero??
        assert pipj == np.linalg.norm(pi - pj, ord=2)

        kij = np.arccos(mimj) / pipj  # radians?
        return kij

    def wi(i_facet, ja_facets, c):
        """
        Returns the weight of a facet i_facet.
        Adds kij of all centroids of the neighbour facets.
        ja_facets = list of centroid indices (face index).
        i_facet is a face index. """

        assert i_facet not in ja_facets
        assert len(ja_facets) == 3

        ki = 0
        for j_facet in ja_facets:
            ki += kij(i_facet, j_facet)
        wi = 1.0 + c*ki

        return wi
    #
    c_ = c  # 2.0  # constant
    vertex_index = 1  # vertex

    umbrella_facets = faceslist_neighbours_of_vertex[vertex_index]  # A list of facets: The indices of faces that vertex vertex_index belongs to.
    print("umbrella_facets: ", umbrella_facets)
    # wa = np.zeros()
    w_list = []
    for i_facet in umbrella_facets:

        print("i_facet", i_facet)
        three_facets = faces_of_faces[i_facet, :]
        print(i_facet, three_facets)
        w = wi(i_facet, three_facets, c_)  # three_facets should be neighbours of the facet i_facet
        # The weight (based on curvature) of neighbour P_i (facet i.e. centroid),
        print("w_i, i=", i_facet, w)
        w_list.append(w)
        # todo: sparse matrix: w[vi=vertex_index, f2=i_facet] = w
        # todo: store in ...
    # exit()
    print "w_list ", w_list

    print("===============")
    #
    # w seems tobe calculated fine. next: store w_i and cache them for adaptive resampling, for which we need to normalise it across the neighbours.
    nfaces = centroids.shape[0]
    wi_total_array = np.zeros((nfaces,))
    for i_facet in range(nfaces):
        three_facets = faces_of_faces[i_facet, :]
        w = wi(i_facet, three_facets, c_)
        wi_total_array[i_facet] = w
    print wi_total_array
    # The weights are prepared. Now let's resample vertices

    vertex_index = 1

    umbrella_facets = np.array(faceslist_neighbours_of_vertex[vertex_index])  # empty
    print "umbrella_facets", umbrella_facets.shape, "****"
    assert np.allclose(wi_total_array[umbrella_facets] - np.array(w_list), 0)

    def lift_verts(verts, centroids):
        new_verts = verts.copy()
        # assign these to a sparse matrix? and  do:  M = M/normalise(M); verts = M * verts
        for vertex_index in range(verts.shape[0]):
            umbrella_facets = np.array(faceslist_neighbours_of_vertex[vertex_index])
            w = wi_total_array[umbrella_facets]
            # w = w * 0 + 1
            w = w / np.sum(w)
            # print w / np.sum(w), w.shape
            new_verts[vertex_index, :] = \
                np.dot(w, centroids[umbrella_facets, 0:3])  # (n) * (n x 3)
        return new_verts

    return lift_verts(verts, centroids)

# evaluate_centroid_gradients


def compute_centroid_gradients(centroids, iobj, normalise=True):
    centroids = centroids[:, :3]
    assert centroids is not None

    check_vector3_vectorized(centroids)
    centroid_gradients = iobj.implicitGradient(centroids)
    assert not np.any(np.isnan(centroid_gradients))
    assert not np.any(np.isinf(centroid_gradients))
    if normalise:
        centroid_normals = normalize_vector3_vectorized(centroid_gradients)
        return centroid_normals
    else:
        return centroid_gradients


def visualise_gradients(mlab, pos, iobj, arrow_size):
    lm = arrow_size  # 1.  # STEPSIZE
    pos3 = pos
#    pos3 = np.concatenate((pos, np.ones((pos.shape[0],1))),axis=1)
    pnormals = - iobj.implicitGradient(pos3)
    pnormals = normalize_vector4_vectorized(pnormals)
    check_vector3_vectorized(pos3)
    xyz = pos3
    uvw = pnormals[:, 0:3] / 2.
    xx, yy, zz = xyz[:, 0], xyz[:, 1], xyz[:, 2]
    uu, vv, ww = uvw[:, 0], uvw[:, 1], uvw[:, 2]
    mlab.quiver3d(xx, yy, zz, uu, vv, ww, color=(0, 0, 0), scale_factor=np.abs(lm), line_width=0.5)


def display_simple_using_mayavi_2(vf_list, pointcloud_list, minmax=(-1, 1), mayavi_wireframe=False, opacity=1.0, separate=True,
 gradients_at=None, gradients_from_iobj=None, pointsizes=None, pointcloud_opacity=1.):
    """Two separate panels"""

    print("Mayavi.")
    sys.stdout.flush()

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
                         opacity=opacities[fi], scale_factor=100.0)
        # opacity = 0.2 #0.1
        # allpoints are plottedon all panels?
        color_list = [(1, 0, 0), (0, 0, 0), (1, 1, 0), (0, 0, 1), (0, 1, 0)]
        i = 0
        for c in pointcloud_list:
            # if separate:
            #    if i != fi:
            #        continue
            # print c[:,0:3]
            mlab.points3d(c[:, 0], c[:, 1], c[:, 2], color=color_list[i], scale_factor=pointsizes[i], opacity=pointcloud_opacity)
            i += 1
        del i

        if minmax is not None:
            (RANGE_MIN, RANGE_MAX) = minmax
            x = np.linspace(RANGE_MIN, RANGE_MAX, 2).reshape(2, 1)
            y = np.zeros((2, 1))
            z = np.zeros((2, 1))

            mlab.plot3d(x, y, z, line_width=3, name="x-axis")
            mlab.plot3d(y, x, z, line_width=3, name="y-axis")
            mlab.plot3d(z, y, x, line_width=3, name="z-axis")

            mlab.text3d(RANGE_MAX, 0, 0, "x", scale=0.3)
            mlab.text3d(0, RANGE_MAX, 0, "y", scale=0.3)
            mlab.text3d(0, 0, RANGE_MAX, "z", scale=0.3)
            mlab.text3d(RANGE_MIN, 0, 0, "-x", scale=0.3)
            mlab.text3d(0, RANGE_MIN, 0, "-y", scale=0.3)
            mlab.text3d(0, 0, RANGE_MIN, "-z", scale=0.3)

    def add_random_interior_points(ax, iobj, avg_edge_len):
        """ Adding random points """
        n = 10000
        import basic_functions

        ampl = avg_edge_len
        # ampl = 2
        x = basic_functions.make_random_vector_vectorized(n, ampl, 1, type="rand", normalize=False)
        v = iobj.implicitFunction(x)
        x_sel = x[v >= 0, :]
        if x_sel.size == 0:
            print("No points")
            return
        ax.points3d(x_sel[:, 0], x_sel[:, 1], x_sel[:, 2], color=(0, 0, 0), scale_factor=0.2)

    if gradients_at is not None:
        verts1, faces1 = vf_list[0]
        avg_edge_len = compute_average_edge_length(verts, faces)
        visualise_gradients(mlab, gradients_at, gradients_from_iobj, avg_edge_len / 20.)
    if gradients_from_iobj is not None:
        add_random_interior_points(mlab, gradients_from_iobj, avg_edge_len)

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
        e1 = np.linalg.norm(expand[:, i1, :] - expand[:, i2, :], axis=1)
        ea_sum += np.mean(e1)

    return ea_sum / 3.

# @profile


def get_A_b(vertex_id, nlist_numpy, centroids, centroid_gradients, qem_origin):

    nai = nlist_numpy
    center_array = centroids[nai, :]

    # note some centers may not be projected successfully in the previous step
    not_projected_successfully = np.isnan(center_array[:].ravel())
    if np.any(not_projected_successfully):
        pass

    normals = centroid_gradients[nai, :]  # why do we have repeats??
    # note : not normalised. But it works.

    norms = np.linalg.norm(normals, ord=2, axis=1)
    # can be either 0, 1 or Nan
    if np.any(norms < 0.000001):  # can be exactly 0.0
        print("Error: bad normal", normals)

    # TH_N = 0.0000001  # 0.000001 = I millions of millimeter = 1 nanometer
    # can be 0,0,0, inf, nonsharp, degenerate, ...
    # degenerate_normals = np.logical_or(np.isnan( np.sum(normals, axis=1)), norms < TH_N )
    assert not np.any(np.isnan(normals))
    assert not np.any(np.isinf(normals))


#    x0 = np.zeros((3, 1))

    assert normals.shape[1] == 3
    # normals = normals   # (4)x4
    # grad = Ax+b
    A = np.zeros((3, 3))
    b = np.zeros((3, 1))
    # assert len(center_array) == len(normals)
    assert normals.shape == center_array.shape
    n_i = normals[:, :, np.newaxis]
    p_i = center_array[:, :, np.newaxis]

    A = np.dot(np.reshape(n_i, (normals.shape[0], 3)).T, np.reshape(n_i, (normals.shape[0], 3)))

    for i in range(normals.shape[0]):

        assert n_i[i].shape == (3, 1)
        nnt = np.dot(n_i[i], np.transpose(n_i[i]))

        assert nnt.shape == (3, 3)

        assert p_i[i].shape == (3, 1)
        b += -np.dot(nnt, p_i[i] - qem_origin)

    return A, b

# @profile


def vertices_apply_qem3(verts, facets, centroids, vertex_neighbours_list, centroid_gradients):
    # based on quadratic_optimise_vertices(self, alpha=1.0)
    assert centroids is not None
    assert vertex_neighbours_list is not None
    assert centroid_gradients is not None

    # alpha = 1.0
    nvert = verts.shape[0]
    assert nvert == len(vertex_neighbours_list)

    result_verts_ranks = np.zeros((nvert,), dtype=int)
    assert verts.shape == (nvert, 3)
    new_verts = np.zeros((nvert, 3))

    for vertex_id in range(nvert):

        vi = vertex_id
        nlist = vertex_neighbours_list[vertex_id]
        nai = np.array(nlist)
        qem_origin = verts[vertex_id, :].reshape(3, 1)*0
        A, b = get_A_b(vi, nai, centroids, centroid_gradients, qem_origin)
    #    A, b = get_A_b(vi, nai, centroids, centroid_gradients)

        u, s, v = np.linalg.svd(A)
        assert np.allclose(A, np.dot(u, np.dot(np.diag(s), v)))
        assert s[0] == np.max(s)

    #    tau = 10. **0.37
        tau = 10. ** 3.
        s[s / s[0] < 1.0/tau] = 0
        # print(s , s[0] , tau)
        rank = np.sum(s / s[0] > 1.0/tau)
    #    print "rank:", rank

        assert np.all(s[:rank]/s[0] >= 1.0/tau)

        x = verts[vi, 0:3, np.newaxis]
        assert x.shape == (3, 1)

        y = np.dot(v, x).copy()
        utb = np.dot(-np.transpose(u), b)

        for i in range(rank):
            assert np.dot(-np.transpose(u), b).shape == (3, 1)
            # print s[i] , 1.0/tau
            # assert s[i] >= 1.0/tau #fails when s[0] is small
            assert s[i]/s[0] >= 1.0/tau
            y[i] = utb[i] / s[i]
        new_x = np.dot(np.transpose(v), y)
        # print(x.ravel(), " -> ", new_x.ravel())
        # print("    delta=", (new_x - x).ravel())

        if not s[0] > 0.000001:
            print("Warning! sigma_1 == 0")
            print(s)
            print("A", A)

            result_verts_ranks[vi] = 0

        assert x.shape == (3, 1)

        new_verts[vi, 0:3] = new_x[:, 0]

        if not np.all(np.abs(utb.ravel()[rank:]) < 0.0001):
            # print("s", s.ravel()/s[0], "   utb", utb.ravel()/s[0])
            pass
        result_verts_ranks[vi] = rank

        # exit()
    print("max rank = ", np.max(result_verts_ranks))
    print("min rank = ", np.min(result_verts_ranks))
    if not np.min(result_verts_ranks) >= 1:
        print("Warning: assertion: np.min(result_verts_ranks) >= 1 failed.")

    if False:
        assert np.min(result_verts_ranks) >= 1
    return new_verts


import mesh_utils


def compute_centroids(verts, facets):
    # see Mesh1::build_centroids
    # self.calculate_face_areas()
    expand = verts[facets, :]
    nfacets = facets.shape[0]
    assert expand.shape == (nfacets, 3, 3)
    assert np.allclose(verts[facets[:], :], expand)
    centroids = np.mean(verts[facets[:], :], axis=1)  # again
    return centroids


def process2_vertex_resampling_relaxation(verts, facets, iobj):
    assert not np.any(np.isnan(verts))
    centroids = compute_centroids(verts, facets)
    centroid_normals_normalized = compute_centroid_gradients(centroids, iobj, normalise=True)

    from mesh_utils import make_neighbour_faces_of_vertex
    faceslist_neighbours_of_vertex = make_neighbour_faces_of_vertex(facets)
#    import ipdb; ipdb.set_trace()
    faces_of_faces = build_faces_of_faces(facets)
    # print faces_of_faces.shape, "*x3"

    new_verts = vertex_resampling(verts, faceslist_neighbours_of_vertex, faces_of_faces, centroids, centroid_normals_normalized, c=2.0)

    return new_verts, facets, centroids  # why does it return facets?


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
        [1., 0., 0.],  # 0
        [0., 1., 0.],  # 1
        [0., 0., 1.],  # 2

        [0.5/(1.+DIP), 0.5/(1.+DIP), DIP/(1.+DIP)],  # 3
        [DIP/(1.+DIP), 0.5/(1.+DIP), 0.5/(1.+DIP)],  # 4
        [0.5/(1.+DIP), DIP/(1.+DIP), 0.5/(1.+DIP)]   # 5
        ])  # .transpose()

    global trace_subdivided_facets
    trace_subdivided_facets = []

    provisional_new_verts_count = 3*len(tobe_subdivided_face_indices)
    provisional_new_facets_count = 3*len(tobe_subdivided_face_indices)
    nverts_old = verts_old.shape[0]
    nfaces_old = facets_old.shape[0]
    new_verts = np.zeros((nverts_old+provisional_new_verts_count, 3), dtype=float)
    new_facets = np.zeros((nfaces_old+provisional_new_facets_count, 3), dtype=int)
    new_verts[:nverts_old, :] = verts_old
    new_facets[:nfaces_old, :] = facets_old

    # on number of added vertices:
    # problem: there may be repeated (Redundant) vertices. (as well as T-junctions)
    # also later check for faces with repeated edges. (which can be another cause of null normals)

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
        # subdiv_centroids = m0123

        vxyz_0123 = np.dot(subdiv_vert_matrix, VVV)  # not efficient
        assert vxyz_0123.shape == (6, 3)

        v012 = facets_old[fi, :].tolist()  # range(0, 3)  #
        v345 = range(new_vertex_counter, new_vertex_counter+3)

        v345_xyz = vxyz_0123[3:6, :]  # only pick the new ones

        assert len(v345) == 3
        new_verts[(new_vertex_counter):(new_vertex_counter+3), :] = v345_xyz

        new_vertex_counter += 3

        # facet's vertex indices
        v012345 = np.array(v012 + v345, dtype=int)
        mini_faces_l = [[0, 3, 5], [3, 1, 4], [5, 4, 2], [3, 4, 5]]  # 0,3,1,4,2,5
        mini_faces = v012345[np.array(mini_faces_l, dtype=int)]

        new_facets[fi, :] = mini_faces[0, :]
        new_facets[new_facet_counter:(new_facet_counter+3), :] = mini_faces[1:(1+3), :]
        assert mini_faces.shape[0] == (1+3)
        trace_subdivided_facets += range(new_facet_counter, (new_facet_counter+3)) + [fi]  # include the face which reuses the old face's index
        # trace_subdivided_facets will contain indices of faces

        new_facet_counter += 3

        if fi % 100 == 0:
            print fi, "\r",
            import sys
            sys.stdout.flush()

    print new_verts.shape[0], new_vertex_counter

    assert new_verts.shape[0] == new_vertex_counter
    assert new_facets.shape[0] == new_facet_counter
    print "v", provisional_new_verts_count+nverts_old, new_vertex_counter
    print "f", provisional_new_facets_count+nfaces_old, new_facet_counter
    assert provisional_new_verts_count+nverts_old == new_vertex_counter
    assert provisional_new_facets_count+nfaces_old == new_facet_counter
    assert len(trace_subdivided_facets) == 0 or np.max(np.array(trace_subdivided_facets)) < new_facet_counter
    return new_verts, new_facets


def compute_facets_subdivision_curvatures(verts, facets, iobj, curvature_epsilon):
    """ Deviation of Mesh from object gradients """

    facet_areas, facet_normals = compute_triangle_areas_old(verts, facets, return_normals=True)

    nf = facets.shape[0]
    assert facet_areas.shape == (nf,)
    assert facet_normals.shape == (nf, 3)

    assert np.all(np.logical_not(np.isnan(facet_areas[np.logical_not(np.isnan(np.linalg.norm(facet_normals, axis=1)))])))

    degenerate_faces = np.isnan(facet_areas)

    # assert np.all(np.isnan(facet_areas[degenerate_faces]))
    # assert np.all(np.logical_not(np.isnan(facet_areas[np.logical_not(degenerate_faces)])))
    # assert np.all(np.isnan(facet_normals[degenerate_faces, :]))
    # assert np.all(np.logical_not(np.isnan(facet_normals[np.logical_not(degenerate_faces), :])))

    centroidmaker_matrix = np.array([
        [1, 0, 0, 1, 0, 1],  # 035
        [0, 1, 0, 1, 1, 0],  # 314
        [0, 0, 1, 0, 1, 1],  # 542
        [0, 0, 0, 1, 1, 1],  # 345
        ]) / 3.

    subdiv_vert_matrix = np.array([
        [1., 0., 0.],  # 0
        [0., 1., 0.],  # 1
        [0., 0., 1.],  # 2

        [0.5, 0.5, 0],  # 3
        [0, 0.5, 0.5],  # 4
        [0.5, 0, 0.5]   # 5
        ])  # .transpose()

    #e_array = np.zeros((nf,))
    e_array = np.ndarray(nf)

    for fi in range(nf):
        if degenerate_faces[fi]:
            e_array[fi] = 0
            continue

        n = facet_normals[fi, :]  # n: (3,)

        triangle = verts[facets[fi, :], :]  # numverts x 3
        assert triangle.shape == (3, 3)
        # print triangle.shape

        assert triangle.shape == (3, 3)
        VVV = triangle  # (nv=3) x 3
        # print np.dot( centroidmaker_matrix, subdiv_vert_matrix)
        # exit()
        m0123 = np.dot(centroidmaker_matrix, np.dot(subdiv_vert_matrix, VVV))
        assert m0123.shape == (4, 3)
        subdiv_centroids = m0123
        # print subdiv_centroids

        mm = - iobj.implicitGradient(subdiv_centroids)
        assert mm.shape == (4, 3)
        nn = np.linalg.norm(mm, axis=1)
        nn_tile = np.tile(nn[:, np.newaxis], (1, 3))
        nn_tile[nn_tile < 0.00000001] = 100000.
    #    mm = mm / np.tile(nn[:,np.newaxis], (1, 3))  # mm: 4 x 3
        mm = mm/nn_tile
        mm = mm.transpose()  # 3x4
        if np.any(np.isnan(n)):
            e = 0.
        else:
            e = facet_areas[fi] * np.sum(1. - np.abs(np.dot(n, mm))) / 4.  # sum(,x4)

        # assert np.all(np.dot(n, mm) > -0.0000001 ), "ingrown normal!"

        # e = np.sum(1 - np.abs(np.dot(n, mm)))   # sum(,x4)   #forgot the abs!
        e_array[fi] = e
        # if e<0:
        #    set_trace()
        if fi % 100 == 0:
            print fi, "\r", ;import sys; sys.stdout.flush()

    # for i in range(nf):
    #     print i, e_array[i]

    num_subivision = len(e_array[e_array > curvature_epsilon])
    need_subidivision = np.ndarray(num_subivision)

    k = 0
    for i in range(nf):
        if e_array[i] > curvature_epsilon:
            need_subidivision[k] = i
            k += 1

    # print str(nf) + "   "
    # print e_array
    l = e_array[np.logical_not(np.isnan(e_array))].tolist()
    l.sort()
    print "curvature: min,max = ", l[0], l[-1]   # 3.80127650325e-08, 0.0240651184551

    return e_array, need_subidivision


# delete some artefacts dues in the sharps part of the mesh
def comparison_verts_new_verts(old_verts, new_verts):
    THL = 10. ** (0.001)  # find experimentally
    nverts = old_verts.shape[0]
    assert old_verts.shape == new_verts.shape

    for i in range(nverts):
        if np.abs(np.linalg.norm(old_verts[i, :] - new_verts[i, :])) > THL:
            new_verts[i, :] = old_verts[i, :]

    return new_verts

@profile
def demo_combination_plus_qem():
    """ Now with QEM """
    curvature_epsilon = 1. / 1000.  # a>eps  1/a > 1/eps = 2000
    # curvature_epsilon = 1. / 10000.
    VERTEX_RELAXATION_ITERATIONS_COUNT = 0
    SUBDIVISION_ITERATIONS_COUNT = 1  # 2  # 5+4

    from example_objects import make_example_vectorized
    object_name = "rdice_vec"  # "sphere_example" #or "rcube_vec" work well #"ell_example1"#"cube_with_cylinders"#"ell_example1"  " #"rdice_vec" #"cube_example"
    iobj = make_example_vectorized(object_name)

    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-3, +5, 0.2)
    if object_name == "cube_with_cylinders" or object_name == "french_fries" or object_name == "rdice_vec" or object_name == "rods" or object_name == "bowl_15_holes":
        VERTEX_RELAXATION_ITERATIONS_COUNT = 1

    if object_name == "cyl4":
        (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-32 / 2, +32 / 2, 1.92 / 4.0)

    elif object_name == "french_fries" or object_name == "rods":
        (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-3, +5, 0.11)

    elif object_name == "bowl_15_holes":
        (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-3, +5, 0.26)

    elif object_name == "cyl3":
        (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-32 / 2, +32 / 2, 1.92 / 4.0)

    elif object_name == "cyl1":
        (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-16, +32, 1.92 * 0.2 * 10 / 2.0)

    # import vectorized, example_objects
    # c2 = vectorized.UnitCube1(1.)
    # def rotate_scale_(iobj, scale, center, angle=0.):
    #     ns = vectorized
    #     import numpy
    #     m = numpy.eye(4)
    #     m[0,0] = 0.1
    #     iobj = ns.Transformed(iobj, m=m)
    #     iobj  \
    #         .resize(scale) \
    #         .move(center[0], center[1], center[2])
    #     if angle != 0.:
    #         iobj.rotate(angle, along=make_vector4(1, 1, 1), units="deg")
    #     return iobj
    #
    # c2 = rotate_scale_(c2, 2., [1,1,1])
    # iobj = vectorized.CrispUnion( example_objects.rcube_vec(1.), c2 )

    from stl_tests import make_mc_values_grid
    gridvals = make_mc_values_grid(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE, old=False)
    verts, facets = vtk_mc(gridvals, (RANGE_MIN, RANGE_MAX, STEPSIZE))
    print("MC calculated.");sys.stdout.flush()

    old_verts, old_facets = verts, facets
    #
    # display_simple_using_mayavi_2( [(verts, facets),(verts, facets), ],
    #    pointcloud_list=[],
    #    mayavi_wireframe=[False, True,], opacity=[1, 1, 0.9], gradients_at=None, separate=False, gradients_from_iobj=None,
    #    minmax=(RANGE_MIN,RANGE_MAX)  )
    # exit()

    total_subdivided_facets = []
    for i in range(SUBDIVISION_ITERATIONS_COUNT):
        e_array, which_facets = compute_facets_subdivision_curvatures(verts, facets, iobj, curvature_epsilon)

        # which_facets = np.arange(facets.shape[0])[e_array > curvature_epsilon]
        print curvature_epsilon, which_facets.shape

        verts4_subdivided, facets3_subdivided = subdivide_multiple_facets(verts, facets, which_facets)
        global trace_subdivided_facets  # third implicit output

        verts, facets = verts4_subdivided, facets3_subdivided
        print("Subdivision applied.");sys.stdout.flush()

        # total_subdivided_facets += trace_subdivided_facets  # old face indices remain valid


    average_edge = compute_average_edge_length(verts, facets)

    old_centroids = np.mean(verts[facets[:], :], axis=1)

    new_centroids = old_centroids.copy()
    # new_centroids = set_centers_on_surface_ohtake_not_correct_points(iobj, new_centroids, average_edge)
    # new_centroids is the output
    # set_centers_on_surface_ohtake(iobj, new_centroids, average_edge*2)
    set_centers_on_surface_ohtake(iobj, new_centroids, average_edge)
    # display_simple_using_mayavi_2( [(verts, facets),(verts, facets), ],
    #    pointcloud_list=[ new_centroids ], pointcloud_opacity=0.2,
    #    mayavi_wireframe=[False, True,], opacity=[1, 1, 0.9], gradients_at=None, separate=False, gradients_from_iobj=None,
    #    minmax=(RANGE_MIN,RANGE_MAX)  )
    # exit()

    # neighbour_faces_of_vertex
    vertex_neighbours_list = mesh_utils.make_neighbour_faces_of_vertex(facets)
    centroid_gradients = compute_centroid_gradients(new_centroids, iobj)

    new_verts_qem = \
        vertices_apply_qem3(verts, facets, new_centroids, vertex_neighbours_list, centroid_gradients)
        # verts = nv1

    alpha = 0.

    new_verts_qem_alpha = (new_verts_qem * alpha + verts * (1-alpha))

    for i in range(VERTEX_RELAXATION_ITERATIONS_COUNT):
        new_verts_qem_alpha, facets_not_used, centroids = process2_vertex_resampling_relaxation(new_verts_qem_alpha, facets, iobj)
        print("Vertex relaxation applied.");sys.stdout.flush()

    chosen_facet_indices = np.array(total_subdivided_facets, dtype=int)

    centroids2, new_centroids2 = old_centroids[chosen_facet_indices], new_centroids[chosen_facet_indices]

    # move the following code into subdivide_multiple_facets() (?)
    if chosen_facet_indices.size == 0:
        chosen_subset_of_facets = np.zeros((0,), dtype=int)
    else:
        chosen_subset_of_facets = facets[chosen_facet_indices, :]

    highlighted_vertices = np.array([131, 71, 132])  # np.arange(100, 200)
    hv = new_verts_qem[highlighted_vertices, :]

    new_verts_final = comparison_verts_new_verts(verts, new_verts_qem)
    # display_simple_using_mayavi_2([(new_verts_final, facets), (new_verts_qem, facets), ],
    #    pointcloud_list=[hv], pointcloud_opacity=0.2,
    #    mayavi_wireframe=[False, True], opacity=[0.2, 0.5, 0.9], gradients_at=None, separate=False, gradients_from_iobj=None,
    #    minmax=(RANGE_MIN, RANGE_MAX))
    # exit()
    #
    display_simple_using_mayavi_2( [(new_verts_qem_alpha, facets),(new_verts_qem, facets), ],
       pointcloud_list=[ hv ], pointcloud_opacity=0.2,
       mayavi_wireframe=[False,False], opacity=[0.4*0, 1, 0.9], gradients_at=None, separate=False, gradients_from_iobj=None,
       minmax=(RANGE_MIN,RANGE_MAX)  )

#from timeit import default_timer as dtimer


if __name__ == '__main__':

#    timer_t1s = dtimer()

    demo_combination_plus_qem()  # subdivision + projection + qem
