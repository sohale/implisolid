from ipdb import set_trace
import profile_support
from vtk_mc import vtk_mc

import sys
import math

import numpy as np
from basic_functions import check_vector3_vectorized, normalize_vector3_vectorized, normalize_vector4_vectorized

mesh_correction = False


def optimised_used():
    global _optimised_used
    _optimised_used = True

    def side_effect():
        global _optimised_used
        _optimised_used = False
        return True
    assert side_effect()

    return _optimised_used


def mysign_np(v):
    ROOT_TOLERANCE = 0.001
    return np.sign(v) * (np.abs(v) > ROOT_TOLERANCE)


def bisection_vectorized5_(iobj, x1_arr, x2_arr, ROOT_TOLERANCE):
    """ based on bisection_vectorized5. Note that this functin assumes there is no root in x1 and x2."""
    check_vector3_vectorized(x1_arr)
    check_vector3_vectorized(x2_arr)
    assert x1_arr.shape[0] == x2_arr.shape[0]
    v1_arr = iobj.implicitFunction(x1_arr)
    v2_arr = iobj.implicitFunction(x2_arr)

    result_x_arr = np.ones(x1_arr.shape)

    n = x1_arr.shape[0]
    active_indices = np.arange(0, n)  # mid
    active_count = n
    solved_count = 0

    x_mid_arr = np.ones((active_count, 3))
    v_mid_arr = np.zeros((active_count,))

    iteration = 1
    while True:
        # print "iteration", iteration
        # assert np.all(mysign_np(v2_arr[:active_count], ROOT_TOLERANCE) * mysign_np(v1_arr[:active_count], ROOT_TOLERANCE) < 0 - EPS)  # greater or equal
        assert np.all(v1_arr[:active_count] < 0-ROOT_TOLERANCE)
        assert active_indices.shape[0] == x1_arr[:active_count].shape[0]
        assert active_indices.shape[0] == x2_arr[:active_count].shape[0]
        x_mid_arr[:active_count] = (x1_arr[:active_count] + x2_arr[:active_count]) / 2.0
        v_mid_arr[:active_count] = iobj.implicitFunction(x_mid_arr[:active_count, :])
        assert active_indices.shape == (active_count,)
        assert active_indices.ndim == 1
        abs_ = np.abs(v_mid_arr[:active_count])
        indices_boundary = np.nonzero(abs_ <= ROOT_TOLERANCE)[0]  # eq
        indices_outside = np.nonzero(v_mid_arr[:active_count] < -ROOT_TOLERANCE)[0]  # gt
        indices_inside = np.nonzero(v_mid_arr[:active_count] > +ROOT_TOLERANCE)[0]  # -v_mid_arr <  ROOT_TOLERANCE
        indices_eitherside = np.nonzero(abs_ > ROOT_TOLERANCE)[0]
        assert indices_boundary.size + indices_inside.size + indices_outside.size == active_count
        assert indices_eitherside.size + indices_boundary.size == active_count
        which_zeroed = active_indices[indices_boundary]  # new start = mid
        found_count = indices_boundary.shape[0]
        solved_count += found_count
        assert active_count-found_count+solved_count == n
        result_x_arr[which_zeroed] = x_mid_arr[indices_boundary]
        assert np.all(indices_boundary < active_count)
        v2_arr[indices_inside] = v_mid_arr[indices_inside]
        x2_arr[indices_inside] = x_mid_arr[indices_inside]
        v1_arr[indices_outside] = v_mid_arr[indices_outside]
        x1_arr[indices_outside] = x_mid_arr[indices_outside]
        assert np.all(indices_outside < active_count)
        assert np.all(indices_inside < active_count)

        # ------ next round: --------
        assert active_count == active_indices.size
        active_indices = active_indices[indices_eitherside]
        assert active_count - found_count == active_indices.size
        old_active_count = active_count
        active_count = active_count - found_count
        assert active_count == indices_eitherside.size
        # again: does this hold again? assert active_count == active_indices.size
        iteration += 1
        assert np.all(indices_eitherside < old_active_count)
        v1_arr[:active_count] = v1_arr[indices_eitherside]
        v2_arr[:active_count] = v2_arr[indices_eitherside]
        x1_arr[:active_count] = x1_arr[indices_eitherside]
        x2_arr[:active_count] = x2_arr[indices_eitherside]

        assert active_indices.shape == v1_arr[:active_count].shape
        assert active_indices.shape[0] == active_count

        del old_active_count

        assert len(active_indices) == active_count
        if len(active_indices) == 0:
            break

    assert active_indices.size == 0
    optimisation_used = optimised_used()
    if not optimisation_used:
        v_arr = iobj.implicitFunction(result_x_arr)
        assert np.all(np.abs(v_arr) < ROOT_TOLERANCE)
    return result_x_arr



def set_centers_on_surface__ohtake_v3s_002(iobj, centroids, average_edge):
    """ see set_centers_on_surface__ohtake() """
    print "Projecting the centroids: new age"

    # print "s", ; flush()

    THRESHOLD_minimum_gradient_len = 0.000001  # kill gradients smaller than this
    THRESHOLD_zero_interval = 0.0001  # f == TH is NOT zero.
    MAX_ITER = 20

    max_dist = average_edge

    x = centroids

    fc_a = iobj.implicitFunction(x)
    g_a = iobj.implicitGradient(x)
    glen_a = np.linalg.norm(g_a, axis=1)
    glen_a[np.abs(glen_a) < THRESHOLD_minimum_gradient_len] = 1.

    g_normalization_factors = 1. / glen_a[:, np.newaxis]
    g_direction_a = g_a * g_normalization_factors

    signs_c = (fc_a > THRESHOLD_zero_interval)*1. - (fc_a < -THRESHOLD_zero_interval)*1.

    x0 = x

    dx1_c = - g_direction_a * signs_c[:, np.newaxis]

    """
    f_a = fc_a # ???
    taubin = f_a/glen_a
    x1_taubin = x0 - g_direction_a * taubin[:, np.newaxis]
    x1_half = x0 + 0.5*dx1_c * step_size
    x1_half_opposite = x0 - 0.5*dx1_c * step_size
    #Opposite search: Ohtake does not loop the opposite direction if it did not find te point in the forward direction.
    # boundary:
    ...
    candidates = [x1]  # [x1, x1_taubin, x1_half, x1_half_opposite]

    for xa in candidates:
        xa4 = augment4(xa)
        f_a = iobj.implicitFunction(xa4)
        signs = (f_a > THRESHOLD_zero_interval)*1. - (f_a < -THRESHOLD_zero_interval)*1.
        signs = signs*step_size
        #zeros, negatives.
        success = f_a * fc_a <= 0.  # May miss forget about accidental zeros. In that case, use the slightely further points.
        assert success.ndim == 1
    """
    step_size = max_dist / 2. * 2.

    alpha_list = []
    assert step_size > 0.001
    while step_size > 0.001:
        step_size = step_size * 0.5
        max_step = min(MAX_ITER, int(math.floor(max_dist/math.fabs(step_size) + 0.001)))
        assert max_step >= 2  # at least one step
        # if max_step
        # violated only at first time but the first point is already done.
        for i in range(1, max_step+1, 2):  # Step size is two, to avoid aready visited points
            alpha = float(i)*step_size
            # print i, alpha/average_edge
            alpha_list += [alpha/average_edge]
            alpha_list += [-alpha/average_edge]
            # alpha is prepared

    # The algorithm
    n = x0.shape[0]
    best_result_x = np.ones((n, 3))
    active_indices = np.arange(0, n, dtype=int)
    active_count = n
    print n
    del n

    still_nonsuccess_indices = active_indices

    print "points: ", active_count, ".",
    already_success = fc_a*0 > 1.  # all False
    success = already_success.copy()  # falses  #.copy() is necessary
    assert not np.any(already_success)

    # print "left(found)",
    for alpha in alpha_list:
            x1_half = x0 + (max_dist*alpha)*dx1_c
            FAST = True
            if FAST:
                active_indices = still_nonsuccess_indices
                # set_trace()
                # Todo: For those that have changed sign, check if they are closer actually.
                f_a = iobj.implicitFunction(x1_half[active_indices, :])
                signs_a = (f_a > THRESHOLD_zero_interval)*1. + (f_a < -THRESHOLD_zero_interval)*(-1.)
                # success = signs_a * signs_c <= 0.
                success0 = signs_a * signs_c[active_indices] <= 0.
                success[:] = False
                # success[success0] = True
                assert np.all(success == False)
                assert np.all(success[active_indices] == False)
                success[active_indices] = success0
                # print "success:", np.sum(success),
            else:
                # Todo: For those that have changed sign, check if they are closer actually.
                f_a = iobj.implicitFunction(x1_half)
                signs_a = (f_a > THRESHOLD_zero_interval)*1. + (f_a < -THRESHOLD_zero_interval)*(-1.)
                success = signs_a * signs_c <= 0.
                # print "success:", np.sum(success),
            assert success.ndim == 1
            # print "success", np.sum(success)
            # print "already_success", np.sum(already_success)

            new_success_indices = np.nonzero(np.logical_and(success, np.logical_not(already_success)))[0]

            still_nonsuccess_indices = np.nonzero(np.logical_and(np.logical_not(success), np.logical_not(already_success)))[0]
            best_result_x[new_success_indices, :] = x1_half[new_success_indices, :]
            # print "new success>>", new_success_indices.size, "<<  ",
            # print "already>>", np.sum(already_success.size), "<<  ",
            # todo: also try som ein already_success and improve by replacing those that are CLOSER.
            # already_success_but_open_to_improvement = ...
            # best_so_far = ...

            # for next round
            already_success = np.logical_or(success, already_success)  # Union
            # print "left:", still_nonsuccess_indices.shape, ".",
            print ("%d(+%d) " % (still_nonsuccess_indices.size, new_success_indices.size)),
            # print "already_success", np.sum(already_success)
            # active_indices = still_nonsuccess_indices
            if still_nonsuccess_indices.shape[0] == 0:
                break
    # if still_nonsuccess_indices.shape[0] > 0:
    best_result_x[still_nonsuccess_indices, :] = x0[still_nonsuccess_indices, :]  # failed to converge

    TEST = False
    if TEST and not optimised_used():

        f1 = iobj.implicitFunction(x0)
        f2 = iobj.implicitFunction(best_result_x)
        s = f1*f2
        print s[still_nonsuccess_indices]
        s[still_nonsuccess_indices] = -1.
        # print s[s > 0]
        assert np.all(s <= +THRESHOLD_zero_interval)  # can contain almost-zeros. Include the ==equality in zero-ness
        print "OK"

    # print "."; flush()

    # ------------
    # Prepare for bisection: By removing zeros and moving negatives to x1 by swapping.

    # collect the zero ones

    f1 = fc_a
    assert np.all(f1 == iobj.implicitFunction(x0), axis=None)
    f2 = iobj.implicitFunction(best_result_x)
    zeros2 = np.abs(f2) <= THRESHOLD_zero_interval
    # Copy the zeros onto results
    zeros1 = np.abs(f1) <= THRESHOLD_zero_interval
    best_result_x[zeros1, :3] = x0[zeros1, :3]
    zeros12 = np.logical_or(zeros1, zeros2)  # output
    assert np.all(np.abs(iobj.implicitFunction(best_result_x[zeros12, :])) <= THRESHOLD_zero_interval)

    ROOT_TOLERANCE = THRESHOLD_zero_interval  # because of the assert

    relevants_boolean = np.logical_and(already_success, np.logical_not(zeros12))

    assert np.all(np.abs(f2[relevants_boolean]) > +THRESHOLD_zero_interval)
    assert np.all(np.abs(f2[zeros12]) <= +THRESHOLD_zero_interval)
#    import ipdb; ipdb.set_trace()
#    assert np.all(mysign_np(f2[relevants_boolean], THRESHOLD_zero_interval) * mysign_np(f1[relevants_boolean], THRESHOLD_zero_interval) < 0)

    x0_v4 = centroids[relevants_boolean, :]
    x2_v4 = best_result_x[relevants_boolean, :]
    f1 = iobj.implicitFunction(x0_v4)
    f2 = iobj.implicitFunction(x2_v4)
    s = f1*f2
    assert np.all(s <= +THRESHOLD_zero_interval)

    # Swap negatives and positives
    swap = f2 < -THRESHOLD_zero_interval
    temp = x2_v4[swap, :]
    x2_v4[swap, :] = x0_v4[swap, :]
    x0_v4[swap, :] = temp
    temp_f = f2[swap]
    f2[swap] = f1[swap]
    f1[swap] = temp_f

    bsresults = bisection_vectorized5_(iobj, x0_v4, x2_v4, ROOT_TOLERANCE)
    assert bsresults.shape[0] == np.sum(relevants_boolean)
    centroids[relevants_boolean, :] = bsresults[:, :]  # x4
    return


def compute_triangle_areas(verts, faces, return_normals=False):
    """ facet_normals: can contain NaN if the area is zero"""

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
        return facet_areas
    else:
        print facet_areas.shape
        assert facet_areas[:, np.newaxis].shape == (nfaces, 1)
        facet_normals = a / np.tile(facet_areas[:, np.newaxis], (1, 3)) / 2.0
        return facet_areas, facet_normals


def build_faces_of_faces(facets):
    """ builds lookup tables. The result if an array of nfaces x 3,
    containing the face index of neighbours of each face.
    Since each face has exactly three neighbours, the size of the result is n x 3."""
    from mesh_utils import make_edge_lookup_sparse
    (edges_of_faces, faces_of_edges, vertpairs_of_edges) = \
        make_edge_lookup_sparse(facets)

    nfaces = facets.shape[0]
    assert edges_of_faces.shape == (nfaces, 3)
    f12 = faces_of_edges[edges_of_faces, :]
    print f12.shape
    assert f12.shape == (nfaces, 3, 2)
    faceindex_eye = np.tile(np.arange(nfaces)[:, np.newaxis, np.newaxis], (1, 3, 2))
    assert faceindex_eye.shape == (nfaces, 3, 2)

    f12 = f12 + 1
    f12[f12 == faceindex_eye+1] = 0

    assert np.allclose(np.prod(f12, axis=2), 0)

    f_uniq = np.sum(f12, axis=2)
    assert np.sum(f_uniq.ravel() == 0) == 0

    return f_uniq - 1


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

        if mimj > 1.0:
            mimj = 1.0
        if mimj < -1.0:
            mimj = -1.0

        pipj = np.linalg.norm(pi - pj)
        if pipj == 0:
            return 0

        assert pipj == np.linalg.norm(pi - pj, ord=2)

        kij = np.arccos(mimj) / pipj
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
    print "w_list ", w_list

    print("===============")
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

        color_list = [(1, 0, 0), (0, 0, 0), (1, 1, 0), (0, 0, 1), (0, 1, 0)]
        i = 0
        for c in pointcloud_list:
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


def get_A_b(vertex_id, nlist_numpy, centroids, centroids_gradients, qem_origin):

    nai = nlist_numpy

    center_array = centroids[nai, :]
    # note some centers may not be projected successfully in the previous step
    not_projected_successfully = np.isnan(center_array[:].ravel())
    if np.any(not_projected_successfully):
        pass

    normals = centroids_gradients[nai, :]  # why do we have repeats??
    # note : not normalised. But it works.new_weight[i])

    norms = np.linalg.norm(normals, ord=2, axis=1)
    # can be either 0, 1 or Nan
    if np.any(norms < 0.000001):  # can be exactly 0.0
        print("Error: bad normal", normals)

    # TH_N = 0.0000001  # 0.000001 = I millions of millimeter = 1 nanometer
    # can be 0,0,0, inf, nonsharp, degenerate, ...
    # degenerate_normals = np.logical_or(np.isnan( np.sum(normals, axis=1)), norms < TH_N )
    assert not np.any(np.isnan(normals))
    assert not np.any(np.isinf(normals))

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


def vertices_apply_qem3(verts, facets, centroids, vertex_neighbours_list, centroids_gradients):
    assert centroids is not None
    assert vertex_neighbours_list is not None
    assert centroids_gradients is not None

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

        A, b = get_A_b(vi, nai, centroids, centroids_gradients, qem_origin)

        u, s, v = np.linalg.svd(A)
        assert np.allclose(A, np.dot(u, np.dot(np.diag(s), v)))
        assert s[0] == np.max(s)

        # experimental value
        tau = 680.  # 10. ** 2.83

        s[s / s[0] < 1.0/tau] = 0
        rank = np.sum(s / s[0] > 1.0/tau)

        assert np.all(s[:rank]/s[0] >= 1.0/tau)

        x = verts[vi, 0:3, np.newaxis]
        assert x.shape == (3, 1)

        y = np.dot(v, x).copy()
        utb = np.dot(-np.transpose(u), b)

        for i in range(rank):
            assert np.dot(-np.transpose(u), b).shape == (3, 1)
            assert s[i]/s[0] >= 1.0/tau
            y[i] = utb[i] / s[i]
        new_x = np.dot(np.transpose(v), y)

        if not s[0] > 0.000001:
            print("Warning! sigma_1 == 0")
            print(s)
            print("A", A)

            result_verts_ranks[vi] = 0

        assert x.shape == (3, 1)

        new_verts[vi, 0:3] = new_x[:, 0]

        if not np.all(np.abs(utb.ravel()[rank:]) < 0.0001):
            pass
        result_verts_ranks[vi] = rank

    print("max rank = ", np.max(result_verts_ranks))
    print("min rank = ", np.min(result_verts_ranks))
    if not np.min(result_verts_ranks) >= 1:
        print("Warning: assertion: np.min(result_verts_ranks) >= 1 failed.")

    if False:
        assert np.min(result_verts_ranks) >= 1
    return new_verts


import mesh_utils


def compute_centroids(verts, facets):
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
    faces_of_faces = build_faces_of_faces(facets)

    new_verts = vertex_resampling(verts, faceslist_neighbours_of_vertex, faces_of_faces, centroids, centroid_normals_normalized, c=2.0)

    return new_verts, facets, centroids


def subdivide_multiple_facets(verts_old, facets_old, tobe_subdivided_face_indices):


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

    facet_areas, facet_normals = compute_triangle_areas(verts, facets, return_normals=True)

    nf = facets.shape[0]
    assert facet_areas.shape == (nf,)
    assert facet_normals.shape == (nf, 3)

    assert np.all(np.logical_not(np.isnan(facet_areas[np.logical_not(np.isnan(np.linalg.norm(facet_normals, axis=1)))])))

    degenerate_faces = np.isnan(facet_areas)

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

    e_array = np.ndarray(nf)

    for fi in range(nf):
        if degenerate_faces[fi]:
            e_array[fi] = 0
            continue

        n = facet_normals[fi, :]  # n: (3,)

        triangle = verts[facets[fi, :], :]  # numverts x 3
        assert triangle.shape == (3, 3)

        assert triangle.shape == (3, 3)
        VVV = triangle

        m0123 = np.dot(centroidmaker_matrix, np.dot(subdiv_vert_matrix, VVV))
        assert m0123.shape == (4, 3)
        subdiv_centroids = m0123

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

        e_array[fi] = e

        if fi % 100 == 0:
            print fi, "\r", ;import sys; sys.stdout.flush()

    num_subivision = len(e_array[e_array > curvature_epsilon])
    need_subidivision = np.ndarray(num_subivision, dtype=int)

    k = 0
    for i in range(nf):
        if e_array[i] > curvature_epsilon:
            need_subidivision[k] = i
            k += 1

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


def simple_histogram(c, title=None, special_values=[]):
    import matplotlib.pyplot as plt
    # special_values

    plt.hist(c, 20*10)
    if title is not None:
        plt.title(title)
    plt.show()

@profile
def demo_combination_plus_qem():
    """ Now with QEM """
    curvature_epsilon = 1. / 1000.  # a>eps  1/a > 1/eps = 2000
    # curvature_epsilon = 1. / 10000.
    VERTEX_RELAXATION_ITERATIONS_COUNT = 0
    SUBDIVISION_ITERATIONS_COUNT = 0  # 2  # 5+4

    from example_objects import make_example_vectorized
    object_name = "cube_with_cylinders"  # "sphere_example" #or "rcube_vec" work well #"ell_example1"#"cube_with_cylinders"#"ell_example1"  " #"rdice_vec" #"cube_example"
    iobj = make_example_vectorized(object_name)

    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-3, +5, 0.2)
    if object_name == "cube_with_cylinders" or object_name == "twist_object" or object_name == "french_fries" or object_name == "rdice_vec" or object_name == "rods" or object_name == "bowl_15_holes":
        VERTEX_RELAXATION_ITERATIONS_COUNT = 1

    if object_name == "cyl4":
        (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-32 / 2, +32 / 2, 1.92 / 4.0)

    elif object_name == "french_fries" or object_name == "rods":
        (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-3, +5, 0.11)  # 0.05

    elif object_name == "bowl_15_holes":
        (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-3, +5, 0.15)

    elif object_name == "cyl3":
        (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-32 / 2, +32 / 2, 1.92 / 4.0)

    elif object_name == "cyl1":
        (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-16, +32, 1.92 * 0.2 * 10 / 2.0)

    from stl_tests import make_mc_values_grid
    gridvals = make_mc_values_grid(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE, old=False)
    verts, facets = vtk_mc(gridvals, (RANGE_MIN, RANGE_MAX, STEPSIZE))
    print("MC calculated.");sys.stdout.flush()


    # display_simple_using_mayavi_2([(verts, facets), ],
    #    pointcloud_list=[],
    #    mayavi_wireframe=[False], opacity=[1], gradients_at=None, separate=False, gradients_from_iobj=None,
    #    minmax=(RANGE_MIN, RANGE_MAX))
    # exit()

    for i in range(VERTEX_RELAXATION_ITERATIONS_COUNT):
        verts, facets_not_used, centroids = process2_vertex_resampling_relaxation(verts, facets, iobj)
        assert not np.any(np.isnan(verts.ravel()))  # fails
        print("Vertex relaxation applied.");sys.stdout.flush()

    average_edge = compute_average_edge_length(verts, facets)

    old_centroids = np.mean(verts[facets[:], :], axis=1)

    new_centroids = old_centroids.copy()

    set_centers_on_surface__ohtake_v3s_002(iobj, new_centroids, average_edge)

    # neighbour_faces_of_vertex
    vertex_neighbours_list = mesh_utils.make_neighbour_faces_of_vertex(facets)
    centroid_gradients = compute_centroid_gradients(new_centroids, iobj)

    new_verts_qem = \
        vertices_apply_qem3(verts, facets, new_centroids, vertex_neighbours_list, centroid_gradients)

    verts_before_qem = verts

    highlighted_vertices = np.array([131, 71, 132])  # np.arange(100, 200)
    hv = new_verts_qem[highlighted_vertices, :]

    total_subdivided_facets = []
    for i in range(SUBDIVISION_ITERATIONS_COUNT):
        e_array, which_facets = compute_facets_subdivision_curvatures(new_verts_qem, facets, iobj, curvature_epsilon)

        # which_facets = np.arange(facets.shape[0])[e_array > curvature_epsilon]
        print "Curvature epsilon:", curvature_epsilon, "which facets need to be subdivided", which_facets.shape

        verts4_subdivided, facets3_subdivided = subdivide_multiple_facets(new_verts_qem, facets, which_facets)
        global trace_subdivided_facets  # third implicit output

        verts, facets = verts4_subdivided, facets3_subdivided
        print("Subdivision applied.");sys.stdout.flush()
        # total_subdivided_facets += trace_subdivided_facets  # old face indices remain valid

        # new_verts_qem_alpha = (new_verts_qem * alpha + verts * (1-alpha))

        highlighted_vertices = np.array([131, 71, 132])  # np.arange(100, 200)
        hv = verts[highlighted_vertices, :]

    chosen_facet_indices = np.array(total_subdivided_facets, dtype=int)

    # move the following code into subdivide_multiple_facets() (?)
    if chosen_facet_indices.size == 0:
        chosen_subset_of_facets = np.zeros((0,), dtype=int)
    else:
        chosen_subset_of_facets = facets[chosen_facet_indices, :]

    # display_simple_using_mayavi_2([(new_verts_final, facets), (new_verts_qem, facets), ],
    #    pointcloud_list=[hv], pointcloud_opacity=0.2,
    #    mayavi_wireframe=[False, True], opacity=[0.2, 0.5, 0.9], gradients_at=None, separate=False, gradients_from_iobj=None,
    #    minmax=(RANGE_MIN, RANGE_MAX))
    # exit()
    #
    display_simple_using_mayavi_2([(verts_before_qem, facets), (new_verts_qem, facets), ],
       pointcloud_list=[hv], pointcloud_opacity=0.2,
       mayavi_wireframe=[False, False], opacity=[0.4*0, 1, 0.9], gradients_at=None, separate=False, gradients_from_iobj=None,
       minmax=(RANGE_MIN, RANGE_MAX))

# from timeit import default_timer as dtimer


if __name__ == '__main__':

#    timer_t1s = dtimer()

    demo_combination_plus_qem()  # subdivision + projection + qem
