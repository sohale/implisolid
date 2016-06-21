from ipdb import set_trace
import profile_support
from vtk_mc import vtk_mc

import sys
import math

import numpy as np
from basic_functions import check_vector3_vectorized, normalize_vector3_vectorized, normalize_vector4_vectorized

mesh_correction = False

B = 1000000L


def optimised_used():
    global _optimised_used
    _optimised_used = True

    def side_effect():
        global _optimised_used
        _optimised_used = False
        return True
    assert side_effect()

    return _optimised_used


def mysign_np(v, ROOT_TOLERANCE=0.001):
    # ROOT_TOLERANCE = 0.001
    return np.sign(v) * (np.abs(v) > ROOT_TOLERANCE)


def isomorphic(a, b):
    if not np.ndim(a) == np.ndim(b):
        return False
    if not a.shape == b.shape:
        return False
    # if not np.all(a == b):
    #    return False
    assert a.dtype.type != np.bool
    assert b.dtype.type != np.bool
    return True


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


def compute_triangle_areas(vertex, faces, return_normals=False):
    """ facet_normals: can contain NaN if the area is zero"""

    DEGENERACY_THRESHOLD = 0.00001
    nfaces = faces.shape[0]
    expand = vertex[faces, :]

    assert expand.shape == (nfaces, 3, 3)
    assert expand[:, 2, :].shape == (nfaces, 3)
    a = np.cross(
        expand[:, 1, :] - expand[:, 0, :],
        expand[:, 2, :] - expand[:, 0, :],
        axis=1)
    faces_areas = np.linalg.norm(a, axis=1, ord=2) / 2.0
    degenerates_count = len(faces_areas[faces_areas < DEGENERACY_THRESHOLD])
    faces_areas[faces_areas < DEGENERACY_THRESHOLD] = np.nan  # -1
    if degenerates_count > 0:
        print("degenerate triangles", degenerates_count)
    if not return_normals:
        return faces_areas
    else:
        print faces_areas.shape
        assert faces_areas[:, np.newaxis].shape == (nfaces, 1)
        facet_normals = a / np.tile(faces_areas[:, np.newaxis], (1, 3)) / 2.0
        return faces_areas, facet_normals


def build_faces_of_faces(faces):
    """ builds lookup tables. The result if an array of nfaces x 3,
    containing the face index of neighbours of each face.
    Since each face has exactly three neighbours, the size of the result is n x 3."""
    from mesh_utils import make_edge_lookup_sparse
    (edges_of_faces, faces_of_edges) = \
        make_edge_lookup_sparse(faces)

    nfaces = faces.shape[0]
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


def vertex_resampling(vertex, faceslist_neighbours_of_vertex, faces_of_faces, centroids, centroid_normals, c=2.0):
    """ faceslist_neighbours_of_vertex: *** """

    def kij(i, j):
        """ Returns (1/r * Theta), a measure of curvature.
        Theta is the angle between two normals at centroids (dual vertex) i, j.
        The 1/r is the inverse of the distance between the pair.

        Notes:
        Normals should be already normalised (centroid normals).
        centroids should be already [projected] on the implicit surface (dual mesh is optimised).
        """

        # i,j are centroids
        assert i != j
        pi, pj = (centroids[i, :], centroids[j, :])

        # based on gradients. normalised.
        mi, mj = (centroid_normals[i, :], centroid_normals[j, :])
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

    def wi(i_faces, ja_faces, c):
        """
        Returns the weight of a facet i_facet.
        Adds kij of all centroids of the neighbour facets.
        ja_facets = list of centroid indices (face index).
        i_facet is a face index. """

        assert i_faces not in ja_faces
        assert len(ja_faces) == 3

        ki = 0
        for j_faces in ja_faces:
            ki += kij(i_faces, j_faces)
        wi = 1.0 + c*ki

        return wi
    #
    c_ = c  # 2.0  # constant


    umbrella_faces = faceslist_neighbours_of_vertex[1]  # A list of facets: The indices of faces that vertex vertex_index belongs to.
    print("umbrella_faces: ", umbrella_faces)
    # wa = np.zeros()
    w_list = []
    for i_faces in umbrella_faces:

        print("i_facet", i_faces)
        three_faces = faces_of_faces[i_faces, :]
        print(i_faces, three_faces)
        w = wi(i_faces, three_faces, c_)  # three_facets should be neighbours of the facet i_facet
        # The weight (based on curvature) of neighbour P_i (facet i.e. centroid),
        print("w_i, i=", i_faces, w)
        w_list.append(w)
    print "w_list ", w_list

    print("===============")
    # w seems tobe calculated fine. next: store w_i and cache them for adaptive resampling, for which we need to normalise it across the neighbours.
    nfaces = centroids.shape[0]
    wi_total_array = np.zeros((nfaces,))
    for i_faces in range(nfaces):
        three_faces = faces_of_faces[i_faces, :]
        w = wi(i_faces, three_faces, c_)
        wi_total_array[i_faces] = w
    print wi_total_array
    # The weights are prepared. Now let's resample vertex

    umbrella_faces = np.array(faceslist_neighbours_of_vertex[1])  # empty
    print "umbrella_faces", umbrella_faces.shape, "****"
    assert np.allclose(wi_total_array[umbrella_faces] - np.array(w_list), 0)

    def lift_vertex(vertex, centroids):
        new_vertex = vertex.copy()
        # assign these to a sparse matrix? and  do:  M = M/normalise(M); vertex = M * vertex
        for i in range(vertex.shape[0]):
            umbrella_faces = np.array(faceslist_neighbours_of_vertex[i])
            w = wi_total_array[umbrella_faces]
            # w = w * 0 + 1
            w = w / np.sum(w)
            # print w / np.sum(w), w.shape
            new_vertex[i, :] = \
                np.dot(w, centroids[umbrella_faces, :])  # (n) * (n x 3)
        return new_vertex

    return lift_vertex(vertex, centroids)


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
        vertex, faces = vf

        if vertex is None:
            continue
        if vertex.size == 0:
            print("Warning: empty vertex")
            continue
        if faces.size == 0:
            print("Warning: no faces")
            continue

        assert vertex.ndim == 2
        assert faces.ndim == 2
        assert vertex.shape == (vertex.shape[0], 3), str(vertex.shape)
        assert faces.shape == (faces.shape[0], 3), str(faces.shape)
        if type(mayavi_wireframe) is list:
            wire_frame1 = mayavi_wireframe[fi]
            assert len(mayavi_wireframe) == len(vf_list)
        else:
            wire_frame1 = mayavi_wireframe
        mlab.triangular_mesh([vert[0] for vert in vertex],
                         [vert[1] for vert in vertex],
                         [vert[2] for vert in vertex],
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
        vertex1, faces1 = vf_list[0]
        avg_edge_len = compute_average_edge_length(vertex, faces)
        visualise_gradients(mlab, gradients_at, gradients_from_iobj, avg_edge_len / 20.)
    if gradients_from_iobj is not None:
        add_random_interior_points(mlab, gradients_from_iobj, avg_edge_len)

    mlab.show()
    return


def compute_average_edge_length(vertex, faces):
    nfaces = faces.shape[0]
    expand = vertex[faces, :]
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


def get_A_b(nlist_numpy, centroids, centroids_gradients, qem_origin):

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


def set_centers_on_surface__ohtake_v3s(iobj, centroids, average_edge, debug_vf=None, mesh_normals=None):
    """ see set_centers_on_surface__ohtake() """
    print "Projecting the centroids: new age"


    THRESHOLD_minimum_gradient_len = 0.000001  # kill gradients smaller than this
    THRESHOLD_zero_interval = 0.0001  # f == TH is NOT zero.
    MAX_ITER = 20
    USE_MESH_NORMALS = True

    max_dist = average_edge

    x = centroids

    fc_a = iobj.implicitFunction(x)
    g_a = iobj.implicitGradient(x)[:, :3]
    glen_a = np.linalg.norm(g_a, axis=1)
    glen_a[np.abs(glen_a) < THRESHOLD_minimum_gradient_len] = 1.

    g_normalization_factors = 1. / glen_a[:, np.newaxis]
    g_direction_a = g_a * g_normalization_factors
    # The directions toward the center

    # step_size = max_dist / 2.  # average_edge / 2.

    signs_c = (fc_a > THRESHOLD_zero_interval)*1. - (fc_a < -THRESHOLD_zero_interval)*1.
    # todo: stop moving when hit zero
    x0_v3 = x[:, :3]
    # Move opposite the direction toward center if the value is positive.

    dx0_c_grad = - g_direction_a * signs_c[:, np.newaxis]

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

    print "Alphas:"
    for i in range(len(alpha_list)):
        print "[%d]%g" % (i, alpha_list[i]),
    print
    print

    # The algorithm
    n = x0_v3.shape[0]
    best_result_x = np.ones((n, 3))
    active_indices = np.arange(0, n, dtype=int)
    active_count = n
    del n

    still_nonsuccess_indices = active_indices

    print "points: ", active_count, ".",
    already_success = fc_a*0 > 1.  # all False
    success = already_success.copy()  # falses  #.copy() is necessary
    assert not np.any(already_success)

    USE_MESH_NORMALS = False
    if USE_MESH_NORMALS:
        assert mesh_normals is not None
        dx0c_mesh_normals = mesh_normals
        assert np.allclose(np.linalg.norm(dx0c_mesh_normals, axis=1), 1.)

    TEST = True  # and not optimised_used():

    # print "left(found)",
    for it in [0, 1] if USE_MESH_NORMALS else [0]:

        if it == 0:
            dxc = dx0_c_grad
            alpha_list1 = alpha_list
        else:
            dxc = dx0c_mesh_normals  # * 0.5
            print
            print "now mesh normals"
            alpha_list1 = alpha_list[:10]

        counter = -1
        for alpha in alpha_list1:
            counter += 1
            x1_half = x0_v3 + (max_dist*alpha)*dxc
            FAST = True
            if FAST:
                active_indices = still_nonsuccess_indices

                # Todo: For those that have changed sign, check if they are closer actually.
                xa4 = x1_half[active_indices, :]
                f_a = iobj.implicitFunction(xa4)
                signs_a = (f_a > THRESHOLD_zero_interval)*1. + (f_a < -THRESHOLD_zero_interval)*(-1.)
                # success = signs_a * signs_c <= 0.
                success0 = signs_a * signs_c[active_indices] <= 0.
                success[:] = False
                assert np.all(success == False)
                assert np.all(success[active_indices] == False)
                success[active_indices] = success0
            else:
                # Todo: For those that have changed sign, check if they are closer actually.
                xa4 = x1_half
                f_a = iobj.implicitFunction(xa4)
                signs_a = (f_a > THRESHOLD_zero_interval)*1. + (f_a < -THRESHOLD_zero_interval)*(-1.)
                success = signs_a * signs_c <= 0.

            assert success.ndim == 1

            new_success_indices = np.nonzero(np.logical_and(success, np.logical_not(already_success)))[0]

            still_nonsuccess_indices = np.nonzero(np.logical_and(np.logical_not(success), np.logical_not(already_success)))[0]
            best_result_x[new_success_indices, :] = x1_half[new_success_indices, :]

            # for next round
            already_success = np.logical_or(success, already_success)  # Union

            print ("[%d](+%d)%d " % (counter, new_success_indices.size, still_nonsuccess_indices.size,)),

            if TEST:
                x2 = best_result_x[already_success, :]
                f2 = iobj.implicitFunction(x2)
                f2[np.abs(f2) < THRESHOLD_zero_interval] = 0.
                x1 = centroids[already_success, :]
                assert x1.shape[1] == 3
                f1 = iobj.implicitFunction(x1)
                f1[np.abs(f1) < THRESHOLD_zero_interval] = 0.

                bad_accepted = (f1*f2)[f1 * f2 > 0]
                if bad_accepted.size > 0:
                        print "bad accepted f1*f2:", (f1*f2)[f1 * f2 > 0]
                # assert np.all((f1*f2)[f1 * f2 > 0])
                assert np.all(f1 * f2 <= 0)
                del f1, f2, x1, x2

            if still_nonsuccess_indices.shape[0] == 0:
                break
    # if still_nonsuccess_indices.shape[0] > 0:
    best_result_x[still_nonsuccess_indices, :] = x0_v3[still_nonsuccess_indices, :]  # failed to converge

    if TEST:
        xa1 = x0_v3
        f1_test = iobj.implicitFunction(xa1)
        xa2 = best_result_x
        f2_test = iobj.implicitFunction(xa2)
        s = f1_test*f2_test

        s[still_nonsuccess_indices] = -1.

        assert np.all(s <= +THRESHOLD_zero_interval)  # can contain almost-zeros. Include the ==equality in zero-ness
        del s
        del f1_test
        del f2_test
        print "OK"

    # ------------
    # Prepare for bisection: By removing zeros and moving negatives to x1 by swapping.

    # collect the zero ones
    xa1 = x0_v3
    f1 = fc_a  # iobj.implicitFunction(xa1)
    assert np.all(f1 == iobj.implicitFunction(xa1), axis=None)
    xa2 = best_result_x
    f2 = iobj.implicitFunction(xa2)
    zeros2_bool = np.abs(f2) <= THRESHOLD_zero_interval
    # Copy the zeros onto results
    zeros1_bool = np.abs(f1) <= THRESHOLD_zero_interval
    best_result_x[zeros1_bool, :3] = x0_v3[zeros1_bool, :3]
    zeros1or2 = np.logical_or(zeros1_bool, zeros2_bool)  # output
    del zeros1_bool, zeros2_bool
    assert np.all(np.abs(iobj.implicitFunction(best_result_x[zeros1or2, :])) <= THRESHOLD_zero_interval)

    ROOT_TOLERANCE = THRESHOLD_zero_interval  # because of the assert

    relevants_bool = np.logical_and(already_success, np.logical_not(zeros1or2))

    assert np.all(np.abs(f2[relevants_bool]) > +THRESHOLD_zero_interval)
    assert np.all(np.abs(f2[zeros1or2]) <= +THRESHOLD_zero_interval)

    assert np.all(mysign_np(f2[relevants_bool], THRESHOLD_zero_interval) * mysign_np(f1[relevants_bool], THRESHOLD_zero_interval) < 0)
    del f1

    x2_relevant_v4 = best_result_x[relevants_bool, :]
    x1_relevant_v4 = centroids[relevants_bool, :]
    f1_relevants = iobj.implicitFunction(x1_relevant_v4)  # for assert only
    f2_relevants = iobj.implicitFunction(x2_relevant_v4)

    assert np.all(f1_relevants*f2_relevants <= +THRESHOLD_zero_interval)
    del f1_relevants

    # Swap negatives and positives
    swap_bool = f2_relevants < -THRESHOLD_zero_interval
    # performance: boolean or indices?
    temp = x2_relevant_v4[swap_bool, :]
    x2_relevant_v4[swap_bool, :] = x1_relevant_v4[swap_bool, :]
    x1_relevant_v4[swap_bool, :] = temp
    del temp

    del f2_relevants

    x_bisect = bisection_vectorized5_(iobj, x1_relevant_v4, x2_relevant_v4, ROOT_TOLERANCE)
    assert np.all(np.abs(iobj.implicitFunction(x_bisect)) < THRESHOLD_zero_interval)
    assert x_bisect.shape[0] == np.sum(relevants_bool)
    centroids[relevants_bool, :] = x_bisect[:, :]  # x4

    print "centroids:", centroids.shape, "best_result_x:", best_result_x.shape, "relevants_bool:", np.sum(relevants_bool), "+", np.sum(np.logical_not(relevants_bool)), "zeros1or2=", np.sum(zeros1or2)
    centroids[zeros1or2, :3] = best_result_x[zeros1or2, :]  # x3

    # relevants_bool = np.logical_and(already_success, np.logical_not(zeros1or2))

    """
    if debug_vf is not None:
        _vs, _fs = debug_vf
        from visual5 import *
        c3 = centroids[zeros1or2, :3]
        #(np.zeros((0, 4)), np.zeros((0, 3), dtype=int))
        display_simple_using_mayavi_2([(_vs, _fs)],
                   pointcloud_list=[centroids[zeros1or2, :]], pointsizes=[0.02], #pointcloud_list=[point_collector.get_as_array()], pointsizes=[0.01],
                   mayavi_wireframe=[False,], opacity=[0.4,],
                   gradients_at=c3,
                   #separate_panels=False,
                   gradients_from_iobj=iobj,
                   #minmax=(RANGE_MIN, RANGE_MAX),
                   #add_noise=[0, 0], noise_added_before_broadcast=True,
                   labels=(centroids, zeros1or2), grad_arrow_len=0.2/2.)

        set_trace()
    """
    return


def vertex_apply_qem3(vertex, faces, centroids, vertex_neighbours_list, centroids_gradients):
    assert centroids is not None
    assert vertex_neighbours_list is not None
    assert centroids_gradients is not None

    nvert = vertex.shape[0]
    assert nvert == len(vertex_neighbours_list)

    result_vertex_ranks = np.zeros((nvert,), dtype=int)
    assert vertex.shape == (nvert, 3)
    new_vertex = np.zeros((nvert, 3))

    for vertex_id in range(nvert):

        vi = vertex_id
        nlist = vertex_neighbours_list[vertex_id]
        nai = np.array(nlist)
        qem_origin = vertex[vertex_id, :].reshape(3, 1)*0

        A, b = get_A_b(nai, centroids, centroids_gradients, qem_origin)
        u, s, v = np.linalg.svd(A)
        assert np.allclose(A, np.dot(u, np.dot(np.diag(s), v)))
        assert s[0] == np.max(s)

        # experimental value
        tau = 680.  # 10. ** 2.83
        s[s / s[0] < 1.0/tau] = 0
        rank = np.sum(s / s[0] > 1.0/tau)

        assert np.all(s[:rank]/s[0] >= 1.0/tau)

        x = vertex[vi, 0:3, np.newaxis]
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

            result_vertex_ranks[vi] = 0

        assert x.shape == (3, 1)

        new_vertex[vi, 0:3] = new_x[:, 0]

        if not np.all(np.abs(utb.ravel()[rank:]) < 0.0001):
            pass
        result_vertex_ranks[vi] = rank

    print("max rank = ", np.max(result_vertex_ranks))
    print("min rank = ", np.min(result_vertex_ranks))
    if not np.min(result_vertex_ranks) >= 1:
        print("Warning: assertion: np.min(result_vertex_ranks) >= 1 failed.")

    if False:
        assert np.min(result_vertex_ranks) >= 1
    return new_vertex


import mesh_utils


def compute_centroids(vertex, faces):
    expand = vertex[faces, :]
    nfacets = faces.shape[0]
    assert expand.shape == (nfacets, 3, 3)
    assert np.allclose(vertex[faces[:], :], expand)
    centroids = np.mean(vertex[faces[:], :], axis=1)  # again
    return centroids


def process2_vertex_resampling_relaxation(vertex, faces, iobj):
    assert not np.any(np.isnan(vertex))
    centroids = compute_centroids(vertex, faces)
    centroid_normals_normalized = compute_centroid_gradients(centroids, iobj, normalise=True)

    from mesh_utils import make_neighbour_faces_of_vertex
    faceslist_neighbours_of_vertex = make_neighbour_faces_of_vertex(faces)
    faces_of_faces = build_faces_of_faces(faces)

    new_vertex = vertex_resampling(vertex, faceslist_neighbours_of_vertex, faces_of_faces, centroids, centroid_normals_normalized, c=2.0)

    return new_vertex, faces, centroids


def subdivide_multiple_faces(vertex_old, faces_old, tobe_subdivided_face_indices):

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

    provisional_new_vertex_count = 3*len(tobe_subdivided_face_indices)
    provisional_new_faces_count = 3*len(tobe_subdivided_face_indices)
    nvertex_old = vertex_old.shape[0]
    nfaces_old = faces_old.shape[0]
    new_vertex = np.zeros((nvertex_old+provisional_new_vertex_count, 3), dtype=float)
    new_faces = np.zeros((nfaces_old+provisional_new_faces_count, 3), dtype=int)
    new_vertex[:nvertex_old, :] = vertex_old
    new_faces[:nfaces_old, :] = faces_old

    # on number of added vertex:
    # problem: there may be repeated (Redundant) vertex. (as well as T-junctions)
    # also later check for faces with repeated edges. (which can be another cause of null normals)

    new_vertex_counter = nvertex_old
    new_faces_counter = nfaces_old
    for subdiv_i in range(len(tobe_subdivided_face_indices)):
        fi = tobe_subdivided_face_indices[subdiv_i]
        oldtriangle = vertex_old[faces_old[fi, :], :]  # numvertex x 3
        assert oldtriangle.shape == (3, 3)
        VVV = oldtriangle  # (nv=3) x 3

        # new verices
        m0123 = np.dot(np.dot(centroidmaker_matrix, subdiv_vert_matrix), VVV)
        assert m0123.shape == (4, 3)

        vxyz_0123 = np.dot(subdiv_vert_matrix, VVV)  # not efficient
        assert vxyz_0123.shape == (6, 3)

        v012 = faces_old[fi, :].tolist()  # range(0, 3)  #
        v345 = range(new_vertex_counter, new_vertex_counter+3)

        v345_xyz = vxyz_0123[3:6, :]  # only pick the new ones

        assert len(v345) == 3
        new_vertex[(new_vertex_counter):(new_vertex_counter+3), :] = v345_xyz

        new_vertex_counter += 3

        # facet's vertex indices
        v012345 = np.array(v012 + v345, dtype=int)
        mini_faces_l = [[0, 3, 5], [3, 1, 4], [5, 4, 2], [3, 4, 5]]  # 0,3,1,4,2,5
        mini_faces = v012345[np.array(mini_faces_l, dtype=int)]

        new_faces[fi, :] = mini_faces[0, :]
        new_faces[new_faces_counter:(new_faces_counter+3), :] = mini_faces[1:(1+3), :]
        assert mini_faces.shape[0] == (1+3)
        trace_subdivided_facets += range(new_faces_counter, (new_faces_counter+3)) + [fi]  # include the face which reuses the old face's index
        # trace_subdivided_facets will contain indices of faces

        new_faces_counter += 3

        if fi % 100 == 0:
            print fi, "\r",
            import sys
            sys.stdout.flush()

    print new_vertex.shape[0], new_vertex_counter

    assert new_vertex.shape[0] == new_vertex_counter
    assert new_faces.shape[0] == new_faces_counter
    print "v", provisional_new_vertex_count+nvertex_old, new_vertex_counter
    print "f", provisional_new_faces_count+nfaces_old, new_faces_counter
    assert provisional_new_vertex_count+nvertex_old == new_vertex_counter
    assert provisional_new_faces_count+nfaces_old == new_faces_counter
    assert len(trace_subdivided_facets) == 0 or np.max(np.array(trace_subdivided_facets)) < new_faces_counter
    return new_vertex, new_faces


def subdivide_multiple_faces_new(vertex_old, faces_old, tobe_subdivided_face_indices, midpoint_map):
    """
    midpoint_map is modified (is input and output).
    midpoint_map is a dictionary that given an edge's unique_int_id, gives you the vertex in the midpoint. It may contain midpoints that are not used anymore.
    Use compute_facets_subdivision_curvatures() to calculate tobe_subdivided_face_indices.
    Does not remove vertex => will be valid. But vertex will change: new elements will be appended to it.
    Returns: new vertex and faces.
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
        [1., 0., 0.],  # 0
        [0., 1., 0.],  # 1
        [0., 0., 1.],  # 2
        #todo: remove unnecessary points

        [0.5/(1.+DIP), 0.5/(1.+DIP), DIP/(1.+DIP)],  # 3
        [DIP/(1.+DIP), 0.5/(1.+DIP), 0.5/(1.+DIP)],  # 4
        [0.5/(1.+DIP), DIP/(1.+DIP), 0.5/(1.+DIP)]   # 5
        ])  # .transpose()

    global trace_subdivided_facets
    trace_subdivided_facets = []

    # Allocate space for the new faces and vertex
    provisional_new_vertex_count = 3*len(tobe_subdivided_face_indices)
    provisional_new_faces_count = 3*len(tobe_subdivided_face_indices)
    nvertex_old = vertex_old.shape[0]
    nfaces_old = faces_old.shape[0]
    new_vertex = np.zeros((nvertex_old + provisional_new_vertex_count, 3), dtype=float)
    new_faces = np.zeros((nfaces_old + provisional_new_faces_count, 3), dtype=int)
    new_vertex[:nvertex_old, :] = vertex_old
    new_faces[:nfaces_old, :] = faces_old

    presubdivision_edges = []
    redundancy_counter = 0

    new_vertex_counter = nvertex_old
    new_faces_counter = nfaces_old
    for subdiv_i in range(len(tobe_subdivided_face_indices)):

        fi = tobe_subdivided_face_indices[subdiv_i]
        triangle_old = vertex_old[faces_old[fi, :], :]  # numvertex x 3
        assert triangle_old.shape == (3, 3)

        # output: v345_xyz
        _vxyz_0123 = np.dot(subdiv_vert_matrix, triangle_old)  # not efficient
        assert _vxyz_0123.shape == (6, 3)
        v345_xyz = _vxyz_0123[3:6, :]  # only pick the new ones
        del _vxyz_0123

        original_face_index = fi
        e0 = new_faces[original_face_index, [0, 1]]
        e1 = new_faces[original_face_index, [1, 2]]
        e2 = new_faces[original_face_index, [2, 0]]
        presubdivision_edges .append(tuple(e0.tolist()))
        presubdivision_edges .append(tuple(e1.tolist()))
        presubdivision_edges .append(tuple(e2.tolist()))

        e012 = np.vstack((e0, e1, e2))
        e012.sort(axis=1)
        BB = np.array([1L, B], dtype=np.int64)
        all_edges_triples = np.dot(e012, BB)
        assert all_edges_triples.dtype == np.int64
        assert all_edges_triples.size == 0 or np.min(all_edges_triples) >= 0
        assert np.max(new_faces, axis=None) < B

        # avoid becasue it is redundant
        avoid_which = np.zeros((3,), dtype=np.bool) + False

        idx_counter = new_vertex_counter
        # actual_mapped_midvertex
        actual_3_vertex = np.zeros((3,), dtype=np.int64)
        for i in range(3):
            if all_edges_triples[i] in midpoint_map:
                avoid_which[i] = True
                actual_3_vertex[i] = midpoint_map[all_edges_triples[i]]
                # mapped_midvertex[i] = -1  # for debug
                redundancy_counter += 1
            else:
                assert avoid_which[i] == False
                # x = mapped_midvertex[i]  # wrong!
                #x = idx_counter
                midpoint_map[all_edges_triples[i]] = idx_counter
                idx_counter += 1
                actual_3_vertex[i] = idx_counter - 1  # the new vertex

        use_which = np.logical_not(avoid_which)
        n1 = new_vertex_counter
        n2 = idx_counter
        assert n2 == n1 + np.sum(use_which)
        new_vertex_counter = n2
        new_vertex[n1:n2, :] = v345_xyz[use_which, :]

        # Output: mini_faces
        # adding new vertex and facets
        # indices of the original vertex.
        _v345 = actual_3_vertex.tolist()
        _v012 = faces_old[fi, :].tolist()
        # facet's vertex indices
        _v012345 = np.array(_v012 + _v345, dtype=int)
        _mini_faces_l = [[0, 3, 5], [3, 1, 4], [5, 4, 2], [3, 4, 5]]  # 0,3,1,4,2,5
        mini_faces = _v012345[np.array(_mini_faces_l, dtype=int)]
        del _v012
        del _v345
        del _v012345

        # facet's vertex indices
        new_faces[original_face_index, :] = mini_faces[0, :]
        new_faces[new_faces_counter:(new_faces_counter+3), :] = mini_faces[1:(1+3), :]
        assert mini_faces.shape[0] == (1+3)
        trace_subdivided_facets += range(new_faces_counter, (new_faces_counter+3)) + [fi]  # include the face which reuses the old face's index
        # trace_subdivided_facets will contain indices of faces

        new_faces_counter += 3

        if subdiv_i % 100 == 0:
            print subdiv_i , "       \r", ;import sys; sys.stdout.flush()

    assert new_vertex.shape[0] - new_vertex_counter == redundancy_counter
    new_vertex = new_vertex[:new_vertex_counter, :]
    # quick_vis(noisy(new_vertex, 0.05), new_facets, range(new_facets.shape[0]))
    assert np.max(new_faces.ravel()) < new_vertex.shape[0]
    assert new_vertex.shape[0] == new_vertex_counter
    assert new_faces.shape[0] == new_faces_counter
    assert provisional_new_vertex_count+nvertex_old-redundancy_counter == new_vertex_counter, "vector consistency"
    assert provisional_new_faces_count+nfaces_old == new_faces_counter, "face consistency"
    assert len(trace_subdivided_facets) == 0 or np.max(np.array(trace_subdivided_facets)) < new_faces_counter

    return new_vertex, new_faces, presubdivision_edges


def compute_facets_subdivision_curvatures(vertex, faces, iobj, curvature_epsilon):
    """ Deviation of Mesh from object gradients """

    face_areas, face_normals = compute_triangle_areas(vertex, faces, return_normals=True)

    nf = faces.shape[0]
    assert face_areas.shape == (nf,)
    assert face_normals.shape == (nf, 3)

    assert np.all(np.logical_not(np.isnan(face_areas[np.logical_not(np.isnan(np.linalg.norm(face_normals, axis=1)))])))

    degenerate_faces = np.isnan(face_areas)

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

        n = face_normals[fi, :]  # n: (3,)

        triangle = vertex[faces[fi, :], :]  # numvertex x 3
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
            e = face_areas[fi] * np.sum(1. - np.abs(np.dot(n, mm))) / 4.  # sum(,x4)

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
def comparison_vertex_new_vertex(old_vertex, new_vertex):
    THL = 10. ** (0.001)  # find experimentally
    nvertex = old_vertex.shape[0]
    assert old_vertex.shape == new_vertex.shape

    for i in range(nvertex):
        if np.abs(np.linalg.norm(old_vertex[i, :] - new_vertex[i, :])) > THL:
            new_vertex[i, :] = old_vertex[i, :]

    return new_vertex


def simple_histogram(c, title=None, special_values=[]):
    import matplotlib.pyplot as plt
    # special_values

    plt.hist(c, 20*10)
    if title is not None:
        plt.title(title)
    plt.show()


def compute_facets_curvatures_vectorized(vertex, faces, iobj):
    face_areas, face_normals = compute_triangle_areas(vertex, faces, return_normals=True)

    nf = faces.shape[0]
    assert face_areas.shape == (nf,)
    assert face_normals.shape == (nf, 3)
    assert np.all(np.logical_not(np.isnan(face_areas[np.logical_not(np.isnan(np.linalg.norm(face_normals, axis=1)))])))
    # some edges are repeated
    degenerate_faces = np.isnan(face_areas)

    nn = np.isnan(face_normals[np.logical_not(degenerate_faces), :])
    assert np.all(np.logical_not(np.isnan(face_areas.ravel()))), "facet_areas: never nan. But can be zero."
    assert np.all(np.isnan(face_areas[degenerate_faces]))
    assert np.all(np.logical_not(np.isnan(face_areas[np.logical_not(degenerate_faces)])))
    assert np.all(np.isnan(face_normals[degenerate_faces, :]))


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
        ])

    triangles = vertex[faces[:, :], :]  # nf x 3 x 3

    if True:
        nf = triangles.shape[0]

        subdivmat = np.dot(centroidmaker_matrix, subdiv_vert_matrix)  # 4x3

        q = np.tensordot(triangles, subdivmat, axes=(1, 1))  # nf x 3 x 4
        del triangles
        q = np.swapaxes(q, 1, 2)  # nf x 4 x 3
        assert q.shape[0]*q.shape[1] == nf*4
        q = q.reshape(q.shape[0]*q.shape[1], 3)

        assert q.shape == (nf*4, 3)
        check_vector3_vectorized(q)
        mm = - iobj.implicitGradient(q); del q
        assert mm.shape == (nf*4, 3)
        mmt = mm.reshape(nf, 4, 3); del mm  # nfx4x3

        mmt_norm = np.linalg.norm(mmt, axis=2, keepdims=True)
        mmt_norm[mmt_norm < 0.00001] = 1.
        import ipdb; ipdb.set_trace()
        mmt_hat = mmt / mmt_norm ; del mmt; del mmt_norm

        n_norm = np.linalg.norm(face_normals, axis=1)
        n_norm[n_norm < 0.00001] = 1.
        n_hat = face_normals / n_norm; del n_norm

        assert mmt_hat.shape[2] == 3
        assert n_hat.shape[1] == 3

        nm = np.sum(n_hat[:, np.newaxis, :] * mmt_hat, axis=2)  # Fx1x3 * Fx4x3 -> Fx4

        curvatures_array = face_areas * (1. - np.sum(np.abs(nm), axis=1) / 4.)
        curvatures_array[degenerate_faces] = 0.
        bad_facets_count = np.sum(degenerate_faces)

        return curvatures_array, bad_facets_count

    curvatures_array = np.zeros((nf,))
    for fi in range(nf):
        n = face_normals[fi, :]  # n: (3,)
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
            print "WARNING: degenerate triangle", fi, " = ", faces[fi, :]
        else:
            assert not degenerate_faces[fi]

        triangle = triangles[fi, :, :]  # numvertex x 3
        assert not np.any(np.isnan(triangle.ravel()))
        assert triangle.shape == (3, 3)  # nv=3
        assert not np.any(np.isnan(triangle.ravel()))
        assert not np.any(np.isnan(subdiv_vert_matrix.ravel()))

        subdivmat = np.dot(centroidmaker_matrix, subdiv_vert_matrix)
        m0123 = np.dot(subdivmat, triangle)
        assert m0123.shape == (4, 3)
        assert not np.any(np.isnan(m0123.ravel()))
        subdiv_centroids = m0123

        assert not np.any(np.isnan(subdiv_centroids.ravel()))

        mm = - iobj.implicitGradient(subdiv_centroids)
        assert mm.shape == (4, 3)
        nn = np.linalg.norm(mm, axis=1)
        nn_tile = np.tile(nn[:, np.newaxis], (1, 3))  # mm: 4 x 3
        if mesh_correction:
            nn_tile[nn_tile < 0.00000001] = 100000.
        mm = mm / nn_tile
        mm = mm.transpose()  # 3x4

        if np.any(np.isnan(n)):
            e = 0.
        else:
            e = face_areas[fi] * np.sum(1. - np.abs(np.dot(n, mm))) / 4.  # sum(,x4)
            if e < 0:
                set_trace()
            assert e >= 0

        curvatures_array[fi] = e

        if not mesh_correction:  # NAN is allowed (and used) for mesh correction
            assert not np.isnan(e)

        if fi % 100 == 0:
            print fi, "*   \r", ;import sys; sys.stdout.flush()

    # The only reason for NaN should be the exactly-zero gradients
    # assert np.sum(np.isnan(curvatures_array)) == 0
    l = curvatures_array[np.logical_not(np.isnan(curvatures_array))].tolist()
    l.sort()
    print "curvature range: min,max = ", l[0], l[-1]   # 3.80127650325e-08, 0.0240651184551
    bad_faces_count = np.sum(degenerate_faces)

    return curvatures_array, bad_faces_count


def propagated_subdiv(faces, subdivided_edges):
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

    assert subdivided_edges.shape[1] == 2
    subdivided_edges.sort(axis=1)

    BB = np.array([1L, B], dtype=np.int64)
    subdivided_edges_codes = np.dot(subdivided_edges, BB)  # can be large
    assert subdivided_edges_codes.dtype == np.int64
    assert subdivided_edges_codes.size == 0 or np.min(subdivided_edges_codes) >= 0
    assert np.max(faces, axis=None) < B

    fc0 = faces.copy()
    f0 = faces[:, np.newaxis, [0, 1]]
    f1 = faces[:, np.newaxis, [1, 2]]
    f2 = faces[:, np.newaxis, [0, 2]]
    f012 = np.concatenate((f0, f1, f2), axis=1)
    f012.sort(axis=2)
    assert np.all(fc0.ravel() == faces.ravel())  # no copy() needed
    all_edges_codes = np.dot(f012, BB)  # *x3

    # now look for subdivided_edges_codes in all_edges_codes

    # todo: refactor: intersec seems to be not used anymore
    intersec = np.intersect1d(all_edges_codes, subdivided_edges_codes)
    # x_ indicates those edges taht are in the mesh's edges, hence remain there.
    # x_ is the bottleneck (both in terms of data and in terms of performance)
    x_ = np.lib.arraysetops.in1d(all_edges_codes, intersec)  # elements of A, A.ravel[x_], that are in B
    assert np.prod(all_edges_codes.shape) == x_.size
    assert np.sum(x_) == intersec.size   # 417
    # assert each_of all_edges_codes.ravel()[x_] in intersec
    edges_need_subdivision = all_edges_codes.ravel()[x_]  # all edges_need_subdivision are in intersec

    sides_booleans_Fx3 = x_.reshape(all_edges_codes.shape)  # *x3
    numsides = np.sum(sides_booleans_Fx3, axis=1)  # numsides_needsubdivision: number of sides that need subdivision. index=face index
    assert sides_booleans_Fx3.shape == (faces.shape[0], 3)

    # now I need those all_edges_codes (i.e. sides_booleans_Fx3) for which numsides==1
    # sides_1 = sides_booleans_Fx3[numsides==1, :] for c==1
    propag_dict = {}
    for c in range(1, 4):
        # Range starts with 1 because we only propagate triangles with subdivided 1,2,3 sides.
        idx = np.nonzero(numsides == c)[0]
        propag_dict[c] = idx
        del idx

    return propag_dict, edges_need_subdivision


def get_edge_code_triples_of_mesh(faces):
    """ Returns an array of (F)x(3), containing the 'edge codes' of sides of the faces of a mesh.
    There are 3 sides for each face.
    An 'edge code' is a long integer (int64) v1+B*v2 where v1,v2 are the indices of the ends (vertex) of the edge, where v1<v2."""
    e0 = faces[:, np.newaxis, [0, 1]]
    e1 = faces[:, np.newaxis, [1, 2]]
    e2 = faces[:, np.newaxis, [2, 0]]   # np view
    e012 = np.concatenate((e0, e1, e2), axis=1)  # n x 3 x 2
    assert e012.base is None  # make sure it's not a view of faces
    e012.sort(axis=2)
    BB = np.array([1L, B], dtype=np.int64)
    all_edges_triples = np.dot(e012, BB)  # n x 3
    assert all_edges_triples.dtype == np.int64
    assert all_edges_triples.size == 0 or np.min(all_edges_triples) >= 0
    assert np.max(faces, axis=None) < B
    assert all_edges_triples.size == 0 or np.min(all_edges_triples) >= 0
    assert all_edges_triples.shape == (faces.shape[0], 3)
    return all_edges_triples


def subdivide_1to2_multiple_faces(faces, edges_with_1_side, midpoint_map, careful_for_twosides=True):
    """list_edges_with_1_side contains the edges only. The face should be extracted in this function.
    returns: faces.
    careful_for_twosides: whe ntwo sides are being asked for subdivision. """
    # todo: copy some code from propagated_subdiv()
    # check which of these edges still exist in faces. (Each should be there only once. In this context.)
    # Some edges_with_1_side may not be in facets. They are already subdivided twice.
    # remove them and add more.
    # refactor the code copied from propagated_subdiv() into function
    # need to also get the new points. oops!! damn.
    # There is a guarantee that all faces that the edges_with_1_side belong to, have exactly one edge from this list.
    # Note that these edges should be already subdivided


    # yes of course all of them are there in it
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

    all_edges_triples = get_edge_code_triples_of_mesh(faces)

    all_edge_triples_ravel = all_edges_triples.ravel()  # is a view
    # print all_edges_triples.shape
    assert all_edges_triples.shape[1] == 3

    assert np.all(all_edge_triples_ravel.reshape(-1, 3) == all_edges_triples, axis=None)

    # intersec = np.intersect1d(all_edge_triples_ravel, edges_with_1_side)

    # index_of_edges_that_subdiv2

    x_ = np.lib.arraysetops.in1d(all_edge_triples_ravel, edges_with_1_side)  # elements of A, A.ravel[x_], that are in B

    x_1x3 = x_.reshape(3, -1)
    assert np.all(x_1x3.ravel() == x_, axis=None)

    x3__b_Fx3 = x_.reshape(-1, 3)

    # how many sides are requested to be subdivided
    facemultiplicity = np.sum(x3__b_Fx3, axis=1)

    # Dont want to subdivide 1->2
    # bad2 = np.all(np.sum(x_.reshape(3, -1), axis=0) > 1)  # bug!
    bad2 = np.nonzero(facemultiplicity > 1)[0]
    if careful_for_twosides:
        if bad2.size > 0:
            print midpoint_map
            print bad2
            print faces[bad2, :]
            print edges_with_1_side
            set_trace()
        assert bad2.size == 0
    del bad2

    if careful_for_twosides:
        assert np.all(facemultiplicity <= 1)
    if careful_for_twosides:
        if not np.all(facemultiplicity <= 1):
            # print facemultiplicity.tolist()
            # a = facemultiplicity
            # print np.nonzero(a > 1)
            print "FAILED"
        assert np.all(facemultiplicity <= 1)
    # print "THIS FAILS"
    del x3__b_Fx3

    # indices of all edges
    face3_idx = np.nonzero(x_)[0]
    assert np.ndim(face3_idx) == 1

    # todo(refactor): use np.argwhere()
    idx_xy = np.unravel_index(face3_idx, all_edges_triples.shape)
    # idx_xy is a tuple
    # Triangles subject to be subdivided:
    problem_face_idx = idx_xy[0]
    problem_side_idx = idx_xy[1]
    # assert np.all(problem_face_idx < 3)

    def has_repeats(x):
        y = x.copy()
        y.sort()
        return np.any(np.diff(y) == 0)

    # has repeats, becasue [2]is not resolved before
    # if has_repeats(problem_face_idx):
    #    intersec = np.intersect1d(all_edge_triples_ravel, edges_with_1_side)
    #    print intersec
    #    print facets
    #    exit()

    if careful_for_twosides:
        assert not has_repeats(problem_face_idx), "triangles subject to be subdivided"
    for ii in range(problem_face_idx.size):
        # assert all_edge_triples_ravel[problem_face_idx[ii], problem_side_idx[ii]] in edges_with_1_side
        assert all_edges_triples[problem_face_idx[ii], problem_side_idx[ii]] in edges_with_1_side
        assert problem_side_idx[ii] < 3
    # all problem_face_idx should be removed

    # The subdivided edge is between v1 and v2.
    vert_idx_1 = problem_side_idx  # The problem_side will be between (v1,v2) vertex. vert_idx_1 is not a vertex index but it is a vertex index within a face i.e. in (0,1,2).
    vert_idx_2 = (problem_side_idx + 1) % 3
    vert_idx_3 = (problem_side_idx + 2) % 3
    v1 = faces[problem_face_idx, vert_idx_1]
    v2 = faces[problem_face_idx, vert_idx_2]
    v3 = faces[problem_face_idx, vert_idx_3]

    # The sides (vertex pairs) that need to be subdivided
    subdivedges_vertex_pairs = np.vstack((v1, v2))  # size: 2 x F

    # subdivedges_vertex_pairs never actually used apart from assertion tests.

    # edge_s_codes = A flat array of all the edge codes (For the sides that should be replaced with the sibdivided ones)
    # ?????????????
    edge_s_codes = all_edge_triples_ravel[x_]  # Intersection from actual edges in mesh and edges requested to get removed/subdivided.
    # subdivedges_vertex_pairs: those edges that*

    # observation: edge_s_codes is (up to morphism) a subset of, but not equal to, subdivedges_vertex_pairs
    if careful_for_twosides:
        assert np.unique(edge_s_codes).size == subdivedges_vertex_pairs.shape[1]  # before applying unique

    # if can tolerate two sides:
    # if not careful_for_twosides:
    # edge_s_codes = np.unique(edge_s_codes)
    if careful_for_twosides:
        assert np.unique(edge_s_codes).size == edge_s_codes.size
    # edge_s_codes are unique but subdivedges_vertex_pairs are not unique

    #####################################################################################################################
    tesort = subdivedges_vertex_pairs.T.copy()
    tesort.sort(axis=1)
    eid9 = np.dot(tesort, np.array([1, B], dtype=np.int64)).copy()
    assert eid9.dtype == np.int64
    assert eid9.size == 0 or np.min(eid9) >= 0
    assert np.max(faces, axis=None) < B
    eid9.sort()
    eid10 = edge_s_codes.copy()
    eid10.sort()
    assert np.all(eid10 == eid9)
    del subdivedges_vertex_pairs

    # map one-to-one between: edge_s_codes and midpoints_third_vertex and (v1, v2, ...)
    midpoints_third_vertex = np.array(map(lambda edgecode: midpoint_map[edgecode], edge_s_codes), dtype=np.int64)
    # print midpoints_third_vertex
    if midpoints_third_vertex.size == 0:
        return faces

    assert isomorphic(midpoints_third_vertex, v1)
    # (v1,v2,v3) -> (v1, midpoints_third_vertex, v3) + (midpoints_third_vertex, v2, v3)
    # (v1,v3, midpoints_third_vertex),  (v2,v3, midpoints_third_vertex)
    new_faces1 = np.vstack(((v1, v3, midpoints_third_vertex))).T  # axis is 0. .T.size = N x 3
    new_faces2 = np.vstack(((v2, v3, midpoints_third_vertex))).T
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
    f_rm = np.delete(faces, todelete, axis=0)
    appended_faces = np.concatenate((f_rm, new_faces), axis=0)

    return appended_faces


def do_subdivision(vertex, faces, iobj, curvature_epsilon, randomized_probability=1.):
    assert not np.any(np.isnan(faces.ravel()))
    assert not np.any(np.isnan(vertex.ravel()))  # fails

    print "computing curvatures"; sys.stdout.flush()
    # curvatures, bad_facets_count = compute_facets_curvatures_vectorized(vertex, facets, iobj)
    curvatures, bad_facets_count = compute_facets_subdivision_curvatures(vertex, faces, iobj, curvature_epsilon)
    print "computing curvatures done."; sys.stdout.flush()

    assert np.sum(np.isnan(curvatures)) == 0, "NaN"
    curvatures[np.isnan(curvatures)] = 0  # treat NaN curvatures as zero curvature => no subdivision

    which_faces = np.arange(faces.shape[0])[curvatures > curvature_epsilon]
    if randomized_probability < 1.:
        n0 = which_faces.shape[0]
        m0 = int(np.ceil(float(randomized_probability)*float(n0)))
        ridx = np.random.choice(n0, m0, replace=False)
        assert len(ridx) == m0
        which_faces = which_faces[ridx]

    print which_faces.shape
    print "applying subdivision on %d triangles." % (int(which_faces.shape[0]))
    midpoint_map = {}
    vertex2, faces2, presubdivision_edges = subdivide_multiple_faces_new(vertex, faces, which_faces, midpoint_map)
    global trace_subdivided_faces  # third implicit output

    list_edges_with_1_side = []
    while True:
        propag_dict, edges_which_in1 = propagated_subdiv(faces2, presubdivision_edges)
        faces_with_2_or_3_sides = np.concatenate((propag_dict[2], propag_dict[3]), axis=0)
        # what if those faces dont exist anymore in the next round?
        list_edges_with_1_side += [edges_which_in1]
        # print facets_with_2_or_3_sides.shape
        if faces_with_2_or_3_sides.size == 0:
            break
        vertex2, faces2, old_edges2 = subdivide_multiple_faces_new(vertex2, faces2, faces_with_2_or_3_sides, midpoint_map)
        presubdivision_edges += old_edges2  # bug fixed!

    # Finished with 2 or 3 sides.

    # Now 1 side:
    # Append all the lists in list_edges_with_1_side
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
    # todo: if length zero dont do it
    assert edges_with_1_side.size == 0 or np.min(edges_with_1_side) > 0

    faces2 = subdivide_1to2_multiple_faces(faces2, edges_with_1_side, midpoint_map)

    ###################
    # check_degenerate_faces(vertex2, facets2, "assert")
    # build_faces_of_faces(facets2)

    print("Subdivision applied.");sys.stdout.flush()
    return vertex2, faces2


@profile
def demo_combination_plus_qem():
    """ Now with QEM """
    curvature_epsilon = 1. / 1000.  # a>eps  1/a > 1/eps = 2000
    # curvature_epsilon = 1. / 10000.
    VERTEX_RELAXATION_ITERATIONS_COUNT = 3

    NUMBER_OF_ITERATION_GLOBAL = 2

    from example_objects import make_example_vectorized
    object_name = "cube_with_cylinders"  # "sphere_example" #or "rcube_vec" work well #"ell_example1"#"cube_with_cylinders"#"ell_example1"  " #"rdice_vec" #"cube_example"
    iobj = make_example_vectorized(object_name)

    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-3, +5, 0.2)
    if object_name == "cube_with_cylinders" or object_name == "twist_object" or object_name == "french_fries" or object_name == "rdice_vec" or object_name == "cyl4" or object_name == "rods" or object_name == "bowl_15_holes":
        VERTEX_RELAXATION_ITERATIONS_COUNT = 1

    if object_name == "cyl4":
        (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-32 / 2, +32 / 2, 1.92 / 4.0)

    elif object_name == "french_fries" or object_name == "rods":
        (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-3, +5, 0.11)  # 0.05

    elif object_name == "bowl_15_holes":
        (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-3, +5, 0.10)

    elif object_name == "cyl3":
        (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-32 / 2, +32 / 2, 1.92 / 4.0)

    elif object_name == "cyl1":
        (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-16, +32, 1.92 * 0.2 * 10 / 2.0)

    elif object_name == "cyl2":
        (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-32, +32, 1.92 / 4.0 * 1.5 / 1.5)

    from stl_tests import make_mc_values_grid
    gridvals = make_mc_values_grid(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE, old=False)
    vertex, faces = vtk_mc(gridvals, (RANGE_MIN, RANGE_MAX, STEPSIZE))
    print("MC calculated."); sys.stdout.flush()

    # **********************
    # ******* n-loop *******
    # **********************
    for i in range(NUMBER_OF_ITERATION_GLOBAL):

        # **********************
        # * vertex relaxation *
        # **********************
        for i in range(VERTEX_RELAXATION_ITERATIONS_COUNT):
            vertex, faces_not_used, centroids = process2_vertex_resampling_relaxation(vertex, faces, iobj)
            assert not np.any(np.isnan(vertex.ravel()))  # fails
            print("vertex relaxation applied.");sys.stdout.flush()

        # ************************
        # ** projection and qem **
        # ************************
        average_edge = compute_average_edge_length(vertex, faces)

        old_centroids = np.mean(vertex[faces[:], :], axis=1)

        new_centroids = old_centroids.copy()

        set_centers_on_surface__ohtake_v3s(iobj, new_centroids, average_edge)

        vertex_neighbours_list = mesh_utils.make_neighbour_faces_of_vertex(faces)
        centroid_gradients = compute_centroid_gradients(new_centroids, iobj)

        new_vertex_qem = \
            vertex_apply_qem3(vertex, faces, new_centroids, vertex_neighbours_list, centroid_gradients)

        vertex_before_qem = vertex


        # ************************
        # **  subdivision **
        # ************************
        pre_subdiv_vf = (new_vertex_qem, faces)
    #    total_subdivided_facets = []

        (vertex, faces) = new_vertex_qem, faces  # ??

        print "subdivision:"
        vertex, faces = do_subdivision(vertex, faces, iobj, curvature_epsilon)
        global trace_subdivided_facets  # third implicit output

        print "subdivision done."

        for use_wireframe in [True, False]:

            display_simple_using_mayavi_2([(vertex, faces), (vertex, faces), (pre_subdiv_vf[0], pre_subdiv_vf[1]), ],
                  pointcloud_list=[],
                  mayavi_wireframe=[False, use_wireframe, True, ], opacity=[0.2, 1, 0.3], gradients_at=None, gradients_from_iobj=None,
                  minmax=(RANGE_MIN,RANGE_MAX))


    # ***************************
    # ******* And finally *******
    # ***************************

    for i in range(VERTEX_RELAXATION_ITERATIONS_COUNT):
        vertex, faces_not_used, centroids = process2_vertex_resampling_relaxation(vertex, faces, iobj)
        assert not np.any(np.isnan(vertex.ravel()))  # fails
        print("vertex relaxation applied.");sys.stdout.flush()

    # ************************
    # ** projection and qem **
    # ************************
    average_edge = compute_average_edge_length(vertex, faces)

    old_centroids = np.mean(vertex[faces[:], :], axis=1)

    new_centroids = old_centroids.copy()

    set_centers_on_surface__ohtake_v3s(iobj, new_centroids, average_edge)

    vertex_neighbours_list = mesh_utils.make_neighbour_faces_of_vertex(faces)
    centroid_gradients = compute_centroid_gradients(new_centroids, iobj)

    new_vertex_qem = \
        vertex_apply_qem3(vertex, faces, new_centroids, vertex_neighbours_list, centroid_gradients)

    vertex_before_qem = vertex

    display_simple_using_mayavi_2([(vertex_before_qem, faces), (new_vertex_qem, faces), ],
       pointcloud_list=[],
       mayavi_wireframe=[False, False], opacity=[0.4*0, 1, 0.9], gradients_at=None, separate=False, gradients_from_iobj=None,
       minmax=(RANGE_MIN, RANGE_MAX))

# from timeit import default_timer as dtimer


if __name__ == '__main__':

    # timer_t1s = dtimer()

    demo_combination_plus_qem()  # subdivision + projection + qem
