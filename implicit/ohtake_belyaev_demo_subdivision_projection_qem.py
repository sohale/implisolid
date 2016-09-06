from ipdb import set_trace
import profile_support
from vtk_mc import vtk_mc

import sys
import math

import numpy as np
from basic_types import check_vector4_vectorized, normalize_vector4_vectorized

def bisection_3_standard(iobj, p1, p2, f1, f2, MAX_ITER):
    TH1 = 0.001
    TH3 = 0.001
    #for i in range(p1.shape[0])
    assert p1.shape[0] == 1
    assert p2.shape[0] == 1

    assert f1 < 0
    assert f1*f2 < 0, "Opposite signs required"
    for j in range(MAX_ITER):
        if np.linalg.norm(p1-p2) < TH1:
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
    #print "Convergence of the bisection did not happen"
    return None, MAX_ITER


#@profile
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

        if f3 < 0.:
            if VERBOSE:
                print "A",
            (p1, f1) = (p3, f3)
        else:
            if VERBOSE:
                print "B",
            (p2, f2) = (p3, f3)

        dt = p2 - p1

        if np.linalg.norm(dt) < TH2:
            return True, p3, None, j

    if VERBOSE:
        print "Convergence based on the proportional method did not happen"

    p,ni = bisection_3_standard(iobj, p1, p2, f1, f2, MAX_ITER*4)
    if p is not None:
        if VERBOSE:
            print "Bisection converged "
        return True, p, None, 0  # very unlikely
    else:
        return False, p1, p2, -ni


VERBOSE = False

#@profile
def search_near_ohtake_old(iobj, start_x, direction, lambda_val, MAX_ITER):  # max_dist
    """Returns either the point, or None, if not found. lambda_val is the expected distance from the surface.
    The resommended value is half of the average edge length, but a better way is half of the MC'step size (because the expected error is half of the MC grid voxel size).
    Remeber: we may be on a totally irrelevant direction here.
    'direction' should be normalised. lambda*direction is used. Note that lambda is negated by default.
    Note: along_1d mode is in fact the same. It's just initialises direction=gradient(start_x).
    Does both searchNearPoint1D and searchNearPoint()
    :param direction: description
    @param np.array direction
    """

    TH1 = 0.001

    TH2_L = 0.00001

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

    exit_A = False
    while True:

        for j in range(MAX_ITER):
            if not along_1d:
                direction = iobj.implicitGradient(start_x)  ## what?! why start_x ??
                dn = np.linalg.norm(direction)
                if dn>0.0:  # 00000001:
                    direction = direction/dn
                else:
                    pass

            p2 = p2 + lambda_ * direction
            p2[:, 3] = 1
            f2 = iobj.implicitFunction(p2)
            eval_count += 1

            if f1*f2 < 0.0:
                # (A)
                exit_A = True
                break
            p1 = p2

        else:

            lambda_ = lambda_ / 2.

            if np.abs(lambda_) < TH2_L:
                print "(B)", eval_count
                return None   #

            p1 = start_x
            assert f0 * f1 >= 0
            f1 = f0
            p2 = p1

        if exit_A:  # f1*f2 < 0.0:
            break

    assert f1*f2 < 0.0
    if f1 > 0:
        (p1, p2) = (p2, p1)
        (f1, f2) = (f2, f1)


    converged, p1, p2, iter_log = bisection_prop_2(iobj, p1, p2, f1, f2, MAX_ITER/2)
    assert f1 < 0
    assert f2 > 0

    if converged:
        assert p2 is None
        return p1
    else:
        return None

def get_data_base(arr):
    """For a given Numpy array, finds the
    base array that "owns" the actual data."""
    base = arr
    while isinstance(base.base, np.ndarray):
        base = base.base
    return base

def arrays_share_data(x, y):
    return get_data_base(x) is get_data_base(y)

@profile
def search_near_ohtake(iobj, start_x, lambda_val, MAX_ITER):  # max_dist
    """Vectorized function of the search_near_ohtake_old who take a vector (start_x)
     and return the vector p1_found if we found the centroids"""

    #lambda should be ~ expected distance?  (that's why it should be half of the average edge size, i.e. half of the MC step size)
    TH1 = 0.001
    TH2_L = 0.00001

    #start_x_new = start_x.copy()

    #p1_found = np.ndarray(start_x.shape)
    p1_found = np.zeros(start_x.shape)
#    p2 = np.ndarray(start_x.shape)
#    f0 = np.ndarray(start_x.shape[0])
#    f2 = np.ndarray(start_x.shape[0])
    direction = np.ndarray(start_x.shape)

    p1 = start_x.copy() #adding the copy : solve one problem of the code
    #p1 = start_x
    f1 = iobj.implicitFunction(p1)

    lambda_ = lambda_val * (-np.sign(f1))

#    counting = 0

    for i in range(start_x.shape[0]):
        xxv = start_x[i,:].reshape(1,4)
        check_vector4_vectorized(xxv)
        direction[i,:] = (iobj.implicitGradient(xxv)) #no need for copy here!
        dn = np.linalg.norm(direction[i,0:3]) #calculation of the normalization

        if  dn > 0.0:
            direction[i,:] = (direction[i,:]/dn)
            direction[i,3] = 1

#        f0 = f1[i]
        eval_count = 0
        eval_count += 1


        if math.fabs(f1[i]) < TH1:
            #print i, "Out"
            p1_found[i,:] = p1[i,:]

        else:

        #    negative_f1 = -1 if f1[i] < 0. else 1
        #    lambda_ = lambda_val * (-negative_f1)  # negation of input is bad practice

            p2 = p1[i,:].copy()
            exit_A = False
            not_found = False
            number_of_iteration = 0
            s1 = start_x[i,:].copy()

            while True:

                # (C) jumps back here.
                for j in range(MAX_ITER):
                    number_of_iteration += 1

                    xxv2 = start_x[i,:].copy().reshape(1,4)
                    check_vector4_vectorized(xxv2)
                    dir1 = iobj.implicitGradient(xxv2) #we have to give a copy of start_x to avoid problems
                    dir1 = dir1[0,:]
                    dn1 = np.linalg.norm(dir1[0:3].copy())

                    if  dn1 > 0.0:
                        dir1 = dir1/dn1
                    dir1[3] = 1
                    dir2 = direction[i,:]


                    # if not np.allclose(s2, s1):
                    #     print i, "222"
                    #     print start_x_new[i].shape, start_x_new[i].dtype
                    #     print s1.shape, s1.dtype
                    if not np.allclose(dir1,dir2):
                        print i,"patate!"
                        print dir1, dir2
                        print dir1/dir1[0], dir2/dir2[0]
                        print "ALLCLOSE failed"
                        exit()


                    p2 = p1[i,:] + lambda_[i] * dir2
            #        p2[i,:] = p1[i,:] + lambda_[i] * direction[i]

                    p2[3] = 1
                    f2 =iobj.implicitFunction(p2.reshape(1,4))
                    eval_count += 1

                    if f1[i]*f2 < 0.0:
                        exit_A = True #(A)
                        break

                    p1[i,:] = p2

                    if number_of_iteration == MAX_ITER:
                        number_of_iteration = 0
                        lambda_[i] = lambda_[i] / 2.

                    #    print lambda_[i], TH2_L

                        if np.abs(lambda_[i]) < TH2_L:
                            print i, "(B)", eval_count
                            p1_found[i,3] = 0
                            not_found = True
                            break
                        else:
                            p1[i,:] = start_x[i,:].copy()
                            # assert f0 * f1[i] >= 0
                            # f1[i] = f0  #This was missing in Ohtake, because the sign of f1 is not expected to change. So I added the assert above.
                            p2 = p1[i,:].copy()

                if exit_A:
                    break
                if not_found: #np.abs(lambda_)< TH2_l
                    break
    #
            if not not_found:

                #import ipdb; ipdb.set_trace()
                assert f1[i]*f2 < 0.0
                if f1[i] > 0:
                    (p1[i,:], p2) = (p2, p1[i,:])
                    (f1[i], f2) = (f2, f1[i])

                assert f1[i] < 0
                assert f2 > 0
                assert math.fabs(f1[i]) < 1000, str(p1[i,:])
                assert math.fabs(f2) < 1000, str(p2)

                converged, p1[i,:], p2, iter_log = bisection_prop_2(iobj, p1[i,:].reshape(1,4), p2.reshape(1,4), f1[i], f2, MAX_ITER/2)
                assert f1[i] < 0
                assert f2 > 0

                if converged:
                    p1_found[i,:] = p1[i,:].copy()

                else:
                    p1_found[:,3] = 0
                    p1_found[i,3] = 0

                #    print p1_found[i,:]
        #    print direction[i,:], np.linalg.norm(direction[i,:])
        # if f1[i] == f0:
        #     counting += 1
    #    print p1_found[i] - start_x[i]
#    print counting
    return p1_found

@profile
def set_centers_on_surface_ohtake(iobj, centroids, average_edge):
#    here we consider that the max_dist is the average_edge and lambda = average_edge/2
#    new function who is a combination of sers_on_surface_ohtake and project_point_bidir_ohtake
    lambda_val = average_edge/2
    check_vector4_vectorized(centroids)
    #definition of the matrix that are gonna be used in the rest of the programm
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

    p1 = search_near_ohtake(iobj, centroids, lambda_val, max_iter)

    for i in range(centroids.shape[0]):

        if np.allclose(p1[i, 3], 1, 0.00000000000001) == True: #make sure that p are found by the program and they respect the condition enforce by check_vector4_vectorized
            p1[i,:].reshape(1,4)
            f1[i] = iobj.implicitFunction(p1[i,:].reshape(1,4))

                # Mirror image: search the opposite way and accept only if it is closer than the best found.
            p2[i,:] = 2*centroids_new[i,:] - p1[i,:] #p2 correspond to S in the paper
            f2[i] = iobj.implicitFunction(p2[i,:].reshape(1,4))
            p[i,:] = p1[i,:]


            if f1[i]*f2[i] < 0:
                direction_3[i,:] = (centroids_new[i,:] - p1[i,:].copy())  # as in Ohtake
                dn_3[i] = np.linalg.norm(direction_3[i,:])
                if dn_3[i] > 0:
                    direction_3[i,:] = direction_3[i,:]/dn_3[i]
                    p3 = search_near_ohtake_old(iobj, centroids[i,:].reshape(1,3), direction_3[i,:].reshape(1,3), lambda_val, max_iter)


                    #no max_dist
                if p3 is not None:
                    if np.allclose(p3[:,3], 1, 0.00000000000001) == True:
                        if np.linalg.norm(centroids[i,:] - p3) > np.linalg.norm(centroids[i,:] - p1[i,:]):
                            p[i,:] = p3

            if np.linalg.norm(centroids[i,:] - p[i,:]) <= average_edge:
                centroids[i,:] = p[i,:]


def display_simple_using_mayavi_(vf_list, pointcloud_list, minmax=(-1,1), mayavi_wireframe=False, opacity=1.0):
    """ opacity can be either a list of a constant """
    from mayavi import mlab

    if type(opacity) is list:
        opacities = opacity  # 1.0
    else:
        opacities = [opacity] + [0.2]*(len(vf_list)-1)  # 1.0, 0.2 #0.1

    i = 0
    for vf in vf_list:
        verts, faces = vf
        if verts is None:
            continue
        if verts.size == 0:
            print("Warning: empty vertices")
            continue
        if faces.size == 0:
            print("Warning: no faces")
            continue
        mlab.triangular_mesh([vert[0] for vert in verts],
                         [vert[1] for vert in verts],
                         [vert[2] for vert in verts],faces,representation="surface" if not mayavi_wireframe else "wireframe",
                         opacity=opacities[i], scale_factor = 100.0)
        i += 1


    color_list = [(1, 0, 0), (0, 0, 0), (1, 1, 0), (0, 0, 1), (0,1,0)]
    i = 0
    for c in pointcloud_list:
        #print c[:,0:3]
        mlab.points3d(c[:, 0], c[:, 1], c[:, 2], color=color_list[i], scale_factor=0.2)
        i+=1
    del i


    #if pointcloud_list[0].shape[0] == pointcloud_list[1].shape[0]:
    #    print "PLOT3d"
    #    c0 = pointcloud_list[0][np.newaxis, :, :]   # 2 x N x 4
    #    c1 = pointcloud_list[1][np.newaxis, :, :]
    #    c01 = np.concatenate((c0, c1), axis=0)
    #    for i in range(c01.shape[1]):
    #        #print c01[:,i,0], c01[:,i,1], c01[:,i,3]
    #        mlab.plot3d(c01[:,i,0],c01[:,i,1],c01[:,i,3])
    #        #exit()
    #        pass
    #    #mlab.plot3d((c01[:,:,0]), (c01[:,:,1]), (c01[:,:,3]))


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

    mlab.show()

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
        separate=True, gradients_at=None, gradients_from_iobj=None, pointsizes=None, pointcloud_opacity=1.):
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
            mlab.points3d(c[:, 0], c[:, 1], c[:, 2], color=color_list[i], scale_factor=pointsizes[i], opacity=pointcloud_opacity )
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
        assert not degenerate_faces[fi]
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
        assert not degenerate_faces[fi]
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

#@profile
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
    n_i = normals[:, 0:3, np.newaxis]
#    print n_i.shape, n_i.T.shape
    p_i = center_array[:, 0:3, np.newaxis]
    # nnt = np.dot(n_i, n_i.T)
#    print repr(n_i) + "+++++++++++++++++++++++++++++++++++++++++++++="
    # import ipdb; ipdb.set_trace()
    A = np.dot(np.reshape(n_i,(normals.shape[0],3)).T, np.reshape(n_i,(normals.shape[0],3)))
#   b = -np.dot(np.diag(np.reshape(n_i,(normals.shape[0],3)),0), np.dot(np.diag(np.reshape(n_i,(normals.shape[0],3)),0), np.diag(p_i.reshape(normals.shape[0],3),0)))
    for i in range(normals.shape[0]):

        assert n_i[i].shape == (3, 1)
        nnt = np.dot(n_i[i], np.transpose(n_i[i]))

        assert nnt.shape == (3, 3)
        #A += nnt
        #It is correct if A contains equal rows. In this case, we have faces that are parallel or on the same plane (e.g. on the same side of a cube)
        assert p_i[i].shape == (3, 1)
        b += -np.dot(nnt, p_i[i] - x0)
        # IN PROGRESS

    return A, b



#@profile
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

@profile
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
        "sphere_example")  #
        # "bowl_15_holes")  # works too. But too many faces => too slow, too much memory. 32K?
    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-3, +5, 0.2)

    import vectorized, example_objects
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

    #exit()
    old_verts, old_facets = verts, facets

    # display_simple_using_mayavi_2( [(verts, facets),(verts, facets), ],
    #    pointcloud_list=[],
    #    mayavi_wireframe=[False, True,], opacity=[1, 1, 0.9], gradients_at=None, separate=False, gradients_from_iobj=None,
    #    minmax=(RANGE_MIN,RANGE_MAX)  )
    # exit()

    for i in range(VERTEX_RELAXATION_ITERATIONS_COUNT):
        verts, facets_not_used, centroids = process2_vertex_resampling_relaxation(verts, facets, iobj)
        print("Vertex relaxation applied.");sys.stdout.flush()


    total_subdivided_facets = []
    for i in range(SUBDIVISION_ITERATIONS_COUNT):
        e_array, bad_facets_count = compute_facets_subdivision_curvatures(verts, facets, iobj)

        e_array[np.isnan(e_array)] = 0  # treat NaN curvatures as zero curvature => no subdivision

        which_facets = np.arange(facets.shape[0])[ e_array > curvature_epsilon ]

        verts4_subdivided, facets3_subdivided = subdivide_multiple_facets(verts, facets, which_facets)
        global trace_subdivided_facets  # third implicit output
        verts, facets = verts4_subdivided, facets3_subdivided
        print("Subdivision applied.");sys.stdout.flush()

        total_subdivided_facets += trace_subdivided_facets  # old face indices remain valid

        for i in range(VERTEX_RELAXATION_ITERATIONS_COUNT):
            verts, facets_not_used, centroids = process2_vertex_resampling_relaxation(verts, facets, iobj)
            print("Vertex relaxation applied.");sys.stdout.flush()

    average_edge = compute_average_edge_length(verts, facets)

    c3 = np.mean(verts[facets[:], :], axis=1)
    old_centroids = np.concatenate((c3, np.ones((c3.shape[0], 1))), axis=1)

    new_centroids = old_centroids.copy()
    set_centers_on_surface_ohtake(iobj, new_centroids, average_edge)
    #new_centroids is the output

    # display_simple_using_mayavi_2( [(verts, facets),(verts, facets), ],
    #    pointcloud_list=[ new_centroids ], pointcloud_opacity=0.2,
    #    mayavi_wireframe=[False, True,], opacity=[1, 1, 0.9], gradients_at=None, separate=False, gradients_from_iobj=None,
    #    minmax=(RANGE_MIN,RANGE_MAX)  )
    # exit()



    # The two CHOICEs are equaivalent. Two rewrite of the same method.
    CHOICE = 1
    if CHOICE == 1:
        #neighbour_faces_of_vertex
        vertex_neighbours_list = mesh_utils.make_neighbour_faces_of_vertex(facets)
    #    import ipdb; ipdb.set_trace()
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

    highlighted_vertices = np.array([131,  71, 132])  # np.arange(100, 200)
    hv = new_verts_qem[highlighted_vertices, :]

    display_simple_using_mayavi_2( [(new_verts_qem_alpha, facets),(new_verts_qem, facets), ],
       pointcloud_list=[ hv ], pointcloud_opacity=0.2,
       mayavi_wireframe=[False,False], opacity=[0.4*0, 1, 0.9], gradients_at=None, separate=False, gradients_from_iobj=None,
       minmax=(RANGE_MIN,RANGE_MAX)  )

#from timeit import default_timer as dtimer


if __name__ == '__main__':

#    timer_t1s = dtimer()

    demo_combination_plus_qem()  # subdivision + projection + qem
