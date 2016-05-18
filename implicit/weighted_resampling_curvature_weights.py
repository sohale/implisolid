import numpy as np
from optimize_dual_mesh import MeshOptimizer
from scipy.sparse import csr_matrix
from optimize_dual_mesh import build_neighbor_list_c2c
from ipdb import set_trace


def initialize():
    opt = MeshOptimizer()
    opt.load_example(res=2.5)
    opt.run(calc_proj=True, calc_opt=False, update_centroids=False)
    # opt.vertices = compute_weighted_resampling(opt)
    print "number of vertices,number of centroids:, " + repr(len(opt.vertices)) + repr(len(opt.centroids))
    return opt


def build_csr_matrix(vertices, neighbors):
    """
    This function builds the CSR sparse matrix of vertices and neighboring_projections
    If Nv is the number of vertices and Np  is the number of projections(or trianges or faces or centroids)
    the sparse matrix will be of size Nv x Np, and its entry in position i,j will be 1 if the i-th vertex
    is neighbor to the j-th projection
    """
    indptr = [0]
    indices = []
    data = []
    index = 0
    # incremental construction of matrix
    # see documentation of scipy.sparse.csr_matrix
    for i in range(len(vertices)):
        lgth = len(neighbors[i])
        indices += neighbors[i]
        index += lgth
        indptr.append(index)
        data += (lgth * [1])
    sparse_mat = csr_matrix((data, indices, indptr))
    return sparse_mat

@profile
def compute_weighted_resampling(opt, sparse_compressed, c=10, num_iters=10):
    """
        This function computes the weighted resampling process for a mesh specified by opt and sparse_compressed


        ARGUMENTS

        opt              : MeshOptimizer object that holds information about the mesh, such as
                           vertices, faces, triangles, implicitFunction, implicitGradient (almost similar to iobj)
        sparse_compressed: is the sparse matrix in CSR format of faces neighbors to vertices
        c                : is the user-specified constant
        num_iters        : the number of iterations for the weighted resampling

        RETURNS

        The function does not return a value, instead it updates the mesh vertices in place.
    """
    sparse = sparse_compressed.toarray()

    nbr_list = build_neighbor_list_c2c(opt.centroids, opt.faces, opt.vertex_neighbours_list)
    # set_trace()
    for i in range(num_iters):
        if i > 0:
            opt.run(calc_proj=True, calc_opt=False, update_centroids=True)
        projections = opt.optimized_dual_mesh
        vertices = opt.vertices
        dim = len(projections)
        # print len(vertices)
        # exit()

        nbrs = projections[nbr_list]           #  faces of faces neighbors
        normals_p = build_normals(opt, projections)  #  normals at projections (calculated once)
        normals_p = normals_p.reshape(dim, 1, 3)
        # set_trace()
        normals_n = normals_p[nbr_list].reshape(dim, 3, 3)  # normals at each neighbor projection (fancy indexing to avoid excessive calculations)
        dott_prod = (normals_n * normals_p)     #
        dott_prod = dott_prod.sum(axis=2)       # m(Pi) * m(Pj)
        arccos = np.arccos(dott_prod)           # arccos(m(Pi) * m(Pj))
        PiPj = np.linalg.norm(projections[:,np.newaxis,:3] - nbrs[:, :, :-1], axis=2)    # ||Pi - Pj||
        k_weights = (arccos / PiPj).sum(axis=1)     # ki = Sum_k ( arccos(m(Pi) * m(Pj))/||Pi - Pj|| )
        weights = 1 + c * k_weights                 # wi = 1 + c*ki
        weights.shape = (1, dim)    # final w_weights
        # median_norm = np.median(weights)    # will be used form normalisation of outliers  TODO: might not be needed after correcting line 60
        # sparse[np.where( sparse > 1000)] = median_norm
        # weights[np.where(weights > 1000)] = median_norm
        sparse_of_weights = (sparse * weights)
        norm_factor = sparse_of_weights.sum(axis=1)
        # set_trace()
        new_verts = np.dot(sparse_of_weights, projections)
        new_verts /= norm_factor.reshape(len(opt.vertices),1)
        from mayavi import mlab

        mlab.figure()
        mlab.triangular_mesh([vert[0] for vert in opt.vertices],
                         [vert[1] for vert in opt.vertices],
                         [vert[2] for vert in opt.vertices],
                         opt.faces, representation="surface")

        mlab.figure()
        mlab.triangular_mesh([vert[0] for vert in new_verts],
                         [vert[1] for vert in new_verts],
                         [vert[2] for vert in new_verts],
                         opt.faces, representation="surface")
    #-note
    #
    #     # set_trace()

        opt.vertices = new_verts
    mlab.show()


@profile
def build_normals(opt, point_matrix):
    dim1, dim2 = point_matrix.shape
    assert not np.any(np.isnan(point_matrix)), "there should not be any NaN values"
    point_matrix = point_matrix.reshape(dim1, 4)
np    n = (opt.gradient(point_matrix))[:, :3].reshape(dim1, 3)
    n /= np.linalg.norm(n, axis=1).reshape(dim1, 1)
    return n

if __name__ == "__main__":
    opt = initialize()
    # w = np.repeat(weights, 456, axis=0)

    sparse = build_csr_matrix(opt.vertices, opt.vertex_neighbours_list)
    compute_weighted_resampling(opt, sparse, c=2, num_iters=1)
