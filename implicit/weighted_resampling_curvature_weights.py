import numpy as np
from optimize_dual_mesh import MeshOptimizer
from scipy.sparse import csr_matrix
from optimize_dual_mesh import build_neighbor_list_c2c
from ipdb import set_trace


def initialize():
    opt = MeshOptimizer()
    opt.load_example(res=0.5)
    opt.run(calc_proj=True, calc_opt=False, update_centroids=False)
    # opt.vertices = compute_weighted_resampling(opt)
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


def compute_weighted_resampling(opt, c=0.1):
    """ c is the user-specified constant """
    projections = opt.optimized_dual_mesh
    vertices = opt.vertices
    dim = len(projections)
    # print len(vertices)
    # exit()
    nbr_list = build_neighbor_list_c2c(opt.centroids, opt.faces, opt.vertex_neighbours_list)
    nbrs = projections[nbr_list]
    normals_p = build_normals(opt, projections).reshape(dim, 1, 3)
    normals_n = normals_p[nbr_list].reshape(dim, 3, 3)
    dott_prod = (normals_n * normals_p)     # element wise mult of normals_p, normals_n
    dott_prod = dott_prod.sum(axis=2)       # sum to get dot product
    arccos = np.arccos(dott_prod)           # vectorized arccos
    PiPj = np.linalg.norm(normals_n - normals_p, axis=2)    # compute norms of PiPj
    k_weights = (arccos / PiPj).sum(axis=1)
    w_weights = 1 + c * k_weights    # final w_weights
    return w_weights.reshape(1, dim)


def build_normals(opt, point_matrix):
    dim1, dim2 = point_matrix.shape
    assert not np.any(np.isnan(point_matrix)), "there should not be any NaN values"
    point_matrix = point_matrix.reshape(dim1, 4)
    n = (opt.gradient(point_matrix))[:, :3].reshape(dim1, 3)
    n /= np.linalg.norm(n, axis=1).reshape(dim1, 1)
    return n

if __name__ == "__main__":
    opt = initialize()
    weights = compute_weighted_resampling(opt)
    # w = np.repeat(weights, 456, axis=0)
    sparse = build_csr_matrix(opt.vertices, opt.vertex_neighbours_list).toarray()
    median_norm = np.median(weights)    # will be used form normalisation of outliers
    sparse[np.where(sparse > 100)] = median_norm
    weights[np.where(weights > 100)] = median_norm
    sparse_of_weights = (sparse * weights)
    norm_factor = sparse_of_weights.sum(axis=1)
    new_verts = np.dot(sparse, opt.optimized_dual_mesh)
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

    mlab.show()
