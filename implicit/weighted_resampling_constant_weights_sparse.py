from scipy.sparse import csr_matrix
import numpy as np
from optimize_dual_mesh import MeshOptimizer
from ipdb import set_trace
from mayavi import mlab


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


def main(plots=False):
    res = 1
    opt = MeshOptimizer()
    opt.load_example(res=res)
    opt.run(calc_opt=True, update_centroids=False)
    sparse = build_csr_matrix(opt.vertices, opt.vertex_neighbours_list)
    row_sizes = sparse.getnnz(axis=1)
    row_sizes = row_sizes.reshape(len(row_sizes), 1)
    # set_trace()
    sparse_array = sparse.toarray()
    # dense_array = sparse.todense()
    #
    new_verts = sparse_array.dot(opt.optimized_dual_mesh) / row_sizes

    # set plot to True to enable visualisation of the resampling results.

    if plots:
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
    # new_verts_dense = np.dot(dense_array, opt.optimized_dual_mesh)
    return new_verts
if __name__ == "__main__":
    main()
