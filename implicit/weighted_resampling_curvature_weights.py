import numpy as np
from optimize_dual_mesh import MeshOptimizer
from scipy.sparse import csr_matrix
from optimize_dual_mesh import build_neighbor_list_c2c
from ipdb import set_trace


def initialize():
    opt = MeshOptimizer()
    opt.load_example(res=0.5)
    opt.run(calc_proj=True, calc_opt=False, update_centroids=False)
    opt.vertices = compute_weighted_resampling(opt)
    return opt


def compute_weighted_resampling(opt):
    projections = opt.optimized_dual_mesh
    vertices = opt.vertices
    nbr_list = build_neighbor_list_c2c(opt.centroids, opt.faces, opt.vertex_neighbours_list)
    nbrs = projections[nbr_list]
    normals_p = build_normals(opt, projections).reshape(len(opt.centroids), 3, 1)
    normals_n = normals_p[nbr_list].reshape(len(opt.centroids), 3, 3)
    dott_prod = (normals_n * normals_p)
    dott_prod = dott_prod.sum(axis=2)
    set_trace()

    
def build_normals(opt, point_matrix):
    dim1, dim2 = point_matrix.shape
    assert not np.any(np.isnan(point_matrix)), "there should not be any NaN values"
    point_matrix = point_matrix.reshape(dim1, 4)
    n = (opt.gradient(point_matrix))[:, :3].reshape(dim1, 3)
    n /= np.linalg.norm(n, axis=1).reshape(dim1, 1)
    return n

if __name__ == "__main__":
    opt = initialize()
    compute_weighted_resampling(opt)
