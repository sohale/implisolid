from .ohtake_belyaev_demo_subdivision_projection_qem import *
import json

mesh_correction = False
writing_test_file = False
B = 1000000L


def puppy_magic(mp5source):

    from mp5toufunc import get_fonuky, get_root_node
    root_node = get_root_node(mp5source)
    iobj = get_fonuky(root_node)

    from stl_tests import make_mc_values_grid
    gridvals = make_mc_values_grid(iobj, -60., 60., 1., old=False)
    vertex, faces = vtk_mc(gridvals, (-60., 60., 1.))
    sys.stdout.flush()



    return m2stl_mesh(vertex, faces)

    # display_simple_using_mayavi_2([(vertex_before_qem, faces), (new_vertex_qem, faces), ],
    #    pointcloud_list=[],
    #    mayavi_wireframe=[False, False], opacity=[0.4*0, 1, 0.9], gradients_at=None, separate=False, gradients_from_iobj=None,
    #    minmax=(RANGE_MIN, RANGE_MAX))


def m2stl_mesh(verts, faces):
    from stl import mesh
    fv = verts[faces, :]

    data = np.zeros(fv.shape[0], dtype=mesh.Mesh.dtype)
    for i in range(fv.shape[0]):
        facet = fv[i]
        data['vectors'][i] = facet

    m = mesh.Mesh(data)
    return m

