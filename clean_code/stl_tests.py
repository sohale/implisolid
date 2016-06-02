import numpy as np
#from stl import mesh


"""Creates and exports STL files based on polygonized implicit objects.
STL files can be used directly in a slicer like Cura for printing.
To install the packages:
pip install numpy-stl
"""


def write_stl(verts, faces, filename):
    pass


def load_stl(stl_fn, save_stl_fn):

    # Using an existing stl file:
    your_mesh = mesh.Mesh.from_file(stl_fn)

    # Or creating a new mesh (make sure not to overwrite the `mesh` import by
    # naming it `mesh`):
    VERTICE_COUNT = 100
    data = np.zeros(VERTICE_COUNT, dtype=mesh.Mesh.dtype)
    your_mesh = mesh.Mesh(data, remove_empty_areas=False)

    # The mesh normals (calculated automatically)
    your_mesh.normals
    # The mesh vectors
    your_mesh.v0, your_mesh.v1, your_mesh.v2
    # Accessing individual points (concatenation of v0, v1 and v2 in triplets)
    assert (your_mesh.points[0][0:3] == your_mesh.v0[0]).all()
    assert (your_mesh.points[0][3:6] == your_mesh.v1[0]).all()
    assert (your_mesh.points[0][6:9] == your_mesh.v2[0]).all()
    assert (your_mesh.points[1][0:3] == your_mesh.v0[1]).all()

    your_mesh.save(save_stl_fn)

#
# def show_stl(stl_fn):
#     from stl import mesh
#     from mpl_toolkits import mplot3d
#     from matplotlib import pyplot
#
#     # Create a new plot
#     figure = pyplot.figure()
#     axes = mplot3d.Axes3D(figure)
#
#     # Load the STL files and add the vectors to the plot
#     your_mesh = mesh.Mesh.from_file(stl_fn)
#     axes.add_collection3d(mplot3d.art3d.Poly3DCollection(your_mesh.vectors))
#
#     # Auto scale to the mesh size
#     scale = your_mesh.points.flatten(-1)
#     axes.auto_scale_xyz(scale, scale, scale)
#
#     # Show the plot to the screen
#     pyplot.show()

# load_stl('some_file.stl', 'new_stl_file.stl')
# show_stl('tests/stl_binary/HalfDonut.stl')


def plot_stlmesh(m):
    from matplotlib import pyplot
    from mpl_toolkits import mplot3d

    # Create a new plot
    figure = pyplot.figure()
    axes = mplot3d.Axes3D(figure)

    # Render the cube faces
    # for m in meshes:
    axes.add_collection3d(mplot3d.art3d.Poly3DCollection(m.vectors))

    # Auto scale to the mesh size
    scale = np.concatenate([m]).flatten(-1)
    axes.auto_scale_xyz(scale, scale, scale)

    pyplot.show()


def display_simple_using_mayavi_vf1(verts, faces, minmax=(-1, 1), mayavi_wireframe=False):
    from mayavi import mlab
    mlab.triangular_mesh([vert[0] for vert in verts],
                         [vert[1] for vert in verts],
                         [vert[2] for vert in verts], faces, representation="surface" if not mayavi_wireframe else "wireframe", opacity=1, scale_factor=100.0)

    (RANGE_MIN, RANGE_MAX) = minmax
    x = np.linspace(RANGE_MIN, RANGE_MAX, 2).reshape(2, 1)
    y = np.zeros((2, 1))
    z = np.zeros((2, 1))

    mlab.plot3d(x, y, z, line_width=3, name="x-axis")
    mlab.plot3d(y, x, z, line_width=3, name="y-axis")
    mlab.plot3d(z, y, x, line_width=3, name="z-axis")

    mlab.show()     # figure=fig,
#


def make_mc_mesh_scikit(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE):
    """ Uses Scikit's MC algorithm,which has minor bugs. """
    rng = np.arange(RANGE_MIN, RANGE_MAX, STEPSIZE)
    import mc_utils
    vgrid = mc_utils.make_mc_values_grid(iobj, rng, old=True)
    from skimage import measure
    verts, faces = measure.marching_cubes(vgrid, 0)
    verts = ((verts) * STEPSIZE + rng[0])
    print("OLD: swapping x,y")
    verts = np.concatenate((verts[:, 1, np.newaxis], verts[:, 0, np.newaxis], verts[:, 2, np.newaxis]), axis=1)
    return verts, faces

# @profile

# todo: use "3" or avoid using old=True


def make_mc_values_grid(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE, old=True):
    rng = np.arange(RANGE_MIN, RANGE_MAX, STEPSIZE)
    import mc_utils
    vgrid = mc_utils.make_grid(iobj, rng, old=old)
    if old:
        return np.swapaxes(vgrid, 0, 1)
    else:
        print ("*********************************************")
        vgrid = np.swapaxes(vgrid, 1, 2)
        vgrid = np.swapaxes(vgrid, 0, 1)
        return vgrid
        # print("no swap")
        # return vgrid

#
# def make_mc_values_grid_mayavi(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE):
#     rng = np.arange(RANGE_MIN, RANGE_MAX, STEPSIZE)
#     import mc_utils
#     vgrid = mc_utils.make_grid(iobj, rng, old=True)
#     return np.swapaxes(vgrid, 0, 1)


def test3():
    exname = "screw3"
    import example_objects
    iobj = example_objects.make_example_vectorized(exname)

    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-2.5, +2.5, 0.1)
    verts, faces = make_mc_mesh_scikit(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE)
    display_simple_using_mayavi_vf1(verts, faces)


# test3()

def m2stl_mesh(verts, faces):
    from stl import mesh
    fv = verts[faces, :]
    print fv.shape

    data = np.zeros(fv.shape[0], dtype=mesh.Mesh.dtype)
    for i in range(fv.shape[0]):
        facet = fv[i]
        data['vectors'][i] = facet

    m = mesh.Mesh(data)
    return m


def test4():
    # simply marching cubes

    exname = "screw3"
    import example_objects
    iobj = example_objects.make_example_vectorized(exname)

    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-2.5, +2.5, 0.1)
    verts, faces = make_mc_mesh_scikit(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE)

    m = m2stl_mesh(verts, faces)
    m.save('m.stl')

    plot_stlmesh(m)

    # display_simple_using_mayavi_vf1(verts, faces)


def test5_screw():
    # simply marching cubes

    exname = "screw3"
    import example_objects
    iobj = example_objects.make_example_vectorized(exname, 8.0)

    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-2.5 * 8, +2.5 * 8, 0.1 * 8)
    verts, faces = make_mc_mesh_scikit(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE)

    # verts = optimise_mesh(verts, faces, iobj)

    m = m2stl_mesh(verts, faces)
    if ACTUALLY_SAVE:
        m.save('m-optim2.stl')

    display_simple_using_mayavi_vf1(verts, faces)

    plot_stlmesh(m)


def test6_blend():
    """Printed nicely."""

    exname = "blend_example2_discs"  # "blend_example2"
    import example_objects
    iobj = example_objects.make_example_vectorized(exname, 8.0)

    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-2. * 8, +4. * 8, 0.4 * 8 / 5)
    verts, faces = make_mc_mesh_scikit(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE)

    # verts = optimise_mesh(verts, faces, iobj)

    m = m2stl_mesh(verts, faces)
    if ACTUALLY_SAVE:
        m.save('implicit6-blend.stl')  # wow

    display_simple_using_mayavi_vf1(verts, faces)

    plot_stlmesh(m)


def test7_dice():
    """ Dice prints well. I used NEtfabb to correct the STL though. """
    """ It is interesting that the result DOES depend on size (Scaling
    everything inslucing the grid step size). It works well when scale is x1
    and works perfect when scale is x2. But if x3, it starts to look rubbish.
    Because the numerical methods use absolute sizes (lambda, etc?).
    """
    # -8.8, 7.2
    rescale = 1.
    dicesize = rescale * 8.
    exname = "udice_vec"  # "blend_example2"
    import example_objects
    iobj = example_objects.make_example_vectorized(exname, dicesize)

    # (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-2.*8, +4.*8, 0.4*8/4)
    # (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-16., +32., 0.8)  #non-spiky
    # (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-10, +9., 0.8)  # spikes at bottom!
    (RANGE_MIN, RANGE_MAX, STEPSIZE) = \
        (-11 * rescale, +10. * rescale, 0.8 * rescale)  # non-spiky!!
    # (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-8, +8., 0.8)  #
    verts, faces = make_mc_mesh_scikit(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE)
    from numerical_utils import average_edge_size
    print "average_edge_size", average_edge_size(verts, faces)
    # When rescale is doubled (=2.), average edge size 0.824 -> 1.647 .
    # Seems the latter is when it works very well.

    print "average_edge_size", average_edge_size(verts, faces)  # 0.86 -> 1.72

    m = m2stl_mesh(verts, faces)
    if ACTUALLY_SAVE:
        m.save('stl/implicit7-dice.stl')  # wow

    display_simple_using_mayavi_vf1(verts, faces)
    print(np.min(verts.ravel()), np.max(verts.ravel()))

    plot_stlmesh(m)


def test8_bigdice():
    """ A larger dice. May need support though.
    Neat, high resolution, a bit heavy. """
    # -8.8, 7.2
    dicesize = 16.
    exname = "udice_vec"  # "blend_example2"
    import example_objects
    iobj = example_objects.make_example_vectorized(exname, dicesize)

    # (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-22, +20., 0.8/2)  # non-spiky!!
    # will look perfect. Neat, high resolution, a bit heavy.
    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-22, +20., 0.8)
    verts, faces = make_mc_mesh_scikit(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE)

    # verts = optimise_mesh(verts, faces, iobj)

    m = m2stl_mesh(verts, faces)
    if ACTUALLY_SAVE:
        m.save('stl/implicit8-bigdice-.stl')  # wow

    display_simple_using_mayavi_vf1(verts, faces)
    print(np.min(verts.ravel()), np.max(verts.ravel()))

    plot_stlmesh(m)


def test9_icesl1():
    """ Comparing with IceSL. """

    import vectorized
    ns = vectorized
    #    #40,25
    #    #s1 = 25/10.
    #    r1 = 40/10.
    #    m = np.eye(4) * r1
    #    m[3, 3] = 1
    #    m[0,3] = 25  #tooth, r1=4.
    #    print(m)

    sc = 1. / 4.
    r1 = 2 * 25 * sc
    m = np.eye(4) * r1
    m[3, 3] = 1
    print(m)

    a = ns.Ellipsoid(m)
    b = ns.UnitCube1(2 * 40. * sc)
    iobj = ns.CrispSubtract(b, a)
    # iobj = ns.CrispUnion(b, a)
    # iobj = ns.Transformed(d)

    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-90 * sc, +90. * sc, 4. * sc / 2.)
    verts, faces = make_mc_mesh_scikit(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE)

    # verts = optimise_mesh(verts, faces, iobj)

    m = m2stl_mesh(verts, faces)
    if ACTUALLY_SAVE:
        m.save('stl/icesl_example1.stl')  # wow

    display_simple_using_mayavi_vf1(verts, faces)
    print(np.min(verts.ravel()), np.max(verts.ravel()))

    plot_stlmesh(m)


if __name__ == '__main__':

    # Slice using: Cura or #https://sketchfab.com/

    ACTUALLY_SAVE = False

    # test4()  #simple screw with raw MC
    # test5_screw()
    # test6_blend()
    # test7_dice()  # +
    # test8_bigdice()
    test9_icesl1()
