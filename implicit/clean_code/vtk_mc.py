
import numpy as np
# from ipdb import set_trace
# see http://nullege.com/codes/show/src@p@y@pyformex-0.9.0@pyformex@plugins@vtk_itf.py/36/vtk.util.numpy_support.vtk_to_numpy

def vtk_mc(gridvals, rrr):
    (RANGE_MIN, RANGE_MAX, STEPSIZE) = rrr

    from vtk import *
    
    data_numpy = gridvals

    # np.save("dice.npy", data_numpy)

    # spacing = STEPSIZE  # 1. #5  # STEPSIZE
    spacing = 1

    x_dim, y_dim, z_dim = data_numpy.shape
    # Flatten the array so it can be loaded to VTK
    data_numpy = data_numpy.ravel(order='C')
    # convert to vtk
    from vtk.util.numpy_support import numpy_to_vtk
    numpy_to_vtk_data = numpy_to_vtk(data_numpy)

    vtk_image_grid = vtkImageData()  # sample , vtk_Data
    vtk_image_grid.SetExtent(0, x_dim-1, 0, y_dim-1, 0, z_dim-1)
    vtk_image_grid.SetSpacing(spacing, spacing, spacing)  # check this
    vtk_image_grid.SetOrigin(0, 0, 0)
    vtk_image_grid.SetNumberOfScalarComponents(1)
    vtk_image_grid.Update()

    vtk_image_grid.GetPointData().SetScalars(numpy_to_vtk_data)

    # Using vtkMarchingCubes class:
    #    #Init marching cubes and set outpout to mapper of class vtkPolyDataMapper()
    #    mcubes = vtkMarchingCubes()
    #    mcubes.SetInput(vtk_image_grid)
    #    mcubes.ComputeNormalsOn()  # ?
    #    mcubes.SetValue(0, 0.0)
    #    #o = mcubes.GetOutputPort()
    #    #source = o
    #    #source=convert2VPD(source)
    #
    #    #mcubes.Update()  # sample

    # contour
    contour = vtkContourFilter()
    contour.SetInput(vtk_image_grid)
    contour.GenerateValues(1, 0, 1)  # (int numContours, double rangeStart, double rangeEnd)
    contour.Update()

    mesh_data = contour.GetOutput()

    vcount = mesh_data.GetNumberOfPoints()
    va = np.zeros((vcount, 3))
    for i in range(vcount):
        p = mesh_data.GetPoint(i)
        va[i, :] = p
        va[i, :] = va[i, :] * STEPSIZE + RANGE_MIN

    # va = va[:, [0,1,2]]

    fcount = mesh_data.GetNumberOfCells()
    fa = np.zeros((fcount, 3), dtype=int)
    for i in range(fcount):
        c = mesh_data.GetCell(i)
        (fa[i, 0], fa[i, 1], fa[i, 2]) = (c.GetPointId(0), c.GetPointId(1), c.GetPointId(2))

    # contour.Delete()

    va = va[:, [1, 0, 2]]

    return va, fa

    """
    #mesh = MeshData()
    print("OK")
    return None

    mcubes_mapper = vtkPolyDataMapper()
    mcubes_mapper.SetInputConnection(o)

    mcubes_mapper.SetInput(mcubes)
    # Then what to do next with mcubes?

    #mcubes_mapper = vtkPolyDataMapper()
    #mcubes_mapper.SetInputConnection(mcubes.GetOutputPort())
    #mcubes_mapper

    from vtk.util.numpy_support import vtk_to_numpy as v2n

    #Coords(convertFromVPD(transformFilter.GetOutput())[0])
    a = convertFromVPD( source )
    print a
    """

from stl_tests import display_simple_using_mayavi_vf1
import sys


def vtk_mc_test():

    # set_trace()
    dicesize = 16.
    exname = "udice_vec"  # "blend_example2"
    import example_objects
    iobj = example_objects.make_example_vectorized(exname, dicesize)
    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-22, +20., 0.8)

    from example_objects import cyl4
    iobj, (RANGE_MIN, RANGE_MAX, STEPSIZE) = cyl4()

    from stl_tests import make_mc_values_grid
    numpy_array = make_mc_values_grid(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE, old=False)
    gridvals = numpy_array

    # gridvals = make_npy_file_dice()
    verts, faces = vtk_mc(gridvals, (RANGE_MIN, RANGE_MAX, STEPSIZE))
    print ("MC calculated")
    sys.stdout.flush()

    from mesh_utils import mesh_invariant
    mesh_invariant(faces)

    # from stl_tests import make_mc_mesh_scikit
    # verts, faces = make_mc_mesh_scikit(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE)

    display_simple_using_mayavi_vf1(verts, faces)

if __name__ == '__main__':
    # cyl_test_example1()
    # vtk_mc()
    vtk_mc_test()
