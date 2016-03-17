
import numpy as np
#from ipdb import set_trace
#see http://nullege.com/codes/show/src@p@y@pyformex-0.9.0@pyformex@plugins@vtk_itf.py/36/vtk.util.numpy_support.vtk_to_numpy

# def cleanVPD(vpd):
#     """Clean the vtkPolydata

#     Clean the vtkPolydata, adjusting connectivity, removing duplicate elements
#     and coords, renumbering the connectivity. This is often needed after
#     setting the vtkPolydata, to make the vtkPolydata fit for use with other
#     operations. Be aware that this operation will change the order and
#     numbering of the original data.

#     Parameters:

#     - `vpd`: a vtkPolydata

#     Returns the cleaned vtkPolydata.
#     """
#     from vtk import vtkCleanPolyData
#     cleaner=vtkCleanPolyData()
#     cleaner.SetInput(vpd)
#     cleaner.Update()
#     return cleaner.GetOutput()


# def convert2VPD(M,clean=False,verbose=False):
#     """Convert pyFormex data to vtkPolyData.

#     Convert a pyFormex Mesh or Coords object into vtkPolyData.
#     This is limited to vertices, lines, and polygons.
#     Lines should already be ordered (with connectedLineElems for instance).

#     Parameters:

#     - `M`: a Mesh or Coords type. If M is a Coords type it will be saved as
#       VERTS. Else...
#     - `clean`: if True, the resulting vtkdata will be cleaned by calling
#       cleanVPD.

#     Returns a vtkPolyData.
#     """
#     from vtk import vtkPolyData,vtkPoints,vtkIdTypeArray,vtkCellArray

#     if verbose:
#         print('STARTING CONVERSION FOR DATA OF TYPE %s '%type(M))

#     if  False: # isinstance(M,Coords):
#         M = Mesh(M,arange(M.ncoords()))

#     Nelems = M.nelems() # Number of elements
#     Ncxel = M.nplex() # # Number of nodes per element

#     # create a vtkPolyData variable
#     vpd=vtkPolyData()

#     # creating  vtk coords
#     pts = vtkPoints()
#     ntype=gnat(pts.GetDataType())
#     coordsv = n2v(asarray(M.coords,order='C',dtype=ntype),deep=1) #.copy() # deepcopy array conversion for C like array of vtk, it is necessary to avoid memry data loss
#     pts.SetNumberOfPoints(M.ncoords())
#     pts.SetData(coordsv)
#     vpd.SetPoints(pts)


#     # create vtk connectivity
#     elms = vtkIdTypeArray()
#     ntype=gnat(vtkIdTypeArray().GetDataType())
#     elmsv = concatenate([Ncxel*ones(Nelems).reshape(-1,1),M.elems],axis=1)
#     elmsv = n2v(asarray(elmsv,order='C',dtype=ntype),deep=1) #.copy() # deepcopy array conversion for C like array of vtk, it is necessary to avoid memry data loss
#     elms.DeepCopy(elmsv)

#     # set vtk Cell data
#     datav = vtkCellArray()
#     datav.SetCells(Nelems,elms)
#     if Ncxel == 1:
#         try:
#             if verbose:
#                 print("setting VERTS for data with %s maximum number of point for cell "%Ncxel)
#             vpd.SetVerts(datav)
#         except:
#             raise ValueError,"Error in saving  VERTS"

#     elif Ncxel == 2:
#         try:
#             if verbose:
#                 print ("setting LINES for data with %s maximum number of point for cell "%Ncxel)
#             vpd.SetLines(datav)
#         except:
#             raise  ValueError,"Error in saving  LINES"

#     else:
#         try:
#             if verbose:
#                 print ("setting POLYS for data with %s maximum number of point for cell "%Ncxel)
#             vpd.SetPolys(datav)
#         except:
#             raise ValueError,"Error in saving  POLYS"

#     vpd.Update()
#     if clean:
#         vpd=cleanVPD(vpd)
#     return vpd


# def convertFromVPD(vpd, verbose=False):
#     """Convert a vtkPolyData into pyFormex objects.
#     Convert a vtkPolyData into pyFormex objects.
#     Parameters:
#     - `vpd`: a vtkPolyData

#     Returns a tuple with points, polygons, lines, vertices numpy arrays.
#     Returns None for the missing data.

#     """
#     pts = polys = lines = verts = None

#     # getting points coords
#     if vpd.GetPoints().GetData().GetNumberOfTuples():
#         ntype = gnat(vpd.GetPoints().GetDataType())
#         pts = asarray(v2n(vpd.GetPoints().GetData()),dtype=ntype)
#         if verbose:
#             print('Saved points coordinates array')

#     # getting Polygons
#     if vpd.GetPolys().GetData().GetNumberOfTuples():
#         ntype=gnat(vpd.GetPolys().GetData().GetDataType())
#         Nplex = vpd.GetPolys().GetMaxCellSize()
#         polys = asarray(v2n(vpd.GetPolys().GetData()),dtype=ntype).reshape(-1,Nplex+1)[:,1:]
#         if verbose:
#             print('Saved polys connectivity array')

#     # getting Lines
#     if vpd.GetLines().GetData().GetNumberOfTuples():
#         ntype=gnat(vpd.GetLines().GetData().GetDataType())
#         Nplex = vpd.GetLines().GetMaxCellSize()
#         lines = asarray(v2n(vpd.GetLines().GetData()),dtype=ntype).reshape(-1,Nplex+1)[:,1:]
#         if verbose:
#             print('Saved lines connectivity array')

#     # getting Vertices
#     if vpd.GetVerts().GetData().GetNumberOfTuples():
#         ntype=gnat(vpd.GetVerts().GetData().GetDataType())
#         Nplex = vpd.GetVerts().GetMaxCellSize()
#         verts = asarray(v2n(vpd.GetVerts().GetData()),dtype=ntype).reshape(-1,Nplex+1)[:,1:]
#         if verbose:
#             print('Saved verts connectivity array')

#     return pts, polys, lines, verts


def vtk_mc(gridvals, rrr):
    (RANGE_MIN, RANGE_MAX, STEPSIZE) = rrr

    from vtk import *

    data_numpy = gridvals

    #np.save("dice.npy", data_numpy)

    #spacing = STEPSIZE  # 1. #5  # STEPSIZE
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

    fcount = mesh_data.GetNumberOfCells()
    fa = np.zeros((fcount, 3), dtype=int)
    for i in range(fcount):
        c = mesh_data.GetCell(i)
        (fa[i, 0], fa[i, 1], fa[i, 2]) = (c.GetPointId(0), c.GetPointId(1), c.GetPointId(2))

    # contour.Delete()

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

    #set_trace()
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

    #gridvals = make_npy_file_dice()
    verts, faces = vtk_mc(gridvals, (RANGE_MIN, RANGE_MAX, STEPSIZE) )
    print ("MC calculated")
    sys.stdout.flush()

    from mesh_utils import mesh_invariant
    mesh_invariant(faces)

    #from stl_tests import make_mc_mesh_scikit
    #verts, faces = make_mc_mesh_scikit(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE)

    display_simple_using_mayavi_vf1(verts, faces)

if __name__ == '__main__':
    #cyl_test_example1()
    #vtk_mc()
    vtk_mc_test()
