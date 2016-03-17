from vtk import *
import numpy as np
from vtk.util.colors import tomato
import os, sys


def make_npy_file_dice():
    dicesize = 16.
    exname = "udice_vec"  # "blend_example2"
    import example_objects
    iobj = example_objects.make_example_vectorized(exname, dicesize)

    from stl_tests import make_mc_values_grid
    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-22, +20., 0.8)
    numpy_array = make_mc_values_grid(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE, old=True)
    return numpy_array


def vtk_show_numpy_array(numpy_array, spacing=5):
    from vtk.util.numpy_support import numpy_to_vtk

    try:
        data_numpy = np.load(numpy_array)
    except:
        data_numpy = numpy_array

    x_dim, y_dim, z_dim = data_numpy.shape
    # Flatten the array so it can be loaded to VTK
    data_numpy = data_numpy.ravel(order='C')
    # keep a reference to the numpy data (not sure if needed)
    b = data_numpy
    #initialize vtkImageData instance
    vtk_Data = vtkImageData()
    vtk_Data.SetExtent( 0, x_dim - 1, 0, y_dim - 1 , 0, z_dim - 1)
    vtk_Data.SetSpacing(spacing, spacing, spacing) # check this
    vtk_Data.SetOrigin(-4, -4, -4)

    vtk_Data.SetNumberOfScalarComponents(1)
    vtk_Data.Update()
    # convert to vtk
    numpy_to_vtk_data = numpy_to_vtk(data_numpy)

    vtk_Data.GetPointData().SetScalars(numpy_to_vtk_data)


    #Init marching cubes and set outpout to mapper of class vtkPolyDataMapper()
    mcubes = vtkMarchingCubes()
    mcubes.SetInput(vtk_Data)
    mcubes.ComputeNormalsOn()
    mcubes.SetValue(0,0.0)
    mcubes_mapper = vtkPolyDataMapper()
    mcubes_mapper.SetInputConnection(mcubes.GetOutputPort())


    mcubes_actor = vtkActor()
    mcubes_actor.SetMapper(mcubes_mapper)
    mcubes_actor.GetProperty().SetOpacity(1)
    mcubes_actor.GetProperty().SetColor(tomato)
    renderer = vtkRenderer()
    camera = renderer.GetActiveCamera()
    camera.SetPosition(0,0,1000)
    renderer.AddActor(mcubes_actor)
    renderer.SetBackground(1,1,1)
    renWin = vtkRenderWindow()
    renWin.AddRenderer(renderer)
    renWin.SetSize(300, 300)

    iren = vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)
    iren.Initialize()
    renWin.Render()
    iren.Start()


if __name__ == "__main__":
    # vtk_show_numpy_array('vgrid_array.npy')
    vtk_show_numpy_array(make_npy_file_dice())
