from vtk import *
import numpy as np
from vtk.util.numpy_support import numpy_to_vtk
from vtk.util.colors import tomato

import os, sys
ROOT_DIR = os.getcwd()
os.listdir(ROOT_DIR)
# Flatten the array so it can be loaded to VTK

data_numpy = np.load('vgrid_array.npy')
data_numpy_shape = data_numpy.shape
data_numpy = data_numpy.ravel(order='C')
b = data_numpy
# initialize vtkImage.Data instance

vtk_Data = vtkImageData()
vtk_Data.SetExtent(0,26,0,26,0,26)
vtk_Data.SetSpacing(10,10,10) # check this
vtk_Data.SetOrigin(0,0,0)

vtk_Data.SetNumberOfScalarComponents(1)
vtk_Data.Update()
# convert to vtk
numpy_to_vtk_data = numpy_to_vtk(data_numpy)

vtk_Data.GetPointData().SetScalars(numpy_to_vtk_data)

# voi = vtkExtractVOI()
# voi.SetInputConnection(vtk_Data.GetO)
#Init marching cubes and set outpout to mapper of class vtkPolyDataMapper()
mcubes = vtkMarchingCubes()
mcubes.SetInput(vtk_Data)
mcubes.SetValue(0,0.0)
mcubes_mapper = vtkPolyDataMapper()
mcubes_mapper.SetInputConnection(mcubes.GetOutputPort())



mcubes_actor = vtkActor()
mcubes_actor.SetMapper(mcubes_mapper)
mcubes_actor.GetProperty().SetOpacity(1)
mcubes_actor.GetProperty().SetColor((1,1,1))
mcubes_actor.GetProperty.SetRepresentation('wireframe')

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
