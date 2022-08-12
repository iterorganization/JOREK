from __future__ import print_function
import vtk
import numpy as np
import jorek.jorek_read_h5 as jorek

f = jorek.fields()
f.read('/tmp/jorek_restart.h5', variables=[2, 5])


grid = f.to_vtk(n_sub=3, phi=[0,180], n_plane=3, bezier=False)

#writer = vtk.vtkXMLUnstructuredGridWriter()
writer = vtk.vtkUnstructuredGridWriter()
#writer.SetDataModeToAscii()
writer.SetFileName('jorek_3d_quadratic.vtk')
writer.SetInputData(grid)
writer.Write()

if False:  # Visualize
    mapper = vtk.vtkDataSetMapper()
    if vtk.VTK_MAJOR_VERSION <= 5:
        mapper.SetInput(grid)
    else:
        mapper.SetInputData(grid)
    # Create color map based on range of first var
    data_range = grid.GetPointData().GetAbstractArray(0).GetRange()
    colormap = vtk.vtkLookupTable()
    colormap.SetTableRange(data_range[0], data_range[1])
    colormap.Build()
    mapper.SetColorModeToDefault()
    mapper.SetScalarRange(data_range)
    mapper.SetScalarVisibility(True)
    mapper.SetLookupTable(colormap)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetPointSize(20)

    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)

    renderer.AddActor(actor)

    renderWindow.Render()
    renderWindowInteractor.Start()
