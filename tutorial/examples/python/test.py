import sys
sys.path.append("../../..")
import plantbox as pb
import numpy as np
import vtk
from vtk.util.colors import *

# VTK right approach is vtkLineSource vtkTubeFilter


def test(rs):
    ren = vtk.vtkRenderer()
    # ren.SetBackground(1., 1., 1.)

    poly = rs.getPolylines()
    a_ = rs.getParameter("radius")

    for i, pl in enumerate(poly):

        col = forest_green  # tomato
        a = a_[i]

        polyData = vtk.vtkPolyData();
        points = vtk.vtkPoints()
        lines = vtk.vtkCellArray()
        # polyData.InsertNextCell(len(pl))
        for j, p in enumerate(pl):
            points.InsertNextPoint(p.x, p.y, p.z)
            lines.InsertCellPoint(j)

        polyData.SetPoints(points)
        polyData.SetLines(lines)

        tubeFilter = vtk.vtkTubeFilter()
        tubeFilter.SetRadius(a)
        tubeFilter.SetInputData(polyData)

        polyMapper = vtk.vtkPolyDataMapper()
        # polyMapper.SetInputData(polyData)
        polyMapper.SetInputConnection(tubeFilter.GetOutputPort())

        cylinderActor = vtk.vtkActor()
        cylinderActor.SetMapper(polyMapper)
        cylinderActor.GetProperty().SetColor(col)
#             o = y - x
#             cylinderActor.SetOrientation(o[0], o[1], o[2])

        ren.AddActor(cylinderActor)

    return ren


if __name__ == '__main__':

    rs = pb.RootSystem()
    path = "../../../modelparameter/rootsystem/"
    name = "Anagallis_femina_Leitner_2010"
    rs.readParameters(path + name + ".xml")
    rs.initialize()
    rs.simulate(30, True)

    ren = test(rs)

    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)

    # Add the actors to the renderer, set the background and size
    renWin.SetSize(200, 200)

    # This allows the interactor to initalize itself. It has to be
    # called before an event loop.
    iren.Initialize()

    # We'll zoom in a little by accessing the camera and invoking a "Zoom"
    # method on it.
    ren.ResetCamera()
    ren.GetActiveCamera().Zoom(1.5)
    renWin.Render()

    # Start the event loop.
    iren.Start()
