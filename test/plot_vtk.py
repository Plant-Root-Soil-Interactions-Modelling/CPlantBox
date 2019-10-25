import vtk
from vtk import *
import sys
sys.path.append("..")
import plantbox as pb


def vtkPoints(p_):
    """ Creates vtkPoints from a list of Vector3d """
    da = vtk.vtkDataArray.CreateDataArray(vtk.VTK_DOUBLE)
    da.SetNumberOfComponents(3)  # vtk point dimension is always 3
    da.SetNumberOfTuples(len(p_))
    for i, p in enumerate(p_):
        da.InsertTuple3(i, p.x, p.y, p.z)
    points = vtk.vtkPoints()
    points.SetData(da)
    return points

def vtkCells(t):
    """ Creates vtkCells from an numpy array """
    cellArray = vtk.vtkCellArray()
    for i, l in enumerate(t):
        tetra = vtk.vtkLine()
        tetra.GetPointIds().SetId(0, l.x)
        tetra.GetPointIds().SetId(1, l.y)        
        cellArray.InsertNextCell(tetra)
    return cellArray


def vtkData(d):
    """ Creates a vtkDataArray from an numpy array, e.g. grid.GetCellData().SetScalars(vtkData(celldata)) """
    da = vtk.vtkDataArray.CreateDataArray(vtk.VTK_DOUBLE)
    da.SetNumberOfComponents(1)
    da.SetNumberOfTuples(len(d))
    for i in range(0, len(d)):
        da.InsertTuple1(i, d[i])
    return da


# Simulate something
rs = pb.RootSystem()
path = "../modelparameter/rootsystem/"
name = "Zea_mays_1_Leitner_2010.xml"  # "Anagallis_femina_Leitner_2010"
rs.readParameters(path+name)
rs.initialize()
rs.simulate(30, True)

# Export to VTK
nodes = rs.getNodes()
segs = rs.getSegments()
types = rs.getParameter("type")
points = vtkPoints(nodes)
cells = vtkCells(segs)
data = vtkData(types)
pd = vtk.vtkPolyData()
pd.SetPoints(points)
pd.SetLines(cells)
pd.GetCellData().SetScalars(data)

# Set the background color
colors = vtk.vtkNamedColors()
bkg = map(lambda x: x / 255.0, [26, 51, 102, 255])
colors.SetColor("BkgColor", *bkg)

# Set up VTK
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputData(pd)
plantActor = vtk.vtkActor()
plantActor.SetMapper(mapper)
plantActor.RotateX(-90.0)

# Set up window with interaction
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
iren = vtk.vtkRenderWindowInteractor()
# iren.SetInteractorStyle(vtk.vtkInteractorStyleUnicam())  # <- better than default, but maybe we find a better one
iren.SetRenderWindow(renWin)

# Add the actors to the renderer, set the background and size
ren.AddActor(plantActor)
ren.SetBackground(colors.GetColor3d("BkgColor"))
renWin.SetSize(800, 800)
renWin.SetWindowName(name)

# This allows the interactor to initalize itself. It has to be called before an event loop.
iren.Initialize()
ren.ResetCamera()
ren.GetActiveCamera().Zoom(1.5)
renWin.Render()
iren.Start()  # Start the event loop.
