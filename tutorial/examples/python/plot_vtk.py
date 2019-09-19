import vtk
from vtk import *
import py_rootbox as rb
from rb_tools import *


def vtkPoints(p):
    """ Creates vtkPoints from an numpy array
    """
    da = vtk.vtkDataArray.CreateDataArray(vtk.VTK_DOUBLE)
    da.SetNumberOfComponents(3)  # vtk point dimension is always 3
    da.SetNumberOfTuples(p.shape[0])
    for i in range(0, p.shape[0]):
        if p.shape[1] == 2:
            da.InsertTuple3(i, p[i, 0], p[i, 1], 0.)
        elif p.shape[1] == 3:
            da.InsertTuple3(i, p[i, 0], p[i, 1], p[i, 2])
    points = vtk.vtkPoints()
    points.SetData(da)
    return points


def vtkCells(t):
    """ Creates vtkCells from an numpy array
    """
    cellArray = vtk.vtkCellArray()
    for vert in t:
        if t.shape[1] == 2:
            tetra = vtk.vtkLine()
        if t.shape[1] == 3:
            tetra = vtk.vtkTriangle()
        elif t.shape[1] == 4:
            tetra = vtk.vtkTetra()
        for i, v in enumerate(vert):
            tetra.GetPointIds().SetId(i, int(v))
        cellArray.InsertNextCell(tetra)
    return cellArray


def vtkData(d):
    """ Creates a vtkDataArray from an numpy array, e.g. grid.GetCellData().SetScalars(vtkData(celldata))
    """
    da = vtk.vtkDataArray.CreateDataArray(vtk.VTK_DOUBLE)
    noc = d.shape[1]
    da.SetNumberOfComponents(noc)
    da.SetNumberOfTuples(d.shape[0])
    for i in range(0, d.shape[0]):
        if  noc == 1:
            da.InsertTuple1(i, d[i, 0])
        elif noc == 2:
            da.InsertTuple2(i, d[i, 0], d[i, 1])
        elif noc == 3:
            da.InsertTuple3(i, d[i, 0], d[i, 1], d[i, 2])
        elif noc == 4:
            da.InsertTuple4(i, d[i, 0], d[i, 1], d[i, 2], d[i, 3])
    return da


# Simulate something
rs = rb.RootSystem()
name = "Zea_mays_1_Leitner_2010"  # "Anagallis_femina_Leitner_2010"
rs.openFile(name)
rs.initialize()
rs.simulate(30, True)

# Export to VTK
nodes = vv2a(rs.getNodes())
segs = seg2a(rs.getSegments())
types = v2ai(rs.getParameter("type"))
print(len(types), types, types.shape)
points = vtkPoints(nodes)
cells = vtkCells(segs)
data = vtkData(types)
pd = vtk.vtkPolyData()
pd.SetPoints(points)
pd.SetLines(cells)
pd.GetCellData().SetScalars(types)

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
iren.SetInteractorStyle(vtk.vtkInteractorStyleUnicam())  # <- better than default, but maybe we find a better one
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
