import plantbox as pb
from vtk_tools import *

import numpy as np
import vtk

""" 
VTK Tools, by Daniel Leitner (refurbished 12/2019) 

for vtk to numpy, and numpy to vtk conversions
reading: vtp, writing: msh, dgf, vtp, rsml
"""


def segs_to_polydata(rs, zoom_factor = 10., param_names = ["radius", "type", "creationTime"]):
    """ Creates vtkPolydata from a RootSystem or Plant using segments 
    @param rs             A RootSystem, Plant, or SegmentAnalyser
    @param zoom_factor    The radial zoom factor, since root are sometimes too thin for vizualisation
    @param param_names    Parameter names of scalar fields, that are copied to the polydata    
    @return A vtkPolydata object of the root system
    """
    if isinstance(rs, pb.Organism):
        ana = pb.SegmentAnalyser(rs)  # for Organism like Plant or RootSystem
    else:
        ana = rs
    nodes = np_convert(ana.nodes)
    segs = np_convert(ana.segments)
    points = vtk_points(nodes)
    cells = vtk_cells(segs)
    pd = vtk.vtkPolyData()
    pd.SetPoints(points)
    pd.SetLines(cells)  # check SetPolys
    for n in param_names:
        param = np.array(ana.getParameter(n))
        if param.shape[0] == segs.shape[0]:
            if n == "radius":
                param *= zoom_factor
            data = vtk_data(param)
            data.SetName(n)
            pd.GetCellData().AddArray(data)
        else:
            print("segs_to_polydata: Warning parameter " + n + " is sikpped because of wrong size", param.shape[0], "instead of", segs.shape[0])

    c2p = vtk.vtkCellDataToPointData()
    c2p.SetPassCellData(True)
    c2p.SetInputData(pd)
    c2p.Update()
    return c2p.GetPolyDataOutput()


def uniform_grid(min_, max_, res):
    """ Creates an uniform grid
    
    @return The vtkUniformGrid
    """
    grid = vtk.vtkUniformGrid()
    grid.SetDimensions(int(res[0]) + 1, int(res[1]) + 1, int(res[2]) + 1)  # points
    grid.SetOrigin(min_[0], min_[1], min_[2])
    s = (max_ - min_) / res
    grid.SetSpacing(s[0], s[1], s[2])
    return grid


def render_window(actor, windowName = "", scalarBar = None):
    """ puts a vtk actor on the stage (an interactive window) @name is the window titel 
    @param actor         the (single) actor
    @param windowName    optional
    @param scalarBar     an optional vtkScalarBarActor
    @return The vtkRenderWindow (not sure if this is ever needed)
    """
    colors = vtk.vtkNamedColors()  # Set the background color
    bkg = map(lambda x: x / 255.0, [26, 51, 102, 255])
    colors.SetColor("BkgColor", *bkg)
    ren = vtk.vtkRenderer()  # Set up window with interaction
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    iren = vtk.vtkRenderWindowInteractor()
    # iren.SetInteractorStyle(vtk.vtkInteractorStyleUnicam())  # <- better than default, but maybe we find a better one
    iren.SetRenderWindow(renWin)
    if isinstance(actor, list):
        actors = actor  # plural
    else:
        actors = [actor]  # army of one
    for a in actors:
        ren.AddActor(a)  # Add the actors to the renderer, set the background and size
    if scalarBar is not None:
        ren.AddActor2D(scalarBar)
    ren.SetBackground(colors.GetColor3d("BkgColor"))
    renWin.SetSize(1000, 1000)
    renWin.SetWindowName(windowName)
    iren.Initialize()  # This allows the interactor to initalize itself. It has to be called before an event loop.
    ren.ResetCamera()
    ren.GetActiveCamera().Zoom(1.5)
    renWin.Render()
    iren.Start()  # Start the event loop.
    return renWin


def write_png(renWin, fileName):
    """" Save the current render window in a png
    @param renWin        the vtkRenderWindow 
    @parma fileName      file name without extension
    """
    windowToImageFilter = vtk.vtkWindowToImageFilter();
    windowToImageFilter.SetInput(renWin)
    windowToImageFilter.SetInputBufferTypeToRGBA()  # also record the alpha (transparency) channel
    windowToImageFilter.ReadFrontBufferOff()  # read from the back buffer
    windowToImageFilter.Update()
    writer = vtk.vtkPNGWriter()
    writer.SetFileName(fileName + ".png")
    writer.SetInputConnection(windowToImageFilter.GetOutputPort())
    writer.Write()


def get_lookup_table():
    """ creates a color lookup table 
    @return A vtkLookupTable
    """
#     # Make the lookup table.
#     lut.SetTableRange(scalarRange) = vtk.vtkColorSeries()
#     # Select a color scheme.
#     # colorSeriesEnum = colorSeries.BREWER_DIVERGING_BROWN_BLUE_GREEN_9
#     # colorSeriesEnum = colorSeries.BREWER_DIVERGING_SPECTRAL_10
#     # colorSeriesEnum = colorSeries.BREWER_DIVERGING_SPECTRAL_3
#     # colorSeriesEnum = colorSeries.BREWER_DIVERGING_PURPLE_ORANGE_9
#     # colorSeriesEnum = colorSeries.BREWER_SEQUENTIAL_BLUE_PURPLE_9
#     # colorSeriesEnum = colorSeries.BREWER_SEQUENTIAL_BLUE_GREEN_9
#     colorSeriesEnum = colorSeries.BREWER_QUALITATIVE_SET3
#     # colorSeriesEnum = colorSeries.CITRUS
#     colorSeries.SetColorScheme(colorSeriesEnum)
#     lut = vtk.vtkLookupTable()
#     lut.SetNumberOfTableValues(16)"radius"
#     colorSeries.BuildLookupTable(lut)
#     # lut.SetNanColor(1, 0, 0, 1)
#     lut.SetTableRange([0, 1])
    lut = vtk.vtkLookupTable()
    lut.SetNumberOfTableValues(16)
    lut.SetHueRange(0.0, 1.0)
    lut.Build()
    return lut


def plot_roots(pd, pname, render = True):
    """ renders the root system in an interactive window 
    @param pd         the polydata representing the root system
    @param pname      parameter name of the data to be visualized
    @param render     render in a new interactive window (default = True)
    @return The vtkActor object
    """
    pd.GetPointData().SetActiveScalars("radius")  # for the the filter
    tubeFilter = vtk.vtkTubeFilter()
    tubeFilter.SetInputData(pd)
    tubeFilter.SetNumberOfSides(9)
    tubeFilter.SetVaryRadiusToVaryRadiusByAbsoluteScalar()
    tubeFilter.Update()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(tubeFilter.GetOutputPort())
    mapper.Update()
    mapper.ScalarVisibilityOn();
    mapper.SetScalarModeToUseCellFieldData()  # Cell is not working
    mapper.SetArrayName(pname)
    mapper.SelectColorArray(pname)
    mapper.UseLookupTableScalarRangeOn()

    plantActor = vtk.vtkActor()
    plantActor.SetMapper(mapper)

    lut = get_lookup_table()
    lut.SetTableRange(pd.GetPointData().GetScalars(pname).GetRange())
    mapper.SetLookupTable(lut)

    scalarBar = vtk.vtkScalarBarActor()
    scalarBar.SetLookupTable(lut)
    scalarBar.SetTitle(pname)
#    textProperty = vtk.vtkTextProperty()
#    scalarBar.SetTitleTextProperty(textProperty)
#    scalarBar.SetLabelTextProperty(textProperty)
#    scalarBar.SetAnnotationTextProperty(textProperty)
#    scalarBar = None

    if render:
        render_window(plantActor, pname, scalarBar)
    return plantActor, scalarBar


def plot_mesh_wireframe(grid, p_name, render = True):
    """ Plots the grid as wireframe 
    @param grid         some vtk grid (structured or unstructured)
    @param p_name       parameter to visualize
        
    """
#     bounds = grid.GetBounds()
#     print("mesh bounds", bounds, "[m]")

    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputData(grid)
    mapper.Update()
    mapper.SetArrayName(p_name)
    mapper.SelectColorArray(p_name)
    mapper.UseLookupTableScalarRangeOn()

    meshActor = vtk.vtkActor()
    meshActor.SetMapper(mapper)
    meshActor.GetProperty().SetRepresentationToWireframe();

    lut = get_lookup_table()
    lut.SetTableRange(grid.GetPointData().GetScalars(p_name).GetRange())
    mapper.SetLookupTable(lut)

    scalarBar = vtk.vtkScalarBarActor()
    scalarBar.SetLookupTable(lut)
    scalarBar.SetTitle(p_name)

    if render:
        render_window(meshActor, p_name, scalarBar)  # todo scalarBar (and plot a thing or two)
    return meshActor, scalarBar


def plot_mesh_cuts(ug, p_name, nz = 7):
    """ """
    bounds = ug.GetBounds()
    print(bounds)
    print("z-axis", bounds[4], bounds[5])
    # z-slices (implicit fucntions)
    planes = []
    for i in range(0, nz):
        p = vtk.vtkPlane()
        z = ((bounds[5] - bounds[4]) / (nz + 1)) * (i + 1)
        print(bounds[4] + z)
        p.SetOrigin(bounds[4] + z, 0, 0)
        p.SetNormal(0, 0, 1)
        planes.append(p)

    # create cutter, mappers, and actors
    actors = []
    for i in range(0, nz):
        cutter = vtk.vtkCutter()
        cutter.SetCutFunction(planes[i])
        cutter.SetInputData(ug)
        # cutter.SetInputConnection(cube.GetOutputPort())
        cutter.Update()
        m = vtk.vtkPolyDataMapper()
        m.SetInputConnection(cutter.GetOutputPort())
        m.Update()
        # create plane actor
        a = vtk.vtkActor()
#         planeActor.GetProperty().SetColor(1.0, 1, 0)
#         planeActor.GetProperty().SetLineWidth(2)
        a.SetMapper(m)
        actors.append(a)
    if render:
        render_window_(actors, "Cuts")
    return actors

