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
    @param min_    minimum of bounding rectangle
    @param max_    maximum of bounding rectangle
    @param res_    cell resolution
    @return The vtkUniformGrid
    """
    grid = vtk.vtkUniformGrid()
    grid.SetDimensions(int(res[0]) + 1, int(res[1]) + 1, int(res[2]) + 1)  # cell to point resolution
    grid.SetOrigin(min_[0], min_[1], min_[2])
    s = (max_ - min_) / res
    grid.SetSpacing(s[0], s[1], s[2])
    return grid


def render_window(actor, title = "", scalarBar = None):
    """ puts a vtk actor on the stage (an interactive window) @name is the window titel 
    
    @param actor         the (single) actor
    @param windowName    optional
    @param scalarBar     an optional vtkScalarBarActor
    @return The vtkRenderWindow (not sure if this is ever needed)
    
    (built in)
    Keypress j / Keypress t: toggle between joystick (position sensitive) and trackball (motion sensitive) styles. In joystick style, motion occurs continuously as long as a mouse button is pressed. In trackball style, motion occurs when the mouse button is pressed and the mouse pointer moves.
    Keypress c / Keypress a: toggle between camera and actor modes. In camera mode, mouse events affect the camera position and focal point. In actor mode, mouse events affect the actor that is under the mouse pointer.
    Button 1: rotate the camera around its focal point (if camera mode) or rotate the actor around its origin (if actor mode). The rotation is in the direction defined from the center of the renderer's viewport towards the mouse position. In joystick mode, the magnitude of the rotation is determined by the distance the mouse is from the center of the render window.
    Button 2: pan the camera (if camera mode) or translate the actor (if actor mode). In joystick mode, the direction of pan or translation is from the center of the viewport towards the mouse position. In trackball mode, the direction of motion is the direction the mouse moves. (Note: with 2-button mice, pan is defined as <Shift>-Button 1.)
    Button 3: zoom the camera (if camera mode) or scale the actor (if actor mode). Zoom in/increase scale if the mouse position is in the top half of the viewport; zoom out/decrease scale if the mouse position is in the bottom half. In joystick mode, the amount of zoom is controlled by the distance of the mouse pointer from the horizontal centerline of the window.
    Keypress 3: toggle the render window into and out of stereo mode. By default, red-blue stereo pairs are created. Some systems support Crystal Eyes LCD stereo glasses; you have to invoke SetStereoTypeToCrystalEyes() on the rendering window.
    Keypress e: exit the application.
    Keypress f: fly to the picked point
    Keypress p: perform a pick operation. The render window interactor has an internal instance of vtkCellPicker that it uses to pick.
    Keypress r: reset the camera view along the current view direction. Centers the actors and moves the camera so that all actors are visible.
    Keypress s: modify the representation of all actors so that they are surfaces.
    Keypress u: invoke the user-defined function. Typically, this keypress will bring up an interactor that you can type commands in. Typing u calls UserCallBack() on the vtkRenderWindowInteractor, which invokes a vtkCommand::UserEvent. In other words, to define a user-defined callback, just add an observer to the vtkCommand::UserEvent on the vtkRenderWindowInteractor object.
    Keypress w: modify the representation of all actors so that they are wireframe.
    
    (additional)
    Keypress g: save as png    
    """
    colors = vtk.vtkNamedColors()  # Set the background color
    bkg = map(lambda x: x / 255.0, [26, 51, 102, 255])
    colors.SetColor("BkgColor", *bkg)
    ren = vtk.vtkRenderer()  # Set up window with interaction
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    iren = vtk.vtkRenderWindowInteractor()
    # iren.SetInteractorStyle(vtk.vtkInteractorStyleUnicam())  # <- better than default, but maybe we find a better one
    if isinstance(actor, list):
        actors = actor  # plural
    else:
        actors = [actor]  # army of one
    for a in actors:
        a.RotateX(-90)
        ren.AddActor(a)  # Add the actors to the renderer, set the background and size
    if scalarBar is not None:
        ren.AddActor2D(scalarBar)
    ren.SetBackground(colors.GetColor3d("Silver"))  #
    renWin.SetSize(1000, 1000)
    renWin.SetWindowName(title)
    ren.ResetCamera()
    ren.GetActiveCamera().Zoom(1.5)
    ren.GetActiveCamera().ParallelProjectionOn()
    iren.SetRenderWindow(renWin)
    renWin.Render()
    iren.Initialize()  # This allows the interactor to initalize itself. It has to be called before an event loop.
    iren.CreateRepeatingTimer(500)  # [ms] 0.5 s in case a timer event is interested
    iren.AddObserver('KeyPressEvent', keypress_callback_, 1.0)
    return iren


def keypress_callback_(obj, ev):
    """ adds the functionality to make a screenshot by pressing 'g' """
    key = obj.GetKeySym()
    if key == 'g':
        renWin = obj.GetRenderWindow()
        file_name = renWin.GetWindowName()
        write_png(renWin, file_name)
        print("saved", file_name + ".png")


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


def get_lookup_table(colorSeriesEnum = 15):
    """ creates a color lookup table 
    @param colorSeriesEnum      the number of the predefined color table, see VTKColorSeriesPatches.html
    @return A vtkLookupTable
    
    @todo I don't know how to modify the number of colors used (or how to interpolate between)
    """
    colorSeries = vtk.vtkColorSeries()
    colorSeries.SetColorScheme(colorSeriesEnum)
    lut = vtk.vtkLookupTable()
    colorSeries.BuildLookupTable(lut, vtk.vtkColorSeries.ORDINAL)
    return lut


def plot_roots(pd, p_name, render = True):
    """ renders the root system in an interactive window 
    @param pd         the polydata representing the root system
    @param p_name      parameter name of the data to be visualized
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
    mapper.SetArrayName(p_name)
    mapper.SelectColorArray(p_name)
    mapper.UseLookupTableScalarRangeOn()

    plantActor = vtk.vtkActor()
    plantActor.SetMapper(mapper)

    lut = get_lookup_table(24)  # 24= Brewer Diverging Brown-Blue-Green (11)
    lut.SetTableRange(pd.GetPointData().GetScalars(p_name).GetRange())
    mapper.SetLookupTable(lut)

    scalarBar = vtk.vtkScalarBarActor()
    scalarBar.SetLookupTable(lut)
    scalarBar.SetTitle(p_name)
    textProperty = vtk.vtkTextProperty()
    textProperty.SetFontSize(1)
    scalarBar.SetTitleTextProperty(textProperty)
    scalarBar.SetLabelTextProperty(textProperty)
#    scalarBar.SetAnnotationTextProperty(textProperty)

    if render:
        render_window(plantActor, pname, scalarBar)
    return plantActor, scalarBar


def plot_mesh(grid, p_name, win_title = "", render = True):
    """ Plots the grid as wireframe
    @param grid         some vtk grid (structured or unstructured)
    @param p_name       parameter to visualize
        
    """
#     bounds = grid.GetBounds()
#     print("mesh bounds", bounds, "[m]")
    if win_title == "":
        win_title = p_name

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
    if p_name != "":
        lut.SetTableRange(grid.GetPointData().GetScalars(p_name).GetRange())
    mapper.SetLookupTable(lut)

    scalarBar = vtk.vtkScalarBarActor()
    scalarBar.SetLookupTable(lut)
    scalarBar.SetTitle(p_name)
#     textProperty = vtk.vtkTextProperty()
#     textProperty.SetFontSize(10)
#     textProperty.UseTightBoundingBoxOff()
#     scalarBar.GetLabelTextProperty().BoldOff()
#     scalarBar.SetTitleTextProperty(textProperty)
#     scalarBar.SetLabelTextProperty(textProperty)
#     scalarBar.GetLabelTextProperty().SetFontSize(1)

    if render:
        render_window(meshActor, win_title, scalarBar)
    return meshActor, scalarBar


def plot_mesh_cuts(pd, p_name, nz = 7):
    """ """

    # Create a cube
    cube = vtk.vtkCubeSource()
    cube.SetXLength(40)
    cube.SetYLength(30)
    cube.SetZLength(20)
#     cubeMapper = vtk.vtkPolyDataMapper()
#     cubeMapper.SetInputConnection(cube.GetOutputPort())

    pd.GetPointData().SetActiveScalars("radius")  # for the the filter
    tubeFilter = vtk.vtkTubeFilter()
    tubeFilter.SetInputData(pd)
    tubeFilter.SetNumberOfSides(9)
    tubeFilter.SetVaryRadiusToVaryRadiusByAbsoluteScalar()
    tubeFilter.Update()

    bounds = pd.GetBounds()
    planes = []
    for i in range(0, nz):  # z-slices (implicit functions)
        p = vtk.vtkPlane()
        z = ((bounds[5] - bounds[4]) / (nz + 1)) * (i + 1)
        p.SetOrigin(0., 0., bounds[4] + z)
        p.SetNormal(0., 0., 1)
        planes.append(p)

    lut = get_lookup_table(24)  # 24= Brewer Diverging Brown-Blue-Green (11)
    # lut.SetTableRange(pd.GetPointData().GetScalars(p_name).GetRange())

    # create cutter, mappers, and actors
    actors = []
    for i in range(0, nz):
        cutter = vtk.vtkCutter()
        # cutter.SetInputData(pd)
        cutter.SetCutFunction(planes[i])
        cutter.SetInputConnection(tubeFilter.GetOutputPort())
        cutter.Update()
        m = vtk.vtkPolyDataMapper()
        m.SetInputConnection(cutter.GetOutputPort())
        m.Update()
        m.SetLookupTable(lut)
        a = vtk.vtkActor()  # create plane actor
        a.GetProperty().SetColor(1.0, 1, 0)
        a.GetProperty().SetLineWidth(2)
        a.SetMapper(m)
        actors.append(a)
    return actors

