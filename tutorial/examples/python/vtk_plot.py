import plantbox as pb
from vtk_tools import *

import time
import numpy as np
import vtk

""" 
VTK Plot, by Daniel Leitner (refurbished 06/2020) 

to make interactive vtk plot of root systems and soil grids
"""


def segs_to_polydata(rs, zoom_factor = 1., param_names = ["age", "radius", "type", "creationTime"]):
    """ Creates vtkPolydata from a RootSystem or Plant using vtkLines to represent the root segments 
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


def render_window(actor, title, scalarBar, bounds):
    """ puts a vtk actor on the stage (renders an interactive window)
    
    @param actor                    a (single) actor, or a list of actors (ensemble)
    @param title                    window title (optional)
    @param scalarBar                one or a list of vtkScalarBarActor (optional)
    @param bounds                   spatial bounds (to set axes actor, and camera position and focal point)
    @return a vtkRenderWindowInteractor     use render_window(...).Start() to start interaction loop, or render_window(...).GetRenderWindow(), to write png
    
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
    Keypress x,y,z,v: various views    
    """
    colors = vtk.vtkNamedColors()  # Set the background color
    ren = vtk.vtkRenderer()  # Set up window with interaction
    ren.SetBackground(colors.GetColor3d("Silver"))

    # Actors
    if isinstance(actor, list):
        actors = actor  # plural
    else:
        actors = [actor]  # army of one
    for a in actors:
        a.GetProperty().BackfaceCullingOff()
        # a.RotateX(-90)  # x y z -> x z y
        ren.AddActor(a)  # Add the actors to the renderer, set the background and size

    if scalarBar:
#         if isinstance(scalarBar, list):
#             c = 0.
#             for sb in scalarBar:  # TODO looks awful
#                 x = sb.GetPosition()
#                 y = (x[0] + c, x[1])
#                 sb.SetPosition(y)
#                 ren.AddActor2D(sb)
#                 c -= 0.2
#         else:
        ren.AddActor2D(scalarBar)

    axes = vtk.vtkAxesActor()
    axes.AxisLabelsOff()  # because i am too lazy to change font size
    translate = vtk.vtkTransform()
    translate.Translate(bounds[0], bounds[2], bounds[4])  # minx, miny, minz
    axes.SetUserTransform(translate)
    ren.AddActor(axes)

    # Camera
    ren.ResetCamera()
    camera = ren.GetActiveCamera()
    camera.ParallelProjectionOn()
    camera.SetFocalPoint([0, 0, 0.5 * (bounds[4] + bounds[5])])
    camera.SetPosition([200, 0, 0.5 * (bounds[4] + bounds[5])])
    camera.SetViewUp(0, 0, 1)
    camera.Azimuth(30)
    camera.Elevation(30)
    camera.OrthogonalizeViewUp()
    camera.SetClippingRange(1, 1000)

    # Render Window
    renWin = vtk.vtkRenderWindow()  # boss
    renWin.SetSize(1000, 1000)
    renWin.SetWindowName(title)
    renWin.AddRenderer(ren)

    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)
    renWin.Render()
    iren.CreateRepeatingTimer(50)  # [ms] 0.5 s in case a timer event is interested
    iren.AddObserver('KeyPressEvent', lambda obj, ev :keypress_callback_(obj, ev, bounds), 1.0)
    iren.Initialize()  # This allows the interactor to initalize itself. It has to be called before an event loop.
    for a in ren.GetActors():
        a.Modified()  #
    renWin.Render()
    return iren


def keypress_callback_(obj, ev, bounds):
    """ adds the functionality to make a screenshot by pressing 'g', 
    and to change view to axis aligned plots (by 'x', 'y', 'z', 'v') """
    key = obj.GetKeySym()
    if key == 'g':
        renWin = obj.GetRenderWindow()
        file_name = renWin.GetWindowName()
        write_png(renWin, file_name)
        print("saved", file_name + ".png")
    if key == 'x' or key == 'y' or key == 'z' or key == 'v':
        renWin = obj.GetRenderWindow()
        ren = renWin.GetRenderers().GetItemAsObject(0)
        camera = ren.GetActiveCamera()
        if key == 'x':
            camera.SetPosition([100, 0, 0.5 * (bounds[4] + bounds[5])])
            camera.SetViewUp(0, 0, 1)
            print("y-z plot")
        if key == 'y':
            camera.SetPosition([0, 100, 0.5 * (bounds[4] + bounds[5])])
            camera.SetViewUp(0, 0, 1)
            print("x-z plot")
        if key == 'z':
            camera.SetPosition([0, 0, 100])  #
            camera.SetViewUp(0, 1, 0)
            print("x-y plot")
        if key == 'v':
            camera.SetPosition([100, 0, 0.5 * (bounds[4] + bounds[5])])
            camera.SetViewUp(0, 0, 1)
            camera.Azimuth(30)
            camera.Elevation(30)
            print("oblique plot")
        camera.OrthogonalizeViewUp()
        renWin.Render()


def write_png(renWin, fileName):
    """" Save the current render window in a png (e.g. from vtkRenderWindowInteractor.GetRenderWindow())
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


def create_lookup_table(tableIdx = 15, numberOfColors = 256):
    """ creates a color lookup table 
    @param tableIdx          index of the predefined color table, see VTKColorSeriesPatches.html
    @param numberOfColors    number of colors interpolated from the predefined table
    @return A vtkLookupTable
    """
    colorSeries = vtk.vtkColorSeries()
    colorSeries.SetColorScheme(tableIdx)
    lut_ = colorSeries.CreateLookupTable(vtk.vtkColorSeries.ORDINAL)
    n = lut_.GetNumberOfTableValues ()

    lut = vtk.vtkLookupTable()
    lut.SetNumberOfTableValues(numberOfColors)
    for i in range(0, numberOfColors):
        psi = (n - 1) * (float(i) / numberOfColors)
        i0 = np.floor(psi)
        theta = psi - i0
        col = (1 - theta) * np.array(lut_.GetTableValue(int(i0))) + theta * np.array(lut_.GetTableValue(int(i0) + 1))
        lut.SetTableValue(i, col)

    return lut


def create_scalar_bar(lut, grid = None, p_name = ""):
    """ creates a vtkScalarBarActor, for a vtkLookupTable, sets hte active scalar to p_name
    @param lut         vtkLookupTable
    @param grid        the grid the scalar bar will be used on (to automatically determine the scalar range)
    @param p_name      name of the cell data or point data, from which the range is determined
    @return a vtkScalarBarActor
    """
    if grid != None and p_name != "":
        range = [0, 1]
        a = grid.GetCellData().GetAbstractArray(p_name)
        if a:
            range = a.GetRange()
            grid.GetCellData().SetActiveScalars(p_name)
        else:
            a = grid.GetPointData().GetAbstractArray(p_name)
            grid.GetPointData().SetActiveScalars(p_name)
            if a:
                range = a.GetRange()
        lut.SetTableRange(range)

    scalarBar = vtk.vtkScalarBarActor()
    scalarBar.SetLookupTable(lut)
    scalarBar.SetTitle(p_name)
    scalarBar.SetDrawAnnotations(False)
    textProperty = vtk.vtkTextProperty()
    textProperty.SetFontSize(30)
    scalarBar.SetAnnotationTextProperty(textProperty)
    scalarBar.SetTitleTextProperty(textProperty)
    scalarBar.SetLabelTextProperty(textProperty)
    scalarBar.AnnotationTextScalingOff()
    scalarBar.SetUnconstrainedFontSize(True)

    return scalarBar


def plot_roots(pd, p_name, win_title = "", render = True):
    """ plots the root system 
    @param pd         the polydata representing the root system (lines, or polylines)
    @param p_name     parameter name of the data to be visualized
    @param win_title  the windows titles (optionally, defaults to p_name)
    @param render     render in a new interactive window (default = True)
    @return a tuple of a vtkActor and the corresponding color bar vtkScalarBarActor
    """
    if isinstance(pd, pb.RootSystem):
        pd = segs_to_polydata(pd, 1., [p_name, "radius"])

    if isinstance(pd, pb.SegmentAnalyser):
        pd = segs_to_polydata(pd, 1., [p_name, "radius"])

    if win_title == "":
        win_title = p_name

    pd.GetPointData().SetActiveScalars("radius")  # for the the filter
    tubeFilter = vtk.vtkTubeFilter()
    tubeFilter.SetInputData(pd)
    tubeFilter.SetNumberOfSides(9)  #
    tubeFilter.SetVaryRadiusToVaryRadiusByAbsoluteScalar()
    tubeFilter.Update()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(tubeFilter.GetOutputPort())
    mapper.Update()
    mapper.ScalarVisibilityOn();
    mapper.SetScalarModeToUseCellFieldData()  # maybe because radius is active scalar in point data?
    mapper.SetArrayName(p_name)
    mapper.SelectColorArray(p_name)
    mapper.UseLookupTableScalarRangeOn()
    plantActor = vtk.vtkActor()
    plantActor.SetMapper(mapper)

    lut = create_lookup_table()  # 24
    scalar_bar = create_scalar_bar(lut, pd, p_name)  # vtkScalarBarActor
    mapper.SetLookupTable(lut)

    if render:
        render_window(plantActor, win_title, scalar_bar, pd.GetBounds()).Start()
    return plantActor, scalar_bar


def plot_mesh(grid, p_name, win_title = "", render = True):
    """ Plots the grid as wireframe
    @param grid         some vtk grid (structured or unstructured)
    @param p_name       parameter to visualize
    @param win_title  the windows titles (optionally, defaults to p_name)
    @param render     render in a new interactive window (default = True)    
    @return a tuple of a vtkActor and the corresponding color bar vtkScalarBarActor
    """
    if win_title == "":
        win_title = p_name

    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputData(grid)
    mapper.Update()
    mapper.SetArrayName(p_name)
    mapper.SelectColorArray(p_name)  # ? to choosvtkScalarBarActore cell data or point data
    mapper.UseLookupTableScalarRangeOn()
    meshActor = vtk.vtkActor()
    meshActor.SetMapper(mapper)
    meshActor.GetProperty().SetRepresentationToWireframe();

    lut = create_lookup_table()
    scalar_bar = create_scalar_bar(lut, grid, p_name)  # vtkScalarBarActor
    mapper.SetLookupTable(lut)

    if render:
        render_window(meshActor, win_title, scalar_bar, grid.GetBounds()).Start()
    return meshActor, scalar_bar


def plot_mesh_cuts(grid, p_name, nz = 3, win_title = "", render = True):
    """ plots orthogonal nz vertical cuts z[:-1] (xy-planes), with z = linspace(min_z, max_z, nz+1), 
    and two additonal sclices at x=0 (yz-plane), y=0 (xz-plane)          
    @param grid         some vtk grid (structured or unstructured)
    @param p_name       parameter to visualize
    @param nz           number of vertical slices
    @param win_title    the windows titles (optionally, defaults to p_name)
    @param render       render in a new interactive window (default = True)
    @return a tuple of a list of vtkActors and a single corresponding color bar vtkScalarBarActor    
    """
    if win_title == "":
        win_title = p_name

    eps = 1.e-2
    planes = []  # create the cut planes
    bounds = grid.GetBounds()
    z = np.linspace(bounds[4] + eps, bounds[5], nz + 1)
    for i in range(0, nz):  # z-slices (implicit functions)
        p = vtk.vtkPlane()
        p.SetOrigin(0, 0, z[i])
        p.SetNormal(0, 0, 1)
        planes.append(p)
    for n in [(1, 0, 0), (0, 1, 0)]:
        p = vtk.vtkPlane()
        p.SetOrigin(bounds[0] + eps, bounds[2] + eps, bounds[4])
        p.SetNormal(n[0], n[1], n[2])
        planes.append(p)

    lut = create_lookup_table()
    scalar_bar = create_scalar_bar(lut, grid, p_name)

    actors = []  # create cutter, mappers, and actors
    for p in planes:
        cutter = vtk.vtkCutter()
        cutter.SetInputData(grid)
        # cutter.SetInputConnection(tubeFilter.GetOutputPort()) # for root system (tube plot)
        cutter.SetCutFunction(p)
        cutter.Update()
        m = vtk.vtkDataSetMapper()
        m.SetInputConnection(cutter.GetOutputPort())
        m.Update()
        m.SetArrayName(p_name)
        m.SelectColorArray(p_name)
        m.UseLookupTableScalarRangeOn()
        m.SetLookupTable(lut)
        m.SetColorModeToMapScalars();
        a = vtk.vtkActor()  # create plane actor
#         a.GetProperty().SetColor(1.0, 1, 0)
#         a.GetProperty().SetLineWidth(2)
        a.SetMapper(m)
        actors.append(a)

    if render:
        render_window(actors, win_title, scalar_bar, grid.GetBounds()).Start()

    return actors, scalar_bar


def plot_roots_and_soil(rs, pname :str, rp, s, periodic :bool, min_b, max_b, cell_number, filename :str):
    """ Plots soil slices and roots, additionally saves both grids as files
    @param rs            some Organism (e.g. RootSystem, MappedRootSystem, ...) or MappedSegments
    @param pname         root and soil parameter that will be visualized ("pressure head", or "water content")
    @param rp            root parameter segment data (will be added)
    @param periodic      if yes the root system will be mapped into the domain 
    @param min_b         minimum of domain boundaries
    @param max_b         maximum of domain boundaries    
    @param cell_number   domain resolution
    @param filename      file name (without extension)
    """
    ana = pb.SegmentAnalyser(rs)
    ana.addData(pname, rp)
    if periodic:
        w = max_b - min_b
        ana.mapPeriodic(w[0], w[1])
    pd = segs_to_polydata(ana, 1., ["radius", "subType", "creationTime", pname])

    soil_grid = uniform_grid(np.array(min_b), np.array(max_b), np.array(cell_number))
    soil_water_content = vtk_data(np.array(s.getWaterContent()))
    soil_water_content.SetName("water content")
    soil_grid.GetCellData().AddArray(soil_water_content)
    soil_pressure = vtk_data(np.array(s.getSolutionHead()))
    soil_pressure.SetName("pressure head")  # in macroscopic soil
    soil_grid.GetCellData().AddArray(soil_pressure)

    rootActor, rootCBar = plot_roots(pd, pname, "", False)
    meshActors, meshCBar = plot_mesh_cuts(soil_grid, pname, 4, "", False)
    lut = meshActors[-1].GetMapper().GetLookupTable()  # same same
    rootActor.GetMapper().SetLookupTable(lut)
    meshActors.extend([rootActor])
    render_window(meshActors, filename, meshCBar, pd.GetBounds()).Start()

    if filename:
        path = "results/"
        write_vtp(path + filename + ".vtp", pd)
        write_vtu(path + filename + ".vtu", soil_grid)


def plot_roots_and_soil_files(filename : str, pname :str):
    """ Plots soil slices and roots from two files (one vtp and one vtu), created by plot_roots_and_soil() 
    @param filename      file name (without extension)
    @param pname         root and soil parameter that will be visualized ("pressure head", or "water content")    
    """
    path = "results/"
    pd = read_vtp(path + filename + ".vtp")
    soil_grid = read_rect_vtu(path + filename + ".vtp")
    rootActor, rootCBar = plot_roots(pd, pname, "", False)
    meshActors, meshCBar = plot_mesh_cuts(soil_grid, pname, 4, "", False)
    lut = meshActors[-1].GetMapper().GetLookupTable()  # same same
    rootActor.GetMapper().SetLookupTable(lut)
    meshActors.extend([rootActor])
    render_window(meshActors, filename, meshCBar).Start()


class AnimateRoots:
    """ class to make an interactive animation """

    def __init__(self, rootsystem = None):
        self.rootsystem = rootsystem
        self.root_name = "subType"
        # self.soil_name = "subType"
        #
        self.min = None
        self.max = None
        self.res = None
        self.soil_data = True  # soil data
        self.cuts = False  # Wireframe, or cuts
        #
        self.actors = []
        self.iren = None
        self.color_bar = None
        self.bounds = None
        self.file = None

    def __del__(self):
        if self.file:
            self.writer.End()

    def start(self, axis = 'x', avi_file = None):
        """ creates plot and adjusts camera """
        self.create_root_actors()
        self.create_soil_actors()
        self.iren = render_window(self.actors, "AnimateRoots", self.color_bar, self.bounds)
        renWin = self.iren.GetRenderWindow()
        ren = renWin.GetRenderers().GetItemAsObject(0)
        camera = ren.GetActiveCamera()
        if axis == 'x':
            camera.SetPosition([100, 0, 0.5 * (self.bounds[4] + self.bounds[5])])
            camera.SetViewUp(0, 0, 1)
        if axis == 'y':
            camera.SetPosition([0, 100, 0.5 * (self.bounds[4] + self.bounds[5])])
            camera.SetViewUp(0, 0, 1)
        if axis == 'z':
            camera.SetPosition([0, 0, 100])
            camera.SetViewUp(0, 1, 0)
        if axis == 'v':
            camera.SetPosition([100, 0, 0.5 * (self.bounds[4] + self.bounds[5])])
            camera.SetViewUp(0, 0, 1)
            camera.Azimuth(30)
            camera.Elevation(30)
        if self.file:
            self.windowToImageFilter = vtk.vtkWindowToImageFilter();
            self.windowToImageFilter.SetInput(renWin)
            self.windowToImageFilter.SetInputBufferTypeToRGB()
            self.windowToImageFilter.ReadFrontBufferOff()  # read from the back buffer
            self.windowToImageFilter.Update()
            w = vtk.vtkOggTheoraWriter()
            w.SetFileName(self.file + ".ogv")
            w.SetInputConnection(self.windowToImageFilter.GetOutputPort())
            # w.SetCompressorFourCC("H264")  # feeling lucky
            w.Start()  #
            self.writer = w

    def update(self):
        """ animation call back function (called every 0.1 second) """
        renWin = self.iren.GetRenderWindow()
        ren = renWin.GetRenderers().GetFirstRenderer()
        for a in self.actors:
            ren.RemoveActor(a)
        for a in ren.GetActors2D():
            ren.RemoveActor2D(a)
        ren.AddActor2D(self.color_bar)
        self.actors = []
        self.create_root_actors()
        self.create_soil_actors()
        for a in self.actors:
            ren.AddActor(a)

        self.iren.Render()
        if self.file:
            self.windowToImageFilter.Modified()
            self.windowToImageFilter.Update()
            self.writer.Write()

    def create_root_actors(self):
        if self.rootsystem:
           pd = segs_to_polydata(self.rootsystem, 1., [self.root_name, "radius"])
           newRootActor, rootCBar = plot_roots(pd, self.root_name, "", False)
           self.actors.append(newRootActor)
           self.color_bar = rootCBar
           self.bounds = pd.GetBounds()

    def create_soil_actors(self):
        if self.soil_data:
            if self.cuts:
                pass
                # meshActor, meshCBar = plot_mesh_cuts(grid, p_name, nz = 3, win_title = "", False):
            else:
                grid = uniform_grid(self.min, self.max, self.res)
                meshActor, meshCBar = plot_mesh(grid, "", "", False)
            self.actors.append(meshActor)
            self.bounds = grid.GetBounds()

