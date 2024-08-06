import plantbox as pb
from visualisation.vtk_tools import *

from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
import time
import numpy as np
import vtk
from IPython.display import Image, display

"""
VTK Plot, by Daniel Leitner (refurbished 06/2020)

to make interactive vtk plot of root systems and soil grids
"""


def plot_leaf(leaf):
    """
        plots a single leaf in an interactive window (for debugging leafs)
    """
    leaf_points = vtk.vtkPoints()
    leaf_polys = vtk.vtkCellArray()  # describing the leaf surface area
    create_leaf(leaf, leaf_points, leaf_polys)
    polyData = vtk.vtkPolyData()
    polyData.SetPoints(leaf_points)
    polyData.SetPolys(leaf_polys)
    colors = vtk.vtkNamedColors()
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polyData)
    actor = vtk.vtkActor()
    actor.SetMapper(mapper);
    actor.GetProperty().SetColor(colors.GetColor3d("Green"))
    render_window([actor], "plot_plant", [], [-10, 10, -10, 10, -10, 10]).Start()


def plot_plant(plant, p_name, vals=[], render = True, interactiveImage = True, filename = "plot_plant", range_ = None):
    """
        @param interactiveImage         make image interactive or static (should be static for google Colab)
        plots a whole plant as a tube plot, and additionally plot leaf surface areas as polygons
    """
    # plant as tube plot
    if not isinstance(p_name, list):
        p_name = [p_name]
    p_name_ = ["radius", "organType", "creationTime"]+ p_name
    pd = segs_to_polydata(plant, 1.,p_name_,  vals)  # poly data
    tube_plot_actor, color_bar, lut = plot_roots(pd, p_name_[-1], "", render = False)

    # leafes as polygons
    leaf_points = vtk.vtkPoints()
    leaf_polys = vtk.vtkCellArray()  # describing the leaf surface area

    leafes = plant.getOrgans(pb.leaf)
    globalIdx_y = []
    for l in leafes:
        globalIdx_y =globalIdx_y + create_leaf_(l, leaf_points, leaf_polys)
    globalIdx_y = np.array(globalIdx_y)
    polyData = vtk.vtkPolyData()
    polyData.SetPoints(leaf_points)
    polyData.SetPolys(leaf_polys)
    numnodes = plant.getNumberOfNodes()
    if (len(vals)>0) and (len(leafes)>0):    
        if (not isinstance(vals, list)) and (not isinstance(vals[0], type(np.array([])))):
            vals = [vals]
        for i in range(len(vals)):
            vals_ = vals[-1-i]
            p_name_leaf = p_name[-1-i]
            #print(p_name_leaf, vals_)
            if len(vals_) == numnodes:
                param = vals_[globalIdx_y] #select data for leaf
            else :
                if len(vals_) == (numnodes-1):
                    #print(vals_,globalIdx_y, type(globalIdx_y))
                    param = vals_[globalIdx_y -1]
            data = vtk_data(param)
            data.SetName(p_name_leaf)
            polyData.GetCellData().AddArray(data)
    else:
        colors = vtk.vtkNamedColors()


    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polyData)
    
    mapper.ScalarVisibilityOn();
    mapper.SetScalarModeToUseCellFieldData()  # maybe because radius is active scalar in point data?
    #mapper.SetArrayName(p_name_)
    mapper.SelectColorArray(p_name_[-1])
    mapper.UseLookupTableScalarRangeOn()
    mapper.SetLookupTable(lut)
    actor = vtk.vtkActor()
    actor.SetMapper(mapper);
    
    if not ((len(vals)>0) and (len(leafes)>0)):
        actor.GetProperty().SetColor(colors.GetColor3d("Green"))

    if render:
        write_vtp(filename+"_leaf.vtp", polyData)
        write_vtp(filename+".vtp", pd)
        #ren = render_window([tube_plot_actor, actor], "plot_plant", color_bar, tube_plot_actor.GetBounds(), interactiveImage)
        #if interactiveImage:
        #    ren.Start()

    return [tube_plot_actor, actor], color_bar


def create_leaf_(leaf, leaf_points, leaf_polys):
    """ used by plot plant """
    offs = leaf_points.GetNumberOfPoints()
    globalIdx_y = [] #index of y node
    for i in range(0, leaf.getNumberOfNodes() - 1):  #

        ln1 = leaf.getLeafVis(i)
        ln2 = leaf.getLeafVis(i + 1)

        if len(ln1) > 0 or len(ln2) > 0:
            n1 = leaf.getNode(i)
            n2 = leaf.getNode(i + 1)
            glIdx = leaf.getNodeId(i +1)
            if len(ln1) == 2 and len(ln2) == 2:  # normal case
                offs = add_quad_(n1, ln1[0], ln2[0], n2, leaf_points, leaf_polys, offs)
                offs = add_quad_(n1, ln1[1], ln2[1], n2, leaf_points, leaf_polys, offs)
                globalIdx_y = globalIdx_y + [glIdx,glIdx]
            elif len(ln1) == 6 and len(ln2) == 6:  # convex case
                offs = add_quad_(n1, ln1[0], ln2[0], n2, leaf_points, leaf_polys, offs)
                offs = add_quad_(ln1[1], ln1[2], ln2[2], ln2[1], leaf_points, leaf_polys, offs)
                offs = add_quad_(n1, ln1[3], ln2[3], n2, leaf_points, leaf_polys, offs)
                offs = add_quad_(ln1[4], ln1[5], ln2[5], ln2[4], leaf_points, leaf_polys, offs)
                globalIdx_y = globalIdx_y + [glIdx,glIdx,glIdx,glIdx]
            elif len(ln1) == 2 and len(ln2) == 6:  # normal to convex case
                x1 = leaf.getLeafVisX(i)
                x2 = leaf.getLeafVisX(i + 1)
                offs = add_quad_(n1, ln1[0], ln2[0], n2, leaf_points, leaf_polys, offs)
                offs = add_quad_(n1, ln1[1], ln2[3], n2, leaf_points, leaf_polys, offs)
                globalIdx_y = globalIdx_y + [glIdx,glIdx]
                if x2[1] <= x1[0]:
                    offs = add_quad_(ln1[0], ln1[0], ln2[1], ln2[2], leaf_points, leaf_polys, offs)
                    offs = add_quad_(ln1[1], ln1[1], ln2[4], ln2[5], leaf_points, leaf_polys, offs)
                    globalIdx_y = globalIdx_y + [glIdx,glIdx]
            elif len(ln1) == 6 and len(ln2) == 2:  # convex to normal case
                x1 = leaf.getLeafVisX(i)
                x2 = leaf.getLeafVisX(i + 1)
                offs = add_quad_(n1, ln1[0], ln2[0], n2, leaf_points, leaf_polys, offs)
                offs = add_quad_(n1, ln1[3], ln2[1], n2, leaf_points, leaf_polys, offs)
                globalIdx_y = globalIdx_y + [glIdx,glIdx]
                if x1[1] <= x2[0]:
                    offs = add_quad_(ln1[1], ln1[2], ln2[0], ln2[0], leaf_points, leaf_polys, offs)
                    offs = add_quad_(ln1[4], ln1[5], ln2[1], ln2[1], leaf_points, leaf_polys, offs)
                    globalIdx_y = globalIdx_y + [glIdx,glIdx]
    return globalIdx_y 

def add_quad_(a, b, c, d, leaf_points, leaf_polys, offs):
    """ used by create_leaf_ """
    q = vtk.vtkPolygon()
    q.GetPointIds().SetNumberOfIds(4)
    for j in range(0, 4):
        q.GetPointIds().SetId(j, offs + j)
    leaf_points.InsertNextPoint(a.x, a.y, a.z)
    leaf_points.InsertNextPoint(b.x, b.y, b.z)
    leaf_points.InsertNextPoint(c.x, c.y, c.z)
    leaf_points.InsertNextPoint(d.x, d.y, d.z)
    leaf_polys.InsertNextCell(q)
    offs += 4
    return offs


def solver_to_polydata(solver, min_, max_, res_):
    """ Creates vtkPolydata from dumux-rosi solver as a structured grid
    @param solver      dumux-rosi solver (RichardsSP, or RichardsNCSP)
    @param min_        minimum of bounding box
    @param max_        maximum of bounding box
    @param res_        resolution, number of cells
    """
    pd = uniform_grid(min_, max_, res_)
    data = solver.getSolutionHead()
    # print("Data range from {:g} to {:g}".format(np.min(data), np.max(data)))
    data_array = vtk_data(data)
    data_array.SetName("pressure head")
    pd.GetCellData().AddArray(data_array)
    return pd


def segs_to_polydata(rs, zoom_factor = 1., param_names = ["age", "radius", "type", "organType" "creationTime"], vals = []):
    """ Creates vtkPolydata from a RootSystem or Plant using vtkLines to represent the root segments
    @param rs             a RootSystem, Plant, or SegmentAnalyser
    @param zoom_factor    a radial zoom factor, since root are sometimes too thin for vizualisation
    @param param_names    parameter names of scalar fields, that are copied to the polydata object
    @return A vtkPolydata object of the root system
    """
    if isinstance(rs, pb.Organism):
        ana = pb.SegmentAnalyser(rs)  # for Organism like Plant or RootSystem
    else:
        ana = rs
    nodes = np_convert(ana.nodes)
    
    
    if len(vals) > 0:
        if isinstance(vals[0], list) or isinstance(vals[0], type(np.array([]))):#len(vals[0]) > 0:
            for i in range(len(vals)) :
                if len(vals[-1-i]) == len(nodes) :
                    ana.addData(param_names[-1-i], vals[-1-i])
        else:
            if len(vals) == len(nodes) :
                ana.addData(param_names[-1], vals)
    
    
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
    c2p = vtk.vtkCellDataToPointData()  # set cell and point data
    c2p.SetPassCellData(True)
    c2p.SetInputData(pd)
    c2p.Update()
    return c2p.GetPolyDataOutput()


def uniform_grid(min_, max_, res):
    """ Creates an uniform grid
    @param min_    minimum of bounding rectangle
    @param max_    maximum of bounding rectangle
    @param res_    cell resolution
    @return A vtkUniformGrid
    """
    grid = vtk.vtkUniformGrid()
    grid.SetDimensions(int(res[0]) + 1, int(res[1]) + 1, int(res[2]) + 1)  # cells to corner points
    grid.SetOrigin(min_[0], min_[1], min_[2])
    s = (max_ - min_) / res
    grid.SetSpacing(s[0], s[1], s[2])
    return grid


def render_window(actor, title, scalarBar, bounds, interactiveImage = True):
    """ puts a vtk actor on the stage (renders an interactive window)

    @param actor                    a (single) actor, or a list of actors (ensemble)
    @param title                    window title
    @param scalarBar                one or a list of vtkScalarBarActor (optional)
    @param bounds                   spatial bounds (to set axes actor, and camera position and focal point)
    @param interactiveImage         make image interactive or static (should be static for google Colab)
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

    #print("render_window, scalarBar", rank)
    if scalarBar:
        if not isinstance(scalarBar, list):
            scalarBar = [scalarBar]
        n = len(scalarBar)
        height = (0.9 / n)
        width = 0.1
        for i, sb in enumerate(scalarBar):
            y = (0.85, 0.1 + i * height)
            sb.SetHeight(0.9 * height)
            sb.SetWidth(width)
            sb.SetPosition(y)
            ren.AddActor2D(sb)

    axes = vtk.vtkAxesActor()
    axes.AxisLabelsOff()  # because i am too lazy to change font size
    translate = vtk.vtkTransform()
    translate.Translate(bounds[0], bounds[2], bounds[4])  # minx, miny, minz
    axes.SetUserTransform(translate)
    ren.AddActor(axes)

    #print("render_window, AddActor", rank)
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

    #print("render_window, Camera", rank)
    # Render Window
    renWin = vtk.vtkRenderWindow()  # boss
    renWin.SetSize(1200, 1000)
    renWin.SetWindowName(title)
    renWin.AddRenderer(ren)
    #print("render_window, interactiveImage?", rank, interactiveImage)

    if interactiveImage:
        iren = vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(renWin)
        #print("render_window, renderA", rank, interactiveImage)
        renWin.Render()
        iren.CreateRepeatingTimer(50)  # [ms] 0.5 s in case a timer event is interested
        iren.AddObserver('KeyPressEvent', lambda obj, ev:keypress_callback_(obj, ev, bounds), 1.0)
        iren.Initialize()  # This allows the interactor to initalize itself. It has to be called before an event loop.
        for a in ren.GetActors():
            a.Modified()  #
        #print("render_window, renderB", rank, interactiveImage)
        renWin.Render()
        return iren
    else:  # necessary?
        renWin.SetOffScreenRendering(1)
        renWin.SetDeviceIndex(0)
        renWin.SetShowWindow(True)
        windowToImageFilter = vtk.vtkWindowToImageFilter()
        windowToImageFilter.SetInput(renWin)
        windowToImageFilter.Update()

        writer = vtk.vtkPNGWriter()
        writer.SetWriteToMemory(1)
        writer.SetInputConnection(windowToImageFilter.GetOutputPort())
        writer.Write()

        # move this somewhere else?
        im = Image(writer.GetResult(), format = "png")
        display(im)


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
        psi = (n - 1) * (float(i) / numberOfColors)  # float(i), float(numberOfColors - i - 1) to reverse
        i0 = np.floor(psi)
        theta = psi - i0
        col = (1 - theta) * np.array(lut_.GetTableValue(int(i0))) + theta * np.array(lut_.GetTableValue(int(i0) + 1))
        lut.SetTableValue(i, col)

    return lut


def create_scalar_bar(lut, grid = None, p_name = "", myRange = None):
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
        if p_name == "organType":  # fix range for for organType
            range = [ 2, 4]
        
        if myRange != None:
            lut.SetTableRange(myRange)
        else:
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


def plot_segments(pd, p_name:str, win_title:str = "", render:bool = True, interactiveImage:bool = True):
    """ see plot roots """
    return plot_roots(pd, p_name, render, interactiveImage)


def plot_roots(pd, p_name:str, win_title:str = "", render:bool = True, interactiveImage:bool = True, 
                myRange = None):
    """ plots the root system
    @param pd                       RootSystem, SegmentAnalyser, or polydata representing the root system (lines, or polylines)
    @param p_name                   parameter name of the data to be visualized
    @param win_title                the windows titles (optionally, defaults to p_name)
    @param render                   render in a new interactive window (default = True)
    @param interactiveImage         make image interactive or static (should be static for google Colab)
    @return a tuple of a vtkActor and the corresponding color bar vtkScalarBarActor
    """
    if isinstance(pd, pb.RootSystem):
        pd = segs_to_polydata(pd, 1., [p_name, "radius", "organType"])

    if isinstance(pd, pb.Plant):
        pd = segs_to_polydata(pd, 1., [p_name, "radius", "organType"])

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
    # mapper.SetArrayName(p_name)
    mapper.SelectColorArray(p_name)
    mapper.UseLookupTableScalarRangeOn()
    plantActor = vtk.vtkActor()
    plantActor.SetMapper(mapper)

    lut = create_lookup_table()  # 24
    scalar_bar = create_scalar_bar(lut, pd, p_name, myRange)  # vtkScalarBarActor
    mapper.SetLookupTable(lut)

    if render:
        ren = render_window(plantActor, win_title, scalar_bar, pd.GetBounds(), interactiveImage)
        if interactiveImage:
            ren.Start()
    return plantActor, scalar_bar, lut


def plot_mesh(grid, p_name, win_title = "", render = True, interactiveImage = True):
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
    mapper.SelectColorArray(p_name)
    mapper.UseLookupTableScalarRangeOn()
    meshActor = vtk.vtkActor()
    meshActor.SetMapper(mapper)
    meshActor.GetProperty().SetRepresentationToWireframe();

    lut = create_lookup_table()
    scalar_bar = create_scalar_bar(lut, grid, p_name)  # vtkScalarBarActor
    mapper.SetLookupTable(lut)

    if render:
        ren = render_window(meshActor, win_title, scalar_bar, grid.GetBounds(), interactiveImage)
        if interactiveImage:
            ren.Start()
    return [meshActor], scalar_bar


def plot_mesh_cuts(grid, p_name, nz = 3, win_title = "", render = True, interactiveImage = True,
                    myRange = None):
    """ plots orthogonal nz vertical cuts z[:-1] (xy-planes), with z = linspace(min_z, max_z, nz+1),
    and two additonal sclices at x=0 (yz-plane), y=0 (xz-plane)
    @param grid         some vtk grid (structured or unstructured)
    @param p_name       parameter to visualize
    @param nz           number of vertical slices
    @param win_title    the windows titles (optionally, defaults to p_name)
    @param render       render in a new interactive window (default = True)
    @return a tuple of a list of vtkActors and a single corresponding color bar vtkScalarBarActor
    """
    #print("plot_mesh_cuts render_begin", rank)
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
    scalar_bar = create_scalar_bar(lut, grid, p_name, myRange)

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
        # m.SetArrayName(p_name)
        m.SelectColorArray(p_name)
        m.UseLookupTableScalarRangeOn()
        m.SetLookupTable(lut)
        m.SetColorModeToMapScalars();
        a = vtk.vtkActor()  # create plane actor
#         a.GetProperty().SetColor(1.0, 1, 0)
#         a.GetProperty().SetLineWidth(2)
        a.SetMapper(m)
        actors.append(a)
    #print("plot_mesh_cuts render", rank)
    if render:
        ren = render_window(actors, win_title, scalar_bar, grid.GetBounds(), interactiveImage)
        if interactiveImage:
            ren.Start()

    return actors, scalar_bar

def plot_roots_and_soil(rs, pname:str, rp, s, periodic:bool, min_b, max_b, cell_number, filename:str = "", 
        sol_ind = 0, interactiveImage = True, extraArray = [], extraArrayName = "",
        oneScalarBar = False):
    """ Plots soil slices and roots, additionally saves both grids as files
    @param rs            some Organism (e.g. RootSystem, MappedRootSystem, ...) or MappedSegments
    @param pname         root and soil parameter that will be visualized ("pressure head", or "water content")
    @param s             soil, of type RichardsSP, or RichardsNCSP
    @param rp            root parameter segment data (will be added, in case SegmentAnalyser is creaeted)
    @param periodic      if yes the root system will be mapped into the domain
    @param min_b         minimum of domain boundaries
    @param max_b         maximum of domain boundaries
    @param cell_number   domain resolution
    @param filename      file name (without extension)
    """
    if rank == 0:
        if isinstance(rs, pb.SegmentAnalyser):
            ana = rs
        else:
            ana = pb.SegmentAnalyser(rs)
            if rank == 0:
                ana.addData(pname, rp)
        if periodic:
            w = np.array(max_b) - np.array(min_b)
            ana.mapPeriodic(w[0], w[1])
        pd = segs_to_polydata(ana, 1., [pname, "radius"])

        pname_mesh = "pressure head"  # pname <------ TODO find better function arguments
    soil_grid = uniform_grid(np.array(min_b), np.array(max_b), np.array(cell_number))
    soil_water_content = vtk_data(list(np.array(s.getWaterContent()).reshape(-1,1)))
    

    soil_water_content.SetName("water content")
    
    
    soil_grid.GetCellData().AddArray(soil_water_content)
    path = s.results_dir#"results/"
    soil_pressure = vtk_data(np.array(s.getSolutionHead()))
    soil_pressure.SetName("pressure head")  # in macroscopic soil
    soil_grid.GetCellData().AddArray(soil_pressure)
    if sol_ind > 0:
        d = vtk_data(np.array(s.getSolution(sol_ind)))
        pname_mesh = "solute" + str(sol_ind)
        d.SetName(pname_mesh)  # in macroscopic soil
        soil_grid.GetCellData().AddArray(d)
    if sol_ind < 0:
        d = vtk_data(extraArray)
        pname_mesh = extraArrayName
        d.SetName(pname_mesh)  # in macroscopic soil
        soil_grid.GetCellData().AddArray(d)
    
    if rank == 0:      
        if oneScalarBar:
            rangeS = [[], []]
            p_nameS = [pname,pname_mesh]
            for gid, grid in enumerate([pd, soil_grid]):   
                barRange = [0, 1]
                a = grid.GetCellData().GetAbstractArray(p_nameS[gid])
                if a:
                    barRange = a.GetRange()
                    grid.GetCellData().SetActiveScalars(p_nameS[gid])
                else:
                    a = grid.GetPointData().GetAbstractArray(p_nameS[gid])
                    grid.GetPointData().SetActiveScalars(p_nameS[gid])
                    if a:
                        barRange = a.GetRange()
                if p_nameS[gid] == "organType":  # fix barRange for for organType
                    barRange = [ 2, 4]
                rangeS[gid] = barRange
            sameRange = [min(rangeS[0][0], rangeS[1][0]), max(rangeS[0][1], rangeS[1][1])]
        else:
            sameRange = None
            
        rootActor, rootCBar, lut_ = plot_roots(pd, pname, "", False, myRange = sameRange)
        meshActors, meshCBar = plot_mesh_cuts(soil_grid, pname_mesh, 7, "", False, myRange = sameRange)
        meshActors.extend([rootActor])
        #ren = render_window(meshActors, filename, [meshCBar, rootCBar], pd.GetBounds(), interactiveImage)
        #if interactiveImage:
        #    ren.Start()

        if filename:
            path = s.results_dir#"results/"
            write_vtp(path + filename + ".vtp", pd)
            write_vtu(path + filename + ".vti", soil_grid)


def plot_roots_and_mesh(rs, pname_root, mesh, pname_mesh, periodic:bool, xx = 1, yy = 1, filename:str = "", interactiveImage = True):
    """ Plots soil slices and roots, additionally saves both grids as files
    @param rs            some Organism (e.g. RootSystem, MappedRootSystem, ...) or MappedSegments
    @param pname_root    root parameter that will be visualized 
    @param mesh          vtk grid
    @param pname_mesh    name of grid cell data array that will be visualized 
    @param periodic      if yes the root system will be mapped into the domain
    @param xx, yy        witdh and height of periodic domain    
    @param filename      file name (without extension)
    
    How to add data to a vtk grid: 
    
    celldata = vtk.vtkDoubleArray()
    celldata.SetName("data")
    celldata.SetNumberOfValues(n)
    for j in range(0, n):
        celldata.SetValue(j, data[j])
    grid.GetCellData().AddArray(celldata)    
    """
    ana = pb.SegmentAnalyser(rs)
    if periodic:
        ana.mapPeriodic(xx, yy)
    pd = segs_to_polydata(ana, 1., ["radius", "subType", "creationTime", pname_root])
    rootActor, rootCBar = plot_roots(pd, pname_root, "", False)
    meshActors, meshCBar = plot_mesh_cuts(mesh, pname_mesh, 7, "", False)
    meshActors.extend([rootActor])
    ren = render_window(meshActors, filename, [meshCBar, rootCBar], pd.GetBounds(), interactiveImage)
    if interactiveImage:
        ren.Start()
    if filename:
        path = "results/"
        write_vtp(path + filename + ".vtp", pd)
        write_vtu(path + filename + ".vtu", soil_grid)


def write_soil(filename, s, min_b, max_b, cell_number, solutes = []):
    """ TODO """
    soil_grid = uniform_grid(np.array(min_b), np.array(max_b), np.array(cell_number))
    soil_water_content = vtk_data(np.array(s.getWaterContent()))
    soil_water_content.SetName("water content")
    soil_grid.GetCellData().AddArray(soil_water_content)
    soil_pressure = vtk_data(np.array(s.getSolutionHead()))
    soil_pressure.SetName("pressure head")  # in macroscopic soil
    soil_grid.GetCellData().AddArray(soil_pressure)
    for i, s_ in enumerate(solutes):
        d = vtk_data(np.array(s.getSolution(i + 1)))
        d.SetName(s_)  # in macroscopic soil
        soil_grid.GetCellData().AddArray(d)
    if rank == 0:
        write_vtu(filename + ".vtu", soil_grid)


def write_plant(filename, plant, add_params = []):
    """ write the plants organ ceneterlines and leafs into two seperate vtp files"""
    params = ["radius", "subType", "organType", "age"]
    params.extend(add_params)
    pd = segs_to_polydata(plant, 1., params)
    write_vtp(filename + ".vtp", pd)

    leaf_points = vtk.vtkPoints()
    leaf_polys = vtk.vtkCellArray()  # describing the leaf surface area
    leafes = plant.getOrgans(pb.leaf)
    for l in leafes:
        create_leaf_(l, leaf_points, leaf_polys)

    pd_leafs = vtk.vtkPolyData()
    pd_leafs.SetPoints(leaf_points)
    pd_leafs.SetPolys(leaf_polys)
    write_vtp(filename + "_leafs.vtp", pd_leafs)


def plot_roots_and_soil_files(filename: str, pname:str, interactiveImage = True):
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
    ren = render_window(meshActors, filename, meshCBar, soil_grid.GetBounds(), interactiveImage)
    if interactiveImage:
        ren.Start()


class AnimateRoots:
    """ class to make an interactive animation (TODO unfinished and doc)"""

    def __init__(self, rootsystem = None):
        self.rootsystem = rootsystem
        self.root_name = "subType"
        # self.soil_name = "subType"
        #
        self.min = None
        self.max = None
        self.res = None
        self.soil_data = True  # soil data
        self.soil = None
        self.cuts = False  # Wireframe, or cuts
        self.plant = False  # use plot_roots or plot_plant
        #
        self.actors = []
        self.iren = None
        self.color_bar = None
        self.bounds = None
        self.avi_name = None
        self.fram_c = 0

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
        if self.avi_name:
            write_png(renWin, self.avi_name + str(self.fram_c))
            print("saved", self.avi_name + str(self.fram_c) + ".png")
            self.fram_c = self.fram_c + 1

    def create_root_actors(self):
        if self.rootsystem:
           pd = segs_to_polydata(self.rootsystem, 1., [self.root_name, "radius"])

           if self.plant:
               newRootActor, rootCBar = plot_plant(self.rootsystem, self.root_name, False)
           else:
               newRootActor, rootCBar = plot_roots(pd, self.root_name, "", False)
           if isinstance(newRootActor, list):
               for a in newRootActor:
                   self.actors.append(a)
           else:
               self.actors.append(newRootActor)
           self.color_bar = rootCBar
           self.bounds = pd.GetBounds()

    def create_soil_actors(self):
        if self.soil_data:
            if self.cuts:
                # meshActor, meshCBar = plot_mesh_cuts(grid, p_name, nz = 3, win_title = "", False):
                meshActor, meshCBar, grid = plot_mesh_yz(self.soil, "pressure head", self.min, self.max, self.res)
                self.color_bar = meshCBar
            else:
                grid = uniform_grid(self.min, self.max, self.res)
                meshActor, meshCBar = plot_mesh(grid, "", "", False)

            self.actors.extend(meshActor)
            self.bounds = grid.GetBounds()

