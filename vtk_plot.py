import plantbox as pb
from vtk_tools import *

import time
import numpy as np
import vtk
import os
import datetime
from vtk.util import numpy_support

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


def plot_plant(plant, p_name, render = True, printout = False, outputDirectory = "output/", timestamp = False, date = "", sim_name = "", GraphicalAccuracy = False, ExtraParam = None, NormalsZValue = [],LeafSurfaceList = [],Thickness = False):
    """
        plots a whole plant as a tube plot, and additionally plot leaf surface areas as polygons 
        @param plant: requires a "plantbox.organism"
        
        @param p_name: parameter desired to show on the graph and color bar.
        
        @param printout: This is the main toggle for the output of the 3D structure of the plant. If not set to True, the following parameters will not have an effect(outputDirectory, timestamp, date, sim_name, GraphicalAccuracy).
        
        @param outputDirectory: sets the path to store the file created by the printout parameter (LeafOBJ and SegsOBJ)
        
        @param timestamp: If set to True, will create permanent copies of the files created by printout in a subfolder called "stored"(it cannot be changed). The "stored" folder will be created bellow the outputDirectory.
        
        @param date: If left blank, the date of the current day and time will be created and become the subfolder name if "timestamp" is set to True.
        If filled in, the name of the subfile will change accordingly. It is possible to create subfolders if filled correctly (ex.: "test/test/", that will create the following folders: output/test/test/. The files will be written in the last folder created).
        
        @param sim_name: If left blank, does nothing. If filled in, it will create a sub folder before the date (output/sim_name/date).
        
        Directory order: outputDirectory/stored/sim_name/date/ 3D files (LeafOBJ and SegsOBJ).
        
        @param GraphicalAccuracy: Allows the code to rerun the creation of the leaves so they look better in the render window when printout is set to True.
        
        @param ExtraParam: Allows the user to add multiple data sets to the plant. The first slot must be the name, the 2nd the stem values (it needs to be filled in), the 3rd is for the leaves and is optionnal (the leaves will look funky if left empty). 
        ex.: ExtraParam =[["Test", Test_stem, Test_leaf], ["Test2", Test2_stem, Test2_leaves]] with Test_stem and Test_leaves matrices of the right sizes (stem: number of segments, leaves: number of quadrilaterals).
        @param Thickness: Creates thickened leaves. This parameter REQUIRES printout to be set to "True" to have any impact. This parameter is not compatible with ExtraParam for the moment.
        
    """
    
    if printout:
        # Checks if the output directory name ends with a "/"
        if not(outputDirectory.endswith("/")):
            outputDirectory = outputDirectory + "/"
        #Create the output folder, if not present
        if not os.path.exists(outputDirectory):
            os.makedirs(outputDirectory)
        
        #Clear the content of the folder of old modelisations in the "output/" folder. ATTENTION, it clears ALL ".obj" file in it.
        #This is made to prevent empty files from being kept from old simulations with a bigger number of leaves or 
        #segments per leaf.
        filelist = [ f for f in os.listdir(outputDirectory) if f.endswith(".obj") ]
        for f in filelist:
            os.remove(os.path.join(outputDirectory, f))
            
        #Creates the date if needed
        if (date == "") and (timestamp == True):            
            x = datetime.datetime.now()
            year = str(x.year)
            month = str(x.month)
            day = str(x.day)
            hour = str(x.hour)
            minute = str(x.minute)
            second = str(x.second) 
            if len(hour) == 1:
                hour = "0" + hour
            if len(minute) == 1:
                minute = "0" + minute
            if len(second) == 1:
                second = "0" + second
            date = year + "-" + month + "-" + day +"_" + hour + "h" + minute + "min"+ second + "s"
    

    
    
    # plant as tube plot
    pd = segs_to_polydata(plant, 1., ["radius", "organType", "creationTime", p_name], printout = printout, outputDirectory = outputDirectory, timestamp = timestamp, date = date , sim_name = sim_name, ExtraParam = ExtraParam )  # poly data
    # global lut, tube_plot_actor, color_bar
    tube_plot_actor, color_bar = plot_roots(pd, p_name, "", render = False)
    lut = None
    if ExtraParam != None:
        for P in ExtraParam:
            if P[0] == p_name:
                if P[2] != []:
                    # print("range", min(P[2]), max(P[2]))
                    tube_plot_actor, color_bar, lut = plot_roots(pd, p_name, "", render = False, ReturnLut = True, RangeExtra = [min(P[2]), max(P[2])])
    

    # leafes as polygons
    leaf_points = vtk.vtkPoints()
    leaf_polys = vtk.vtkCellArray()  # describing the leaf surface area
    
    leafes = plant.getOrgans(pb.leaf)
        
    
    if GraphicalAccuracy and printout and render:
        UIVal = 2
    else:
        UIVal = 1
        
    LinearExtrusion = []
    for i in range(UIVal): #This will run TWICE the creation of the leaf (green quadrilaterals) if set to 2. The first step creates individuals segments the second creates all the segments appended in the "polydata".
        global LeafSegIDListList
        LeafSegIDListList = []
        global QuadCounterListList
        QuadCounterListList = []
        # print("Leaf iteration",i) #Let's you know at which run you are at (If printout set to "True", the leafes won't appear on the render_window, it needs a second run to redraw them correctly.)
        if (i == 0) and printout: #doesn't overwrite printout if set to False at the beginning
            printout = True 
        elif i==1:
            printout = False
        
        for k, l in enumerate(leafes):
            LeafSegIDList = []
            QuadCounterList = []
            if Thickness and printout:
                LinearExtrusion = create_leaf(l, leaf_points, leaf_polys, LeafSegIDList = LeafSegIDList, QuadCounterList = QuadCounterList, k = k, printout = printout, outputDirectory = outputDirectory, timestamp = timestamp, date = date, sim_name = sim_name, NormalsZValue = NormalsZValue, LeafSurfaceList = LeafSurfaceList, Thickness = Thickness, OutPut = LinearExtrusion)
            else:
                create_leaf(l, leaf_points, leaf_polys, LeafSegIDList = LeafSegIDList, QuadCounterList = QuadCounterList, k = k, printout = printout, outputDirectory = outputDirectory, timestamp = timestamp, date = date, sim_name = sim_name, NormalsZValue = NormalsZValue, LeafSurfaceList = LeafSurfaceList)
            
            #Debug:
            # print("Current leaf radius", l.getParameter("radius"))
            # print(k)#Current leaf number
            # print("LeafSegIDList:",LeafSegIDList,"Number of Quadrilaterals added per segment:",QuadCounterList)
            # print("Length of the above lists (both should be the same amount):", len(LeafSegIDList),";",len(QuadCounterList)) #The lengths should be the same.
            LeafSegIDListList.extend(LeafSegIDList)
            QuadCounterListList.extend(QuadCounterList)
        if i==1:
            printout = True
    # print("LeafSegIDListList, size:", len(LeafSegIDListList), "IDs:", LeafSegIDListList)
    # print("QuadCounterListList, size:", len(QuadCounterListList), "IDs:", QuadCounterListList)
    # print(LeafSurfaceList) #Computes the surface for each quadrilateral forming the each segment of each leaf
    print("Number of extruded leafes' segments (should be higher then 0):", len(LinearExtrusion))

    polyData = vtk.vtkPolyData()
    polyData.SetPoints(leaf_points)
    polyData.SetPolys(leaf_polys)

    vtkAppendPolyData = vtk.vtkAppendPolyData()
    if Thickness:
        # vtkAppendPolyData.FastDelete()
        for elem in LinearExtrusion:
            vtkAppendPolyData.AddInputData(elem.GetOutput())
            vtkAppendPolyData.Update()
    
    
    ThicknessPolyData = []
    if ExtraParam != None:
        for P in ExtraParam:
            data1 = vtk.vtkDataArray.CreateDataArray(vtk.VTK_DOUBLE)
            data1.SetNumberOfComponents(1)  # number of components
            data1.SetNumberOfTuples(1)
            # print("params length", len(P[1]))
            for d in range(len(P[2])):
                if Thickness:
                    # float(P[2][d])
                    data1.InsertTuple1(d*2+1, float(P[2][d])) #upper surface of the leaf
                    data1.InsertTuple1(d*2, float(P[2][d]))   #lower surface 
                    for j in range (4):
                        data1.InsertTuple1(d*4 + j + len(P[2])*2, float(P[2][d]))
                else:
                    data1.InsertTuple1(d, float(P[2][d]))
                # data1.InsertTuple1(d, float(P[2][d]))
                
            data1.SetName(P[0])
            polyData.GetCellData().AddArray(data1)
            if Thickness:
                ThicknessPolyData = vtkAppendPolyData.GetOutput()
                ThicknessPolyData.GetCellData().AddArray(data1)
    

    mapper = vtk.vtkPolyDataMapper()
    if Thickness and ExtraParam != None:
        mapper.SetInputData(ThicknessPolyData)
    elif Thickness and not(ExtraParam != None):
        mapper.SetInputData(vtkAppendPolyData.GetOutput())
    else:
        mapper.SetInputData(polyData)
    
    if ExtraParam != None:
        for P in ExtraParam:
            if P[0] == p_name:
                if P[2] != []:
                    mapper.ScalarVisibilityOn();
                    mapper.SetScalarModeToUseCellFieldData()  # maybe because radius is active scalar in point data?
                    # mapper.SetArrayName(p_name)
                    mapper.SelectColorArray(p_name)
                    mapper.UseLookupTableScalarRangeOn()
                    mapper.SetLookupTable(lut)
    
    actor = vtk.vtkActor()
    actor.SetMapper(mapper);
    
    colors = vtk.vtkNamedColors()
    if ExtraParam == None:
        actor.GetProperty().SetColor(colors.GetColor3d("Green"))
        
    if render:
        # render_window([tube_plot_actor], "plot_plant", color_bar, tube_plot_actor.GetBounds()).Start()      #shows only plant
        render_window([tube_plot_actor, actor], "plot_plant", color_bar, tube_plot_actor.GetBounds()).Start() #shows plant and leaves
    return [tube_plot_actor, actor], color_bar


def create_leaf(leaf, leaf_points, leaf_polys, meshing = False, LeafSegIDList = [], QuadCounterList = [],k = 0, printout = False, outputDirectory = "output/", timestamp = False, date = "", sim_name = "", NormalsZValue = [], LeafSurfaceList = [], Thickness = False, OutPut = []):
    
    offs = leaf_points.GetNumberOfPoints()
    if printout:
        if not(outputDirectory.endswith("/")):
                outputDirectory = outputDirectory + "/"
        if not(sim_name.endswith("/")) and (sim_name != ""): #if the simulation name is empty, it doesn't add a "/"
                sim_name = sim_name + "/"
    

    for i in range(0, leaf.getNumberOfNodes() - 1):  
        linearExtrusion1 = vtk.vtkLinearExtrusionFilter()
        if printout:
            offs = 0
            leaf_points.Reset()
            leaf_polys.Reset()
            
        
        
        QuadCounter = 0
        ln1 = leaf.getLeafVis(i)
        ln2 = leaf.getLeafVis(i + 1)
        IdSeg = (leaf.getNodeId(i+1)-1)  #allows the retrieval of the relative node IDs from the organ.
        LeafSegIDList.append(IdSeg)      #allows the creation of a list of the IDs of the used segment in leaf creation.

        if len(ln1) > 0 or len(ln2) > 0:
            n1 = leaf.getNode(i)
            n2 = leaf.getNode(i + 1)  #exports the 3D position of the node

#             if len(ln1) > 0 and len(ln2) == 0:
#                 print(" 2 -> 0")
#             if len(ln1) == 0 and len(ln2) > 0:
#                 print(" 0 -> 2")

            if len(ln1) == 2 and len(ln2) == 2:  # normal case
                offs = add_quad_(n1, ln1[0], ln2[0], n2, leaf_points, leaf_polys, offs)
                offs = add_quad_(n1, ln1[1], ln2[1], n2, leaf_points, leaf_polys, offs)
                QuadCounter = QuadCounter + 2 
            elif len(ln1) == 6 and len(ln2) == 6:  # convex case
                offs = add_quad_(n1, ln1[0], ln2[0], n2, leaf_points, leaf_polys, offs)
                offs = add_quad_(ln1[1], ln1[2], ln2[2], ln2[1], leaf_points, leaf_polys, offs)
                offs = add_quad_(n1, ln1[3], ln2[3], n2, leaf_points, leaf_polys, offs)
                offs = add_quad_(ln1[4], ln1[5], ln2[5], ln2[4], leaf_points, leaf_polys, offs)
                QuadCounter = QuadCounter + 4
            elif len(ln1) == 2 and len(ln2) == 6:  # normal to convex case
                x1 = leaf.getLeafVisX(i)
                x2 = leaf.getLeafVisX(i + 1)
                offs = add_quad_(n1, ln1[0], ln2[0], n2, leaf_points, leaf_polys, offs)
                offs = add_quad_(n1, ln1[1], ln2[3], n2, leaf_points, leaf_polys, offs)
                QuadCounter = QuadCounter + 2
                if x2[1] <= x1[0]:
                    offs = add_quad_(ln1[0], ln1[0], ln2[1], ln2[2], leaf_points, leaf_polys, offs)
                    offs = add_quad_(ln1[1], ln1[1], ln2[4], ln2[5], leaf_points, leaf_polys, offs)
                    QuadCounter = QuadCounter + 2
            elif len(ln1) == 6 and len(ln2) == 2:  # convex to normal case
                x1 = leaf.getLeafVisX(i)
                x2 = leaf.getLeafVisX(i + 1)
                offs = add_quad_(n1, ln1[0], ln2[0], n2, leaf_points, leaf_polys, offs)
                offs = add_quad_(n1, ln1[3], ln2[1], n2, leaf_points, leaf_polys, offs)
                QuadCounter = QuadCounter + 2
                if x1[1] <= x2[0]:
                    offs = add_quad_(ln1[1], ln1[2], ln2[0], ln2[0], leaf_points, leaf_polys, offs)
                    offs = add_quad_(ln1[4], ln1[5], ln2[1], ln2[1], leaf_points, leaf_polys, offs)
                    QuadCounter = QuadCounter + 2
            QuadCounterList.append(QuadCounter)
            
            # print(dir(leaf_points),"leaf_points")
            # numpy_support.vtk_to_numpy(polyData.GetPoints().GetData()) #Allows to show the inside of the points position
            
            #printing the positions, before using them. If uncommented, breaks the visual for some reason.
            # for i in range(len(ln1)):
            #     print("ln1",ln1[i],"ln2",ln2[i])
            if printout:
                polyData = vtk.vtkPolyData()
                polyData.SetPoints(leaf_points)
                polyData.SetPolys(leaf_polys)
                
                #Surface computation
                qualityFilter = vtk.vtkMeshQuality()
                qualityFilter.SetInputData(polyData) 
                # qualityFilter.SetTriangleQualityMeasureToArea()
                qualityFilter.SetQuadQualityMeasureToArea()
                qualityFilter.Update()
                # print("t",np.array(qualityFilter.GetOutput().GetCellData().GetArray("Quality"))) #output of vector with the same amount of inputs as the number of quadrilateral in the segment.
                LeafSurfaceList.append(np.array(qualityFilter.GetOutput().GetCellData().GetArray("Quality")))
                # print(LeafSurfaceList)
                
                
                #Normal to the leaf surface
                normal = vtk.vtkPolyDataNormals()
                normal.AddInputData(polyData)
                normal.ComputePointNormalsOn()
                normal.SplittingOn()
                normal.Update()
                NormalData = normal.GetOutput().GetPointData().GetNormals()
                NormalVector = np.array(NormalData)
                # print(type(NormalVector),NormalVector.size,NormalVector) #debug
                # print(NormalVector.size)
                
                # for t in range(0,len(NormalVector),4): #each surface is defined by for points thus 4 normals are computed per surface. We jump by a step of 4 to go to the next quadrilateral. If the leaf has a hand like shape, it could have more than 2 quadrilaterals per segment (due to the dents).
                #     NormalsZValue.append(abs(NormalVector[t][2])) # adds a value for each quadrilateral of the segment
                # print(1-abs(NormalVector[t][0]))
                if NormalVector.size > 1:
                    NormalsZValue.append(abs(NormalVector[0][2]))# Extract the Z value [2] and only adds the first quadrilateral value for the whole segment[0].
                elif NormalVector.size == 1: #happens when the surface is a line and thus no normal can be computed
                    NormalsZValue.append(float(0))

                
                writer = vtk.vtkOBJWriter()
                
                
                if not os.path.exists(outputDirectory):
                    os.makedirs(outputDirectory)
                writer.SetFileName(outputDirectory + "LeafOBJ_{}_{}".format(k,i) + ".obj")
                writer.SetInputData(polyData)
                writer.Update()
                writer.Write()
                #If the timestamp is active, the output files will printout twice. Once in the output with the same format
                # and the second time in the "stored" folder with the date-time.
                if timestamp:
                    if not(sim_name.endswith("/")) and (sim_name != ""): #if the simulation name is empty, it doesn't add a "/"
                        sim_name = sim_name + "/"
                    if not(date.endswith("/")):
                        date = date + "/"
                    path = outputDirectory + "stored/" + sim_name + date + "LeafOBJ/"
                    if not os.path.exists(path):
                        os.makedirs(path)
                    writer.SetFileName(path + "LeafOBJ_{}_{}".format(k,i) + ".obj")
                    writer.SetInputData(polyData)
                    writer.Update()
                    writer.Write()
                    
                #Possible creation of thick leaves.  
                # linearExtrusion1 = vtk.vtkLinearExtrusionFilter()
                
                linearExtrusion1.AddInputData(polyData)
                linearExtrusion1.SetScaleFactor(leaf.getParameter("radius")) #Set the thickness of the leave based on the radius of the segment corresponding
                linearExtrusion1.SetExtrusionTypeToNormalExtrusion()
                linearExtrusion1.SetVector(NormalVector[0]) #Sets the vector. Since one segment has the same orientation, the first quadrilateral's normal is used
                OutPut.append(linearExtrusion1)
                # print(linearExtrusion1.GetScaleFactor())
                
                # print(dir(linearExtrusion1))  #shows in the console all the functions available in this object
                
                writer = vtk.vtkOBJWriter()
                # print(dir(writer)) #Debug
                writer.SetFileName(outputDirectory + "LeafThicknessOBJ{}_{}".format(k,i) + ".obj")
                writer.SetInputConnection(linearExtrusion1.GetOutputPort())
                writer.Write()
                if timestamp:
                    if not(sim_name.endswith("/")) and (sim_name != ""): #if the simulation name is empty, it doesn't add a "/"
                        sim_name = sim_name + "/"
                    if not(date.endswith("/")):
                        date = date + "/"
                    path = outputDirectory + "stored/" + sim_name + date + "LeafThicknessOBJ/"
                    if not os.path.exists(path):
                        os.makedirs(path)
                    writer.SetFileName(path + "LeafThicknessOBJ{}_{}".format(k,i) + ".obj")
                    writer.SetInputConnection(linearExtrusion1.GetOutputPort())
                    writer.Update()
                    writer.Write()
    #Debug:   
    # print("Leaf number:",k,"Number of segments:",i+1)
    if Thickness == True and printout == True:
        return OutPut


def add_quad_(a, b, c, d, leaf_points, leaf_polys, offs, meshing = False):
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
    if meshing == False:
        return offs
    elif meshing == True:
        return offs, mesh
    else:
        return offs


def solver_to_polydata(solver, min_, max_, res_):
    """ Creates vtkPolydata from dumux-rosi solver as a structured grid
    @param solver
    @param min_ 
    @param max_ 
    @param res_
    """
    pd = uniform_grid(min_, max_, res_)
    data = solver.getSolutionHead()
    # print("Data range from {:g} to {:g}".format(np.min(data), np.max(data)))
    data_array = vtk_data(data)
    data_array.SetName("pressure head")
    pd.GetCellData().AddArray(data_array)
    return pd


def segs_to_polydata(rs, zoom_factor = 1., param_names = ["age", "radius", "type", "organType" "creationTime"], printout = False, outputDirectory = "output/", timestamp = False, date = "", sim_name = "", ExtraParam = None ):
    """ Creates vtkPolydata from a RootSystem or Plant using vtkLines to represent the root segments 
    @param rs             a RootSystem, Plant, or SegmentAnalyser
    @param zoom_factor    a radial zoom factor, since root are sometimes too thin for vizualisation
    @param param_names    parameter names of scalar fields, that are copied to the polydata object   
    @return A vtkPolydata object of the root system
    """
    if isinstance(rs, pb.Organism):
        ana = pb.SegmentAnalyser(rs)  # for Organism like Plant or RootSystem
        if ExtraParam != None:
            for P in ExtraParam:
                ana.addData(P[0], P[1])
    else:
        ana = rs
    # print("param_names", param_names) #it shows that only the required name is added to the polydata (saves data storage).
    nodes = np_convert(ana.nodes)
    segs = np_convert(ana.segments)        
    points = vtk_points(nodes)
    cells = vtk_cells(segs)
    pd = vtk.vtkPolyData()
    pd.SetPoints(points)
    pd.SetLines(cells)  # SetPolys not working
    for n in param_names:
        param = np.array(ana.getParameter(n))
        if param.shape[0] == segs.shape[0]:
            if n == "radius":
                param *= zoom_factor
                # print(param)
            data = vtk_data(param)
            data.SetName(n)
            pd.GetCellData().AddArray(data)
        else:
            print("segs_to_polydata: Warning parameter " + n + " is skipped because of wrong size", param.shape[0], "instead of", segs.shape[0])
    c2p = vtk.vtkCellDataToPointData()  # set cell and point data
    c2p.SetPassCellData(True)
    c2p.SetInputData(pd)
    c2p.Update()
    
#     if printout: #this will duplicate the work but maintains the ability to render the structure.
#         if not(outputDirectory.endswith("/")):
#             outputDirectory = outputDirectory + "/"
#         #Create the output folder, if not present
#         if not os.path.exists(outputDirectory):
#             os.makedirs(outputDirectory)
#         # Checks if the simulation name ends with a "/"
#         if not(sim_name.endswith("/")) and (sim_name != ""): #if the simulation name is empty, it doesn't add a "/"
#             sim_name = sim_name + "/"
        
#         for i, seg in enumerate(segs):
#             point = vtk_points(np.array([nodes[seg[0]],nodes[seg[1]]]))
#             cell = vtk_cells(np.array([[0,1]]))
#             radius = np.array(ana.getParameter("radius"))[i]
#             radius *= zoom_factor
#             # print(radius)

#             ##Code from vtk_tools.py -> vtk_data. vtk_data can't work with size "1" arrays as it is not read as an array but as a singular item of its type.
#             data1 = vtk.vtkDataArray.CreateDataArray(vtk.VTK_DOUBLE)
#             data1.SetNumberOfComponents(1)  # number of components
#             data1.SetNumberOfTuples(1)
#             data1.InsertTuple1(0, radius)
#             data1.SetName("radius")
#             ##

#             pd1 = vtk.vtkPolyData()
#             pd1.SetPoints(point)
#             pd1.SetLines(cell)
#             pd1.GetCellData().AddArray(data1)

#             c2p1 = vtk.vtkCellDataToPointData()  # set cell and point data
#             c2p1.SetPassCellData(True)
#             c2p1.SetInputData(pd1)
#             c2p1.Update()
#             pd2 = c2p1.GetPolyDataOutput()
#             pd2.GetPointData().SetActiveScalars("radius")

#             tubeFilter = vtk.vtkTubeFilter()
#             tubeFilter.SetInputData(pd2)
#             tubeFilter.SetNumberOfSides(9)
#             tubeFilter.SetVaryRadiusToVaryRadiusByAbsoluteScalar()
#             tubeFilter.CappingOn()
#             tubeFilter.Update()

#             writer = vtk.vtkOBJWriter()
#             writer.SetFileName(outputDirectory + "SegsOBJ_{}".format(str(i)) + ".obj")
#             writer.SetInputConnection(tubeFilter.GetOutputPort())
#             writer.Write()
#             if timestamp:
#                 if not(date.endswith("/")):
#                     date = date + "/"
#                 path = outputDirectory + "stored/" + sim_name + date + "SegsOBJ/"
#                 if not os.path.exists(path):
#                     os.makedirs(path)
#                 writer.SetFileName(path + "SegsOBJ_{}".format(str(i)) + ".obj")
#                 writer.SetInputConnection(tubeFilter.GetOutputPort())
#                 writer.Update()
#                 writer.Write()
            

#         if ("Oeuf" in sim_name)  and ("Pâques" in date): #Oeuf de Pâques(French)
#             import webbrowser
#             webbrowser.open('https://www.youtube.com/watch?v=dQw4w9WgXcQ')
            
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


def render_window(actor, title, scalarBar, bounds):
    """ puts a vtk actor on the stage (renders an interactive window)
    
    @param actor                    a (single) actor, or a list of actors (ensemble)
    @param title                    window title 
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
    renWin.SetSize(1200, 1000)
    renWin.SetWindowName(title)
    renWin.AddRenderer(ren)

    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)
    renWin.Render()
    iren.CreateRepeatingTimer(50)  # [ms] 0.5 s in case a timer event is interested
    iren.AddObserver('KeyPressEvent', lambda obj, ev:keypress_callback_(obj, ev, bounds), 1.0)
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
    if isinstance( tableIdx, int): 
        colorSeries.SetColorScheme(tableIdx)
    if isinstance (tableIdx, str):
        test = "vtk.vtkColorSeries." + tableIdx #"BREWER_SEQUENTIAL_YELLOW_ORANGE_BROWN_9"#"BREWER_SEQUENTIAL_BLUE_GREEN_9"
        # print(eval(test))
        colorSeries.SetColorScheme(eval(test)) #transform the string into the command
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


def create_scalar_bar(lut, grid = None, p_name = "", RangeOverwrite = None):
    """ creates a vtkScalarBarActor, for a vtkLookupTable, sets the active scalar to p_name
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
        if p_name == "organType":  # fix range for organType
            range = [ 2, 4]
        if RangeOverwrite != None:
            range = RangeOverwrite
        # range = [0, 150]  #bypassing for testing
        lut.SetTableRange(range)
    # print("range", range, type(range))

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
    # scalarBar.SetNumberOfLabels(3) #by default, it's 5

    return scalarBar


def plot_roots(pd, p_name:str, win_title:str = "", render:bool = True, ReturnLut = False, RangeExtra = None):
    """ plots the root system 
    @param pd         RootSystem, SegmentAnalyser, or polydata representing the root system (lines, or polylines)
    @param p_name     parameter name of the data to be visualized
    @param win_title  the windows titles (optionally, defaults to p_name)
    @param render     render in a new interactive window (default = True)
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
    # write_vtp("StructureVTP.vtp",pd)
    
    tubeFilter = vtk.vtkTubeFilter()
    tubeFilter.SetInputData(pd)
    tubeFilter.SetNumberOfSides(9)
    tubeFilter.SetVaryRadiusToVaryRadiusByAbsoluteScalar()
    tubeFilter.Update()
    
    #extract the entire surface of the segments of the plant. (This does not include the "leaf" in green by default).
    writer = vtk.vtkOBJWriter()
    writer.SetFileName("extractSurface.obj")
    writer.AddInputConnection(tubeFilter.GetOutputPort())
    writer.Write()
    
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
    
    lut = create_lookup_table()
    if ReturnLut and RangeExtra != None:
        lut = create_lookup_table(tableIdx = "BREWER_SEQUENTIAL_YELLOW_ORANGE_BROWN_9") # it is from https://vtk.org/doc/nightly/html/vtkColorSeries_8h_source.html
        a = pd.GetCellData().GetAbstractArray(p_name)
        RangeTemp = a.GetRange()
        RangeLow = min(RangeExtra[0], RangeTemp[0])
        RangeHigh =  max(RangeExtra[1], RangeTemp[1])
        RangeOverwrite = (RangeLow, RangeHigh)
        scalar_bar = create_scalar_bar(lut, pd, p_name, RangeOverwrite)  # vtkScalarBarActor
    else:
        lut = create_lookup_table()
        scalar_bar = create_scalar_bar(lut, pd, p_name)
    mapper.SetLookupTable(lut)

    if render:
        render_window(plantActor, win_title, scalar_bar, pd.GetBounds()).Start()
    if ReturnLut:
        return plantActor, scalar_bar, lut
    else:
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
    # mapper.SetArrayName(p_name)
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
    return [meshActor], scalar_bar


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

    if render:
        render_window(actors, win_title, scalar_bar, grid.GetBounds()).Start()

    return actors, scalar_bar


def plot_mesh_yz(s, p_name , min_b, max_b, cell_number, x = -4, render = False):
    """ TODO """
    grid = uniform_grid(np.array(min_b), np.array(max_b), np.array(cell_number))
    soil_water_content = vtk_data(np.array(s.getWaterContent()))
    soil_water_content.SetName("water content")
    grid.GetCellData().AddArray(soil_water_content)
    soil_pressure = vtk_data(np.array(s.getSolutionHead()))
    soil_pressure.SetName("pressure head")  # in macroscopic soil
    grid.GetCellData().AddArray(soil_pressure)
    bounds = grid.GetBounds()
    p = vtk.vtkPlane()
    p.SetOrigin(x, 0, 0)
    p.SetNormal(1, 0, 0)

    lut = create_lookup_table()
    scalar_bar = create_scalar_bar(lut, grid, p_name)

    actors = []  # create cutter, mappers, and actors
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

    return actors, scalar_bar, grid


def plot_roots_and_soil(rs, pname:str, rp, s, periodic:bool, min_b, max_b, cell_number, filename:str, sol_ind = 0):
    """ Plots soil slices and roots, additionally saves both grids as files
    @param rs            some Organism (e.g. RootSystem, MappedRootSystem, ...) or MappedSegments
    @param pname         root and soil parameter that will be visualized ("pressure head", or "water content")
    @param s
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
        w = np.array(max_b) - np.array(min_b)
        ana.mapPeriodic(w[0], w[1])
    pd = segs_to_polydata(ana, 1., ["radius", "subType", "creationTime", pname])

    pname_mesh = pname
    soil_grid = uniform_grid(np.array(min_b), np.array(max_b), np.array(cell_number))
    soil_water_content = vtk_data(np.array(s.getWaterContent()))
    soil_water_content.SetName("water content")
    soil_grid.GetCellData().AddArray(soil_water_content)
    soil_pressure = vtk_data(np.array(s.getSolutionHead()))
    soil_pressure.SetName("pressure head")  # in macroscopic soil
    soil_grid.GetCellData().AddArray(soil_pressure)
    if sol_ind > 0:
        d = vtk_data(np.array(s.getSolution(sol_ind)))
        pname_mesh = "solute" + str(sol_ind)
        d.SetName(pname_mesh)  # in macroscopic soil
        soil_grid.GetCellData().AddArray(d)

    rootActor, rootCBar = plot_roots(pd, pname, "", False)
    meshActors, meshCBar = plot_mesh_cuts(soil_grid, pname_mesh, 7, "", False)
    lut = meshActors[-1].GetMapper().GetLookupTable()  # same same
    rootActor.GetMapper().SetLookupTable(lut)
    meshActors.extend([rootActor])
    render_window(meshActors, filename, meshCBar, pd.GetBounds()).Start()

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
    write_vtu(filename + ".vtu", soil_grid)


def plot_roots_and_soil_files(filename: str, pname:str):
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
    render_window(meshActors, filename, meshCBar, soil_grid.GetBounds()).Start()


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

