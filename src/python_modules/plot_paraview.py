
from paraview.simple import *  # use paraview python interpreter 


def plot_roots(rs):

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    # uncomment following to set a specific view size
    # renderView1.ViewSize = [1976, 908]
        
    layout1 = GetLayout()  # get layout
    
    polylines = rs.getPolylines()

    print(len(polylines))

    polyLineSources = []

    for pl in polylines: 
    
    # for pl in [polylines[0], polylines[1], polylines[2]]:
    # for i in range(0, 100):
        # pl = polylines[i]
        
        pl_ = []
        for p in pl:
            pl_.extend([p.x, p.y, p.z])  # convert from Vector3 to list        
        
        if len(pl_) > 2:
            s = PolyLineSource()
            polyLineSources.append(s)        
            s.Points = list(pl_)
            
#             polyLineSource1Display = Show(polyLineSource, renderView1)  # , 'GeometryRepresentation'
#             polyLineSource1Display.Representation = 'Surface'
#             polyLineSource1Display.ColorArrayName = [None, '']
#             polyLineSource1Display.OSPRayScaleFunction = 'PiecewiseFunction'
#             polyLineSource1Display.SelectOrientationVectors = 'None'
#             polyLineSource1Display.ScaleFactor = 0.4
#             polyLineSource1Display.SelectScaleArray = 'None'
#             polyLineSource1Display.GlyphType = 'Arrow'
#             polyLineSource1Display.GlyphTableIndexArray = 'None'
#             polyLineSource1Display.GaussianRadius = 0.02
#             polyLineSource1Display.SetScaleArray = [None, '']
#             polyLineSource1Display.ScaleTransferFunction = 'PiecewiseFunction'
#             polyLineSource1Display.OpacityArray = [None, '']
#             polyLineSource1Display.OpacityTransferFunction = 'PiecewiseFunction'
#             polyLineSource1Display.DataAxesGrid = 'GridAxesRepresentation'
#             polyLineSource1Display.PolarAxes = 'PolarAxesRepresentation'        

    # reset view to fit data
    renderView1.ResetCamera()
    
    # get the material library
    materialLibrary1 = GetMaterialLibrary()
    
    # update the view to ensure updated data information
    renderView1.Update()
