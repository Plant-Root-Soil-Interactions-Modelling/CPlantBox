# trace generated using paraview version 5.7.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find source
dic = GetSources()
keys = list(dic.keys())
vtp_file = FindSource(keys[0][0])

# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(Input=vtp_file)
cellDatatoPointData1.CellDataArraytoprocess = ['creationTime', 'radius', 'subType']

# Properties modified on cellDatatoPointData1
cellDatatoPointData1.PassCellData = 1

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1970, 907]

# show data in view
cellDatatoPointData1Display = Show(cellDatatoPointData1, renderView1)

# trace defaults for the display properties.
cellDatatoPointData1Display.Representation = 'Surface'
cellDatatoPointData1Display.ColorArrayName = [None, '']
cellDatatoPointData1Display.OSPRayScaleArray = 'creationTime'
cellDatatoPointData1Display.OSPRayScaleFunction = 'PiecewiseFunction'
cellDatatoPointData1Display.SelectOrientationVectors = 'None'
cellDatatoPointData1Display.ScaleFactor = 2.6079999923706056
cellDatatoPointData1Display.SelectScaleArray = 'None'
cellDatatoPointData1Display.GlyphType = 'Arrow'
cellDatatoPointData1Display.GlyphTableIndexArray = 'None'
cellDatatoPointData1Display.GaussianRadius = 0.13039999961853027
cellDatatoPointData1Display.SetScaleArray = ['POINTS', 'creationTime']
cellDatatoPointData1Display.ScaleTransferFunction = 'PiecewiseFunction'
cellDatatoPointData1Display.OpacityArray = ['POINTS', 'creationTime']
cellDatatoPointData1Display.OpacityTransferFunction = 'PiecewiseFunction'
cellDatatoPointData1Display.DataAxesGrid = 'GridAxesRepresentation'
cellDatatoPointData1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
cellDatatoPointData1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 19.76020050048828, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
cellDatatoPointData1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 19.76020050048828, 1.0, 0.5, 0.0]

# hide data in view
Hide(vtp_file, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Tube'
tube1 = Tube(Input=cellDatatoPointData1)
tube1.Scalars = ['POINTS', 'creationTime']
tube1.Vectors = [None, '1']
tube1.Radius = 0.26079999923706054

# set active source
SetActiveSource(tube1)

# show data in view
tube1Display = Show(tube1, renderView1)

# trace defaults for the display properties.
tube1Display.Representation = 'Surface'
tube1Display.ColorArrayName = [None, '']
tube1Display.OSPRayScaleArray = 'TubeNormals'
tube1Display.OSPRayScaleFunction = 'PiecewiseFunction'
tube1Display.SelectOrientationVectors = 'None'
tube1Display.ScaleFactor = 2.621436905860901
tube1Display.SelectScaleArray = 'None'
tube1Display.GlyphType = 'Arrow'
tube1Display.GlyphTableIndexArray = 'None'
tube1Display.GaussianRadius = 0.13107184529304505
tube1Display.SetScaleArray = ['POINTS', 'TubeNormals']
tube1Display.ScaleTransferFunction = 'PiecewiseFunction'
tube1Display.OpacityArray = ['POINTS', 'TubeNormals']
tube1Display.OpacityTransferFunction = 'PiecewiseFunction'
tube1Display.DataAxesGrid = 'GridAxesRepresentation'
tube1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
tube1Display.ScaleTransferFunction.Points = [-0.9999293684959412, 0.0, 0.5, 0.0, 0.9999293684959412, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
tube1Display.OpacityTransferFunction.Points = [-0.9999293684959412, 0.0, 0.5, 0.0, 0.9999293684959412, 1.0, 0.5, 0.0]

# hide data in view
Hide(cellDatatoPointData1, renderView1)

# Properties modified on tube1
tube1.Scalars = ['POINTS', 'radius']
tube1.Vectors = ['POINTS', '1']
tube1.VaryRadius = 'By Absolute Scalar'

# show data in view
tube1Display = Show(tube1, renderView1)

# reset view to fit data
renderView1.ResetCamera()

# hide data in view
Hide(cellDatatoPointData1, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(tube1Display, ('POINTS', 'subType'))

# rescale color and/or opacity maps used to include current data range
tube1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
tube1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'subType'
subTypeLUT = GetColorTransferFunction('subType')

# get opacity transfer function/opacity map for 'subType'
subTypePWF = GetOpacityTransferFunction('subType')

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [-51.36486420658536, -23.859532037280484, -1.980878887889233]
renderView1.CameraFocalPoint = [-2.781435012817383, 0.2826528549194336, -16.05103600025177]
renderView1.CameraViewUp = [0.2135400977846104, 0.1342352329335891, 0.9676680881776583]
renderView1.CameraParallelScale = 21.2379345104894

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
