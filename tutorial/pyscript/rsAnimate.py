# trace generated using paraview version 5.7.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find source
tube1 = FindSource('Tube1')

# create a new 'Threshold'
threshold1 = Threshold(Input=tube1)
threshold1.Scalars = ['POINTS', 'creationTime']
threshold1.ThresholdRange = [0.0, 19.76020050048828]

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [2172, 906]

# hide data in view
Hide(tube1, renderView1)

# set active source
SetActiveSource(threshold1)

# show data in view
threshold1Display = Show(threshold1, renderView1)

# get color transfer function/color map for 'subType'
subTypeLUT = GetColorTransferFunction('subType')

# get opacity transfer function/opacity map for 'subType'
subTypePWF = GetOpacityTransferFunction('subType')

# trace defaults for the display properties.
threshold1Display.Representation = 'Surface'
threshold1Display.ColorArrayName = ['POINTS', 'subType']
threshold1Display.LookupTable = subTypeLUT
threshold1Display.OSPRayScaleArray = 'TubeNormals'
threshold1Display.OSPRayScaleFunction = 'PiecewiseFunction'
threshold1Display.SelectOrientationVectors = 'None'
threshold1Display.ScaleFactor = 2.610576033592224
threshold1Display.SelectScaleArray = 'None'
threshold1Display.GlyphType = 'Arrow'
threshold1Display.GlyphTableIndexArray = 'None'
threshold1Display.GaussianRadius = 0.1305288016796112
threshold1Display.SetScaleArray = ['POINTS', 'TubeNormals']
threshold1Display.ScaleTransferFunction = 'PiecewiseFunction'
threshold1Display.OpacityArray = ['POINTS', 'TubeNormals']
threshold1Display.OpacityTransferFunction = 'PiecewiseFunction'
threshold1Display.DataAxesGrid = 'GridAxesRepresentation'
threshold1Display.PolarAxes = 'PolarAxesRepresentation'
threshold1Display.ScalarOpacityFunction = subTypePWF
threshold1Display.ScalarOpacityUnitDistance = 1.8903992383989479

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
threshold1Display.ScaleTransferFunction.Points = [-0.9999293684959412, 0.0, 0.5, 0.0, 0.9999293684959412, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
threshold1Display.OpacityTransferFunction.Points = [-0.9999293684959412, 0.0, 0.5, 0.0, 0.9999293684959412, 1.0, 0.5, 0.0]

# show color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, True)

# reset view to fit data
renderView1.ResetCamera()

# get animation track
threshold1ThresholdBetweenTrack = GetAnimationTrack('ThresholdBetween', index=1, proxy=threshold1)

# create keyframes for this animation track

# create a key frame
keyFrame6441 = CompositeKeyFrame()

# create a key frame
keyFrame6442 = CompositeKeyFrame()
keyFrame6442.KeyTime = 1.0
keyFrame6442.KeyValues = [19.76020050048828]

# initialize the animation track
threshold1ThresholdBetweenTrack.KeyFrames = [keyFrame6441, keyFrame6442]

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# Properties modified on animationScene1
animationScene1.EndTime = 30.0

# Properties modified on animationScene1
animationScene1.NumberOfFrames = 300

animationScene1.Play()

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [-73.91243369541311, -35.063920045750486, 4.5490810279582305]
renderView1.CameraFocalPoint = [-2.781435012817383, 0.2826528549194336, -16.05103600025177]
renderView1.CameraViewUp = [0.2135400977846104, 0.1342352329335891, 0.9676680881776583]
renderView1.CameraParallelScale = 21.2379345104894

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).