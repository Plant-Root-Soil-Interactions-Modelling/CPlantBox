#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find source
tube1 = FindSource('Tube1')

# set active source
SetActiveSource(tube1)

# create a new 'Threshold'
threshold1 = Threshold(Input=tube1)
threshold1.Scalars = ['POINTS', 'order']
threshold1.ThresholdRange = [0.0, 3.0]

# Properties modified on threshold1
threshold1.Scalars = ['CELLS', 'creationTime']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [930, 504]
typeLUT = GetColorTransferFunction('organtype')
# show data in view
threshold1Display = Show(threshold1, renderView1)
# trace defaults for the display properties.
threshold1Display.ColorArrayName = ['CELLS', 'organType']
threshold1Display.LookupTable = typeLUT
threshold1Display.ScalarOpacityUnitDistance = 0.4049429502949595

# hide data in view
tube1 = FindSource('Tube1')
Hide(tube1, renderView1)

# show color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, True)

# Properties modified on threshold1
threshold1.ThresholdRange = [0.0, 49.24366475880146]

# get animation track
threshold1ThresholdBetweenTrack = GetAnimationTrack('ThresholdBetween', index=1, proxy=threshold1)

# create keyframes for this animation track

# create a key frame
keyFrame8395 = CompositeKeyFrame()
keyFrame8395.KeyValues = [0.24258199334144592]

# create a key frame
keyFrame8396 = CompositeKeyFrame()
keyFrame8396.KeyTime = 1.0
keyFrame8396.KeyValues = [60.0]

# initialize the animation track
threshold1ThresholdBetweenTrack.KeyFrames = [keyFrame8395, keyFrame8396]

# get animation scene
animationScene1 = GetAnimationScene()

# Properties modified on animationScene1
animationScene1.NumberOfFrames = 300

# Properties modified on animationScene1
animationScene1.AnimationTime = 0.0
