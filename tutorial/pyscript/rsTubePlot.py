#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find source
dic = GetSources();
keys = dic.keys();
sorghumvtp = FindSource(keys[0][0])

# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(Input=sorghumvtp)

# Properties modified on cellDatatoPointData1
cellDatatoPointData1.PassCellData = 1

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
cellDatatoPointData1Display = Show(cellDatatoPointData1, renderView1)
# trace defaults for the display properties.
cellDatatoPointData1Display.ColorArrayName = [None, '']

# hide data in view
Hide(sorghumvtp, renderView1)

# create a new 'Tube'
tube1 = Tube(Input=cellDatatoPointData1)
tube1.Scalars = ['POINTS', 'radius']
tube1.Vectors = [None, '1']
tube1.Radius = 0.32146133262802323

# Properties modified on tube1
tube1.Vectors = [None, '']
tube1.VaryRadius = 'By Absolute Scalar'

# show data in view
tube1Display = Show(tube1, renderView1)
# trace defaults for the display properties.
tube1Display.ColorArrayName = [None, '']

# hide data in view
Hide(cellDatatoPointData1, renderView1)

# set scalar coloring
ColorBy(tube1Display, ('CELLS', 'subType'))

# rescale color and/or opacity maps used to include current data range
tube1Display.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
tube1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'type'
typeLUT = GetColorTransferFunction('subType')

# get opacity transfer function/opacity map for 'type'
typePWF = GetOpacityTransferFunction('subType')

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
typeLUT.ApplyPreset('blot', True)
