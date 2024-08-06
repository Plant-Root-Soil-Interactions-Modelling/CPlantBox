from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()
renderView1 = GetActiveViewOrCreate('RenderView')

obj1 = Box()

obj1.XLength = 1.4
obj1.YLength = 27
obj1.ZLength = 27
obj1.Center = [0.,0., -13.5]

obj1Display = Show(obj1,renderView1)
obj1Display.Opacity = 0.1
obj1Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
