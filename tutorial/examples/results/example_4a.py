from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()
renderView1 = GetActiveViewOrCreate('RenderView')

obj1 = Box()

obj1.XLength = 30
obj1.YLength = 30
obj1.ZLength = 2
obj1.Center = [0.,0., -1]

obj1Display = Show(obj1,renderView1)
obj1Display.Opacity = 0.1
obj1Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj1)
obj2= Transform(Input=obj1)
obj2.Transform = 'Transform'
obj2.Transform.Translate = [0,0,-16]
obj2.Transform.Rotate = [0,0,0]

obj2Display = Show(obj2,renderView1)
obj2Display.Opacity = 0.1
obj2Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
