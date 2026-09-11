from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()
renderView1 = GetActiveViewOrCreate('RenderView')

obj1 = Cylinder()

obj1.Resolution = 50
obj1.Height = 20
obj1.Radius = 2.5
obj1.Center = [0., -10,0.]

obj2 = Transform(Input=obj1)
obj2.Transform = 'Transform'
obj2.Transform.Rotate = [90.0, 0.0, 0.0]

obj1Display = Show(obj2,renderView1)
obj1Display.Opacity = 0.2
obj1Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj2)
obj3= Transform(Input=obj2)
obj3.Transform = 'Transform'
obj3.Transform.Translate = [10,0,-20]
obj3.Transform.Rotate = [0,0,0]

obj3Display = Show(obj3,renderView1)
obj3Display.Opacity = 0.1
obj3Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
