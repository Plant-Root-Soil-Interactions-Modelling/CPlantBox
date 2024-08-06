from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()
renderView1 = GetActiveViewOrCreate('RenderView')

obj1 = Cylinder()

obj1.Resolution = 50
obj1.Height = 40
obj1.Radius = 5
obj1.Center = [0., -20,0.]

obj2 = Transform(Input=obj1)
obj2.Transform = 'Transform'
obj2.Transform.Rotate = [90.0, 0.0, 0.0]

obj1Display = Show(obj2,renderView1)
obj1Display.Opacity = 0.2
obj1Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
