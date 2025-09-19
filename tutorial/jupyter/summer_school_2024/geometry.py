from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()
renderView1 = GetActiveViewOrCreate('RenderView')

# GROUP 
obj1 = Box()

obj1.XLength = 20
obj1.YLength = 10
obj1.ZLength = 120
obj1.Center = [0.,0., -60]

obj1Display = Show(obj1,renderView1)
obj1Display.Opacity = 0.1
obj1Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
obj2 = Cylinder()

obj2.Resolution = 50
obj2.Height = 120
obj2.Radius = 2
obj2.Center = [0., -60,0.]

obj3 = Transform(Input=obj2)
obj3.Transform = 'Transform'
obj3.Transform.Rotate = [90.0, 0.0, 0.0]

obj2Display = Show(obj3,renderView1)
obj2Display.Opacity = 0.2
obj2Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
obj4 = Cylinder()

obj4.Resolution = 50
obj4.Height = 120
obj4.Radius = 2
obj4.Center = [0., -60,0.]

obj5 = Transform(Input=obj4)
obj5.Transform = 'Transform'
obj5.Transform.Rotate = [90.0, 0.0, 0.0]

obj4Display = Show(obj5,renderView1)
obj4Display.Opacity = 0.2
obj4Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj5)
obj6= Transform(Input=obj5)
obj6.Transform = 'Transform'
obj6.Transform.Translate = [4,0,0]
obj6.Transform.Rotate = [0,0,0]

obj6Display = Show(obj6,renderView1)
obj6Display.Opacity = 0.1
obj6Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#

obj7= GroupDatasets(Input=[obj1,obj3,obj6]) 
obj7Display = Show(obj7, renderView1)
obj7Display.Opacity = 0.1
obj7Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
