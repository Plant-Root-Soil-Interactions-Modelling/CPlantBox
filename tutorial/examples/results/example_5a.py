from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()
renderView1 = GetActiveViewOrCreate('RenderView')

# GROUP 
obj1 = Box()

obj1.XLength = 10
obj1.YLength = 20
obj1.ZLength = 50
obj1.Center = [0.,0., -25]

obj1Display = Show(obj1,renderView1)
obj1Display.Opacity = 0.1
obj1Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj1)
obj2= Transform(Input=obj1)
obj2.Transform = 'Transform'
obj2.Transform.Translate = [-4.99,0,0]
obj2.Transform.Rotate = [0,0,0]

obj2Display = Show(obj2,renderView1)
obj2Display.Opacity = 0.1
obj2Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
obj3 = Box()

obj3.XLength = 10
obj3.YLength = 20
obj3.ZLength = 50
obj3.Center = [0.,0., -25]

obj3Display = Show(obj3,renderView1)
obj3Display.Opacity = 0.1
obj3Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj3)
obj4= Transform(Input=obj3)
obj4.Transform = 'Transform'
obj4.Transform.Translate = [4.99,0,0]
obj4.Transform.Rotate = [0,0,0]

obj4Display = Show(obj4,renderView1)
obj4Display.Opacity = 0.1
obj4Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#

obj5= GroupDatasets(Input=[obj2,obj4]) 
obj5Display = Show(obj5, renderView1)
obj5Display.Opacity = 0.1
obj5Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
