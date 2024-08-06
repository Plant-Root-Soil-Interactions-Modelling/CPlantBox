from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()
renderView1 = GetActiveViewOrCreate('RenderView')

# GROUP 
obj1 = Box()

obj1.XLength = 22
obj1.YLength = 20
obj1.ZLength = 5
obj1.Center = [0.,0., -2.5]

obj1Display = Show(obj1,renderView1)
obj1Display.Opacity = 0.1
obj1Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
obj2 = Box()

obj2.XLength = 10
obj2.YLength = 20
obj2.ZLength = 35
obj2.Center = [0.,0., -17.5]

obj2Display = Show(obj2,renderView1)
obj2Display.Opacity = 0.1
obj2Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj2)
obj3= Transform(Input=obj2)
obj3.Transform = 'Transform'
obj3.Transform.Translate = [-6,0,-5]
obj3.Transform.Rotate = [0,0,0]

obj3Display = Show(obj3,renderView1)
obj3Display.Opacity = 0.1
obj3Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
obj4 = Box()

obj4.XLength = 10
obj4.YLength = 20
obj4.ZLength = 35
obj4.Center = [0.,0., -17.5]

obj4Display = Show(obj4,renderView1)
obj4Display.Opacity = 0.1
obj4Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj4)
obj5= Transform(Input=obj4)
obj5.Transform = 'Transform'
obj5.Transform.Translate = [6,0,-5]
obj5.Transform.Rotate = [0,0,0]

obj5Display = Show(obj5,renderView1)
obj5Display.Opacity = 0.1
obj5Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#

obj6= GroupDatasets(Input=[obj1,obj3,obj5]) 
obj6Display = Show(obj6, renderView1)
obj6Display.Opacity = 0.1
obj6Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
