from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()
renderView1 = GetActiveViewOrCreate('RenderView')

# GROUP 
obj1 = Cylinder()

obj1.Resolution = 50
obj1.Height = 160
obj1.Radius = 2.1
obj1.Center = [0., -80,0.]

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
obj3.Transform.Translate = [0,0,0]
obj3.Transform.Rotate = [0,0,0]

obj3Display = Show(obj3,renderView1)
obj3Display.Opacity = 0.2
obj3Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
obj4 = Cylinder()

obj4.Resolution = 50
obj4.Height = 160
obj4.Radius = 2.1
obj4.Center = [0., -80,0.]

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
obj6.Transform.Translate = [1,1,0]
obj6.Transform.Rotate = [0,0,0]

obj6Display = Show(obj6,renderView1)
obj6Display.Opacity = 0.2
obj6Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
obj7 = Cylinder()

obj7.Resolution = 50
obj7.Height = 160
obj7.Radius = 2.1
obj7.Center = [0., -80,0.]

obj8 = Transform(Input=obj7)
obj8.Transform = 'Transform'
obj8.Transform.Rotate = [90.0, 0.0, 0.0]

obj7Display = Show(obj8,renderView1)
obj7Display.Opacity = 0.2
obj7Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj8)
obj9= Transform(Input=obj8)
obj9.Transform = 'Transform'
obj9.Transform.Translate = [2,2,0]
obj9.Transform.Rotate = [0,0,0]

obj9Display = Show(obj9,renderView1)
obj9Display.Opacity = 0.2
obj9Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#

obj10= GroupDatasets(Input=[obj3,obj6,obj9]) 
obj10Display = Show(obj10, renderView1)
obj10Display.Opacity = 0.2
obj10Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
