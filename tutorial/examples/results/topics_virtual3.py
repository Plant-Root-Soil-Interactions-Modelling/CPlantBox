from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()
renderView1 = GetActiveViewOrCreate('RenderView')

# GROUP 
obj1 = Box()

obj1.XLength = 96
obj1.YLength = 126
obj1.ZLength = 130
obj1.Center = [0.,0., -65]

obj1Display = Show(obj1,renderView1)
obj1Display.Opacity = 0.1
obj1Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
# GROUP 
obj2 = Cylinder()

obj2.Resolution = 50
obj2.Height = 96
obj2.Radius = 3
obj2.Center = [0., -48,0.]

obj3 = Transform(Input=obj2)
obj3.Transform = 'Transform'
obj3.Transform.Rotate = [90.0, 0.0, 0.0]

obj2Display = Show(obj3,renderView1)
obj2Display.Opacity = 0.2
obj2Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj3)
obj4= Transform(Input=obj3)
obj4.Transform = 'Transform'
obj4.Transform.Translate = [48,0,0]
obj4.Transform.Rotate = [0,90,0]

obj4Display = Show(obj4,renderView1)
obj4Display.Opacity = 0.1
obj4Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj4)
obj5= Transform(Input=obj4)
obj5.Transform = 'Transform'
obj5.Transform.Translate = [0,-30,-10]
obj5.Transform.Rotate = [0,0,0]

obj5Display = Show(obj5,renderView1)
obj5Display.Opacity = 0.1
obj5Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
obj6 = Cylinder()

obj6.Resolution = 50
obj6.Height = 96
obj6.Radius = 3
obj6.Center = [0., -48,0.]

obj7 = Transform(Input=obj6)
obj7.Transform = 'Transform'
obj7.Transform.Rotate = [90.0, 0.0, 0.0]

obj6Display = Show(obj7,renderView1)
obj6Display.Opacity = 0.2
obj6Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj7)
obj8= Transform(Input=obj7)
obj8.Transform = 'Transform'
obj8.Transform.Translate = [48,0,0]
obj8.Transform.Rotate = [0,90,0]

obj8Display = Show(obj8,renderView1)
obj8Display.Opacity = 0.1
obj8Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj8)
obj9= Transform(Input=obj8)
obj9.Transform = 'Transform'
obj9.Transform.Translate = [0,-18,-20]
obj9.Transform.Rotate = [0,0,0]

obj9Display = Show(obj9,renderView1)
obj9Display.Opacity = 0.1
obj9Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
obj10 = Cylinder()

obj10.Resolution = 50
obj10.Height = 96
obj10.Radius = 3
obj10.Center = [0., -48,0.]

obj11 = Transform(Input=obj10)
obj11.Transform = 'Transform'
obj11.Transform.Rotate = [90.0, 0.0, 0.0]

obj10Display = Show(obj11,renderView1)
obj10Display.Opacity = 0.2
obj10Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj11)
obj12= Transform(Input=obj11)
obj12.Transform = 'Transform'
obj12.Transform.Translate = [48,0,0]
obj12.Transform.Rotate = [0,90,0]

obj12Display = Show(obj12,renderView1)
obj12Display.Opacity = 0.1
obj12Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj12)
obj13= Transform(Input=obj12)
obj13.Transform = 'Transform'
obj13.Transform.Translate = [0,-6,-40]
obj13.Transform.Rotate = [0,0,0]

obj13Display = Show(obj13,renderView1)
obj13Display.Opacity = 0.1
obj13Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
obj14 = Cylinder()

obj14.Resolution = 50
obj14.Height = 96
obj14.Radius = 3
obj14.Center = [0., -48,0.]

obj15 = Transform(Input=obj14)
obj15.Transform = 'Transform'
obj15.Transform.Rotate = [90.0, 0.0, 0.0]

obj14Display = Show(obj15,renderView1)
obj14Display.Opacity = 0.2
obj14Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj15)
obj16= Transform(Input=obj15)
obj16.Transform = 'Transform'
obj16.Transform.Translate = [48,0,0]
obj16.Transform.Rotate = [0,90,0]

obj16Display = Show(obj16,renderView1)
obj16Display.Opacity = 0.1
obj16Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj16)
obj17= Transform(Input=obj16)
obj17.Transform = 'Transform'
obj17.Transform.Translate = [0,6,-60]
obj17.Transform.Rotate = [0,0,0]

obj17Display = Show(obj17,renderView1)
obj17Display.Opacity = 0.1
obj17Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
obj18 = Cylinder()

obj18.Resolution = 50
obj18.Height = 96
obj18.Radius = 3
obj18.Center = [0., -48,0.]

obj19 = Transform(Input=obj18)
obj19.Transform = 'Transform'
obj19.Transform.Rotate = [90.0, 0.0, 0.0]

obj18Display = Show(obj19,renderView1)
obj18Display.Opacity = 0.2
obj18Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj19)
obj20= Transform(Input=obj19)
obj20.Transform = 'Transform'
obj20.Transform.Translate = [48,0,0]
obj20.Transform.Rotate = [0,90,0]

obj20Display = Show(obj20,renderView1)
obj20Display.Opacity = 0.1
obj20Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj20)
obj21= Transform(Input=obj20)
obj21.Transform = 'Transform'
obj21.Transform.Translate = [0,18,-80]
obj21.Transform.Rotate = [0,0,0]

obj21Display = Show(obj21,renderView1)
obj21Display.Opacity = 0.1
obj21Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
obj22 = Cylinder()

obj22.Resolution = 50
obj22.Height = 96
obj22.Radius = 3
obj22.Center = [0., -48,0.]

obj23 = Transform(Input=obj22)
obj23.Transform = 'Transform'
obj23.Transform.Rotate = [90.0, 0.0, 0.0]

obj22Display = Show(obj23,renderView1)
obj22Display.Opacity = 0.2
obj22Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj23)
obj24= Transform(Input=obj23)
obj24.Transform = 'Transform'
obj24.Transform.Translate = [48,0,0]
obj24.Transform.Rotate = [0,90,0]

obj24Display = Show(obj24,renderView1)
obj24Display.Opacity = 0.1
obj24Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj24)
obj25= Transform(Input=obj24)
obj25.Transform = 'Transform'
obj25.Transform.Translate = [0,30,-120]
obj25.Transform.Rotate = [0,0,0]

obj25Display = Show(obj25,renderView1)
obj25Display.Opacity = 0.1
obj25Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#

obj26= GroupDatasets(Input=[obj5,obj9,obj13,obj17,obj21,obj25]) 
obj26Display = Show(obj26, renderView1)
obj26Display.Opacity = 0.1
obj26Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#

obj27= GroupDatasets(Input=[obj1,obj26]) 
obj27Display = Show(obj27, renderView1)
obj27Display.Opacity = 0.1
obj27Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
