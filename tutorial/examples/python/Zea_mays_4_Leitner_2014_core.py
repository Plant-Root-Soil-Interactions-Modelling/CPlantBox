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
obj3.Transform.Translate = [27,30,0]
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
obj6.Transform.Translate = [27,42,0]
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
obj9.Transform.Translate = [27,54,0]
obj9.Transform.Rotate = [0,0,0]

obj9Display = Show(obj9,renderView1)
obj9Display.Opacity = 0.2
obj9Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
obj10 = Cylinder()

obj10.Resolution = 50
obj10.Height = 160
obj10.Radius = 2.1
obj10.Center = [0., -80,0.]

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
obj12.Transform.Translate = [27,66,0]
obj12.Transform.Rotate = [0,0,0]

obj12Display = Show(obj12,renderView1)
obj12Display.Opacity = 0.2
obj12Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
obj13 = Cylinder()

obj13.Resolution = 50
obj13.Height = 160
obj13.Radius = 2.1
obj13.Center = [0., -80,0.]

obj14 = Transform(Input=obj13)
obj14.Transform = 'Transform'
obj14.Transform.Rotate = [90.0, 0.0, 0.0]

obj13Display = Show(obj14,renderView1)
obj13Display.Opacity = 0.2
obj13Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj14)
obj15= Transform(Input=obj14)
obj15.Transform = 'Transform'
obj15.Transform.Translate = [27,78,0]
obj15.Transform.Rotate = [0,0,0]

obj15Display = Show(obj15,renderView1)
obj15Display.Opacity = 0.2
obj15Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
obj16 = Cylinder()

obj16.Resolution = 50
obj16.Height = 160
obj16.Radius = 2.1
obj16.Center = [0., -80,0.]

obj17 = Transform(Input=obj16)
obj17.Transform = 'Transform'
obj17.Transform.Rotate = [90.0, 0.0, 0.0]

obj16Display = Show(obj17,renderView1)
obj16Display.Opacity = 0.2
obj16Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj17)
obj18= Transform(Input=obj17)
obj18.Transform = 'Transform'
obj18.Transform.Translate = [45,30,0]
obj18.Transform.Rotate = [0,0,0]

obj18Display = Show(obj18,renderView1)
obj18Display.Opacity = 0.2
obj18Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
obj19 = Cylinder()

obj19.Resolution = 50
obj19.Height = 160
obj19.Radius = 2.1
obj19.Center = [0., -80,0.]

obj20 = Transform(Input=obj19)
obj20.Transform = 'Transform'
obj20.Transform.Rotate = [90.0, 0.0, 0.0]

obj19Display = Show(obj20,renderView1)
obj19Display.Opacity = 0.2
obj19Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj20)
obj21= Transform(Input=obj20)
obj21.Transform = 'Transform'
obj21.Transform.Translate = [45,42,0]
obj21.Transform.Rotate = [0,0,0]

obj21Display = Show(obj21,renderView1)
obj21Display.Opacity = 0.2
obj21Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
obj22 = Cylinder()

obj22.Resolution = 50
obj22.Height = 160
obj22.Radius = 2.1
obj22.Center = [0., -80,0.]

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
obj24.Transform.Translate = [45,54,0]
obj24.Transform.Rotate = [0,0,0]

obj24Display = Show(obj24,renderView1)
obj24Display.Opacity = 0.2
obj24Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
obj25 = Cylinder()

obj25.Resolution = 50
obj25.Height = 160
obj25.Radius = 2.1
obj25.Center = [0., -80,0.]

obj26 = Transform(Input=obj25)
obj26.Transform = 'Transform'
obj26.Transform.Rotate = [90.0, 0.0, 0.0]

obj25Display = Show(obj26,renderView1)
obj25Display.Opacity = 0.2
obj25Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj26)
obj27= Transform(Input=obj26)
obj27.Transform = 'Transform'
obj27.Transform.Translate = [45,66,0]
obj27.Transform.Rotate = [0,0,0]

obj27Display = Show(obj27,renderView1)
obj27Display.Opacity = 0.2
obj27Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
obj28 = Cylinder()

obj28.Resolution = 50
obj28.Height = 160
obj28.Radius = 2.1
obj28.Center = [0., -80,0.]

obj29 = Transform(Input=obj28)
obj29.Transform = 'Transform'
obj29.Transform.Rotate = [90.0, 0.0, 0.0]

obj28Display = Show(obj29,renderView1)
obj28Display.Opacity = 0.2
obj28Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj29)
obj30= Transform(Input=obj29)
obj30.Transform = 'Transform'
obj30.Transform.Translate = [45,78,0]
obj30.Transform.Rotate = [0,0,0]

obj30Display = Show(obj30,renderView1)
obj30Display.Opacity = 0.2
obj30Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
obj31 = Cylinder()

obj31.Resolution = 50
obj31.Height = 160
obj31.Radius = 2.1
obj31.Center = [0., -80,0.]

obj32 = Transform(Input=obj31)
obj32.Transform = 'Transform'
obj32.Transform.Rotate = [90.0, 0.0, 0.0]

obj31Display = Show(obj32,renderView1)
obj31Display.Opacity = 0.2
obj31Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj32)
obj33= Transform(Input=obj32)
obj33.Transform = 'Transform'
obj33.Transform.Translate = [63,30,0]
obj33.Transform.Rotate = [0,0,0]

obj33Display = Show(obj33,renderView1)
obj33Display.Opacity = 0.2
obj33Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
obj34 = Cylinder()

obj34.Resolution = 50
obj34.Height = 160
obj34.Radius = 2.1
obj34.Center = [0., -80,0.]

obj35 = Transform(Input=obj34)
obj35.Transform = 'Transform'
obj35.Transform.Rotate = [90.0, 0.0, 0.0]

obj34Display = Show(obj35,renderView1)
obj34Display.Opacity = 0.2
obj34Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj35)
obj36= Transform(Input=obj35)
obj36.Transform = 'Transform'
obj36.Transform.Translate = [63,42,0]
obj36.Transform.Rotate = [0,0,0]

obj36Display = Show(obj36,renderView1)
obj36Display.Opacity = 0.2
obj36Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
obj37 = Cylinder()

obj37.Resolution = 50
obj37.Height = 160
obj37.Radius = 2.1
obj37.Center = [0., -80,0.]

obj38 = Transform(Input=obj37)
obj38.Transform = 'Transform'
obj38.Transform.Rotate = [90.0, 0.0, 0.0]

obj37Display = Show(obj38,renderView1)
obj37Display.Opacity = 0.2
obj37Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj38)
obj39= Transform(Input=obj38)
obj39.Transform = 'Transform'
obj39.Transform.Translate = [63,54,0]
obj39.Transform.Rotate = [0,0,0]

obj39Display = Show(obj39,renderView1)
obj39Display.Opacity = 0.2
obj39Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
obj40 = Cylinder()

obj40.Resolution = 50
obj40.Height = 160
obj40.Radius = 2.1
obj40.Center = [0., -80,0.]

obj41 = Transform(Input=obj40)
obj41.Transform = 'Transform'
obj41.Transform.Rotate = [90.0, 0.0, 0.0]

obj40Display = Show(obj41,renderView1)
obj40Display.Opacity = 0.2
obj40Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj41)
obj42= Transform(Input=obj41)
obj42.Transform = 'Transform'
obj42.Transform.Translate = [63,66,0]
obj42.Transform.Rotate = [0,0,0]

obj42Display = Show(obj42,renderView1)
obj42Display.Opacity = 0.2
obj42Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
obj43 = Cylinder()

obj43.Resolution = 50
obj43.Height = 160
obj43.Radius = 2.1
obj43.Center = [0., -80,0.]

obj44 = Transform(Input=obj43)
obj44.Transform = 'Transform'
obj44.Transform.Rotate = [90.0, 0.0, 0.0]

obj43Display = Show(obj44,renderView1)
obj43Display.Opacity = 0.2
obj43Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj44)
obj45= Transform(Input=obj44)
obj45.Transform = 'Transform'
obj45.Transform.Translate = [63,78,0]
obj45.Transform.Rotate = [0,0,0]

obj45Display = Show(obj45,renderView1)
obj45Display.Opacity = 0.2
obj45Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#

obj46= GroupDatasets(Input=[obj3,obj6,obj9,obj12,obj15,obj18,obj21,obj24,obj27,obj30,obj33,obj36,obj39,obj42,obj45]) 
obj46Display = Show(obj46, renderView1)
obj46Display.Opacity = 0.2
obj46Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
