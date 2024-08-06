from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()
renderView1 = GetActiveViewOrCreate('RenderView')

# GROUP 
obj1 = Box()

obj1.XLength = 200
obj1.YLength = 200
obj1.ZLength = 700
obj1.Center = [0.,0., -350]

obj1Display = Show(obj1,renderView1)
obj1Display.Opacity = 0.2
obj1Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
# GROUP 
# GROUP 
obj2 = Cylinder()

obj2.Resolution = 50
obj2.Height = 700
obj2.Radius = 3.2
obj2.Center = [0., -350,0.]

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
obj4.Transform.Translate = [350,0,0]
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
obj6.Height = 700
obj6.Radius = 3.2
obj6.Center = [0., -350,0.]

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
obj8.Transform.Translate = [350,0,0]
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
obj10.Height = 700
obj10.Radius = 3.2
obj10.Center = [0., -350,0.]

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
obj12.Transform.Translate = [350,0,0]
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
obj14.Height = 700
obj14.Radius = 3.2
obj14.Center = [0., -350,0.]

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
obj16.Transform.Translate = [350,0,0]
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
obj18.Height = 700
obj18.Radius = 3.2
obj18.Center = [0., -350,0.]

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
obj20.Transform.Translate = [350,0,0]
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
obj22.Height = 700
obj22.Radius = 3.2
obj22.Center = [0., -350,0.]

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
obj24.Transform.Translate = [350,0,0]
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

Hide(obj26)
obj27= Transform(Input=obj26)
obj27.Transform = 'Transform'
obj27.Transform.Translate = [0,105,0]
obj27.Transform.Rotate = [0,0,0]

obj27Display = Show(obj27,renderView1)
obj27Display.Opacity = 0.1
obj27Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
# GROUP 
obj28 = Cylinder()

obj28.Resolution = 50
obj28.Height = 700
obj28.Radius = 3.2
obj28.Center = [0., -350,0.]

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
obj30.Transform.Translate = [350,0,0]
obj30.Transform.Rotate = [0,90,0]

obj30Display = Show(obj30,renderView1)
obj30Display.Opacity = 0.1
obj30Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj30)
obj31= Transform(Input=obj30)
obj31.Transform = 'Transform'
obj31.Transform.Translate = [0,-30,-10]
obj31.Transform.Rotate = [0,0,0]

obj31Display = Show(obj31,renderView1)
obj31Display.Opacity = 0.1
obj31Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
obj32 = Cylinder()

obj32.Resolution = 50
obj32.Height = 700
obj32.Radius = 3.2
obj32.Center = [0., -350,0.]

obj33 = Transform(Input=obj32)
obj33.Transform = 'Transform'
obj33.Transform.Rotate = [90.0, 0.0, 0.0]

obj32Display = Show(obj33,renderView1)
obj32Display.Opacity = 0.2
obj32Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj33)
obj34= Transform(Input=obj33)
obj34.Transform = 'Transform'
obj34.Transform.Translate = [350,0,0]
obj34.Transform.Rotate = [0,90,0]

obj34Display = Show(obj34,renderView1)
obj34Display.Opacity = 0.1
obj34Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj34)
obj35= Transform(Input=obj34)
obj35.Transform = 'Transform'
obj35.Transform.Translate = [0,-18,-20]
obj35.Transform.Rotate = [0,0,0]

obj35Display = Show(obj35,renderView1)
obj35Display.Opacity = 0.1
obj35Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
obj36 = Cylinder()

obj36.Resolution = 50
obj36.Height = 700
obj36.Radius = 3.2
obj36.Center = [0., -350,0.]

obj37 = Transform(Input=obj36)
obj37.Transform = 'Transform'
obj37.Transform.Rotate = [90.0, 0.0, 0.0]

obj36Display = Show(obj37,renderView1)
obj36Display.Opacity = 0.2
obj36Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj37)
obj38= Transform(Input=obj37)
obj38.Transform = 'Transform'
obj38.Transform.Translate = [350,0,0]
obj38.Transform.Rotate = [0,90,0]

obj38Display = Show(obj38,renderView1)
obj38Display.Opacity = 0.1
obj38Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj38)
obj39= Transform(Input=obj38)
obj39.Transform = 'Transform'
obj39.Transform.Translate = [0,-6,-40]
obj39.Transform.Rotate = [0,0,0]

obj39Display = Show(obj39,renderView1)
obj39Display.Opacity = 0.1
obj39Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
obj40 = Cylinder()

obj40.Resolution = 50
obj40.Height = 700
obj40.Radius = 3.2
obj40.Center = [0., -350,0.]

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
obj42.Transform.Translate = [350,0,0]
obj42.Transform.Rotate = [0,90,0]

obj42Display = Show(obj42,renderView1)
obj42Display.Opacity = 0.1
obj42Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj42)
obj43= Transform(Input=obj42)
obj43.Transform = 'Transform'
obj43.Transform.Translate = [0,6,-60]
obj43.Transform.Rotate = [0,0,0]

obj43Display = Show(obj43,renderView1)
obj43Display.Opacity = 0.1
obj43Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
obj44 = Cylinder()

obj44.Resolution = 50
obj44.Height = 700
obj44.Radius = 3.2
obj44.Center = [0., -350,0.]

obj45 = Transform(Input=obj44)
obj45.Transform = 'Transform'
obj45.Transform.Rotate = [90.0, 0.0, 0.0]

obj44Display = Show(obj45,renderView1)
obj44Display.Opacity = 0.2
obj44Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj45)
obj46= Transform(Input=obj45)
obj46.Transform = 'Transform'
obj46.Transform.Translate = [350,0,0]
obj46.Transform.Rotate = [0,90,0]

obj46Display = Show(obj46,renderView1)
obj46Display.Opacity = 0.1
obj46Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj46)
obj47= Transform(Input=obj46)
obj47.Transform = 'Transform'
obj47.Transform.Translate = [0,18,-80]
obj47.Transform.Rotate = [0,0,0]

obj47Display = Show(obj47,renderView1)
obj47Display.Opacity = 0.1
obj47Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
obj48 = Cylinder()

obj48.Resolution = 50
obj48.Height = 700
obj48.Radius = 3.2
obj48.Center = [0., -350,0.]

obj49 = Transform(Input=obj48)
obj49.Transform = 'Transform'
obj49.Transform.Rotate = [90.0, 0.0, 0.0]

obj48Display = Show(obj49,renderView1)
obj48Display.Opacity = 0.2
obj48Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj49)
obj50= Transform(Input=obj49)
obj50.Transform = 'Transform'
obj50.Transform.Translate = [350,0,0]
obj50.Transform.Rotate = [0,90,0]

obj50Display = Show(obj50,renderView1)
obj50Display.Opacity = 0.1
obj50Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj50)
obj51= Transform(Input=obj50)
obj51.Transform = 'Transform'
obj51.Transform.Translate = [0,30,-120]
obj51.Transform.Rotate = [0,0,0]

obj51Display = Show(obj51,renderView1)
obj51Display.Opacity = 0.1
obj51Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#

obj52= GroupDatasets(Input=[obj31,obj35,obj39,obj43,obj47,obj51]) 
obj52Display = Show(obj52, renderView1)
obj52Display.Opacity = 0.1
obj52Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
# GROUP 
obj53 = Cylinder()

obj53.Resolution = 50
obj53.Height = 700
obj53.Radius = 3.2
obj53.Center = [0., -350,0.]

obj54 = Transform(Input=obj53)
obj54.Transform = 'Transform'
obj54.Transform.Rotate = [90.0, 0.0, 0.0]

obj53Display = Show(obj54,renderView1)
obj53Display.Opacity = 0.2
obj53Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj54)
obj55= Transform(Input=obj54)
obj55.Transform = 'Transform'
obj55.Transform.Translate = [350,0,0]
obj55.Transform.Rotate = [0,90,0]

obj55Display = Show(obj55,renderView1)
obj55Display.Opacity = 0.1
obj55Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj55)
obj56= Transform(Input=obj55)
obj56.Transform = 'Transform'
obj56.Transform.Translate = [0,-30,-10]
obj56.Transform.Rotate = [0,0,0]

obj56Display = Show(obj56,renderView1)
obj56Display.Opacity = 0.1
obj56Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
obj57 = Cylinder()

obj57.Resolution = 50
obj57.Height = 700
obj57.Radius = 3.2
obj57.Center = [0., -350,0.]

obj58 = Transform(Input=obj57)
obj58.Transform = 'Transform'
obj58.Transform.Rotate = [90.0, 0.0, 0.0]

obj57Display = Show(obj58,renderView1)
obj57Display.Opacity = 0.2
obj57Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj58)
obj59= Transform(Input=obj58)
obj59.Transform = 'Transform'
obj59.Transform.Translate = [350,0,0]
obj59.Transform.Rotate = [0,90,0]

obj59Display = Show(obj59,renderView1)
obj59Display.Opacity = 0.1
obj59Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj59)
obj60= Transform(Input=obj59)
obj60.Transform = 'Transform'
obj60.Transform.Translate = [0,-18,-20]
obj60.Transform.Rotate = [0,0,0]

obj60Display = Show(obj60,renderView1)
obj60Display.Opacity = 0.1
obj60Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
obj61 = Cylinder()

obj61.Resolution = 50
obj61.Height = 700
obj61.Radius = 3.2
obj61.Center = [0., -350,0.]

obj62 = Transform(Input=obj61)
obj62.Transform = 'Transform'
obj62.Transform.Rotate = [90.0, 0.0, 0.0]

obj61Display = Show(obj62,renderView1)
obj61Display.Opacity = 0.2
obj61Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj62)
obj63= Transform(Input=obj62)
obj63.Transform = 'Transform'
obj63.Transform.Translate = [350,0,0]
obj63.Transform.Rotate = [0,90,0]

obj63Display = Show(obj63,renderView1)
obj63Display.Opacity = 0.1
obj63Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj63)
obj64= Transform(Input=obj63)
obj64.Transform = 'Transform'
obj64.Transform.Translate = [0,-6,-40]
obj64.Transform.Rotate = [0,0,0]

obj64Display = Show(obj64,renderView1)
obj64Display.Opacity = 0.1
obj64Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
obj65 = Cylinder()

obj65.Resolution = 50
obj65.Height = 700
obj65.Radius = 3.2
obj65.Center = [0., -350,0.]

obj66 = Transform(Input=obj65)
obj66.Transform = 'Transform'
obj66.Transform.Rotate = [90.0, 0.0, 0.0]

obj65Display = Show(obj66,renderView1)
obj65Display.Opacity = 0.2
obj65Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj66)
obj67= Transform(Input=obj66)
obj67.Transform = 'Transform'
obj67.Transform.Translate = [350,0,0]
obj67.Transform.Rotate = [0,90,0]

obj67Display = Show(obj67,renderView1)
obj67Display.Opacity = 0.1
obj67Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj67)
obj68= Transform(Input=obj67)
obj68.Transform = 'Transform'
obj68.Transform.Translate = [0,6,-60]
obj68.Transform.Rotate = [0,0,0]

obj68Display = Show(obj68,renderView1)
obj68Display.Opacity = 0.1
obj68Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
obj69 = Cylinder()

obj69.Resolution = 50
obj69.Height = 700
obj69.Radius = 3.2
obj69.Center = [0., -350,0.]

obj70 = Transform(Input=obj69)
obj70.Transform = 'Transform'
obj70.Transform.Rotate = [90.0, 0.0, 0.0]

obj69Display = Show(obj70,renderView1)
obj69Display.Opacity = 0.2
obj69Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj70)
obj71= Transform(Input=obj70)
obj71.Transform = 'Transform'
obj71.Transform.Translate = [350,0,0]
obj71.Transform.Rotate = [0,90,0]

obj71Display = Show(obj71,renderView1)
obj71Display.Opacity = 0.1
obj71Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj71)
obj72= Transform(Input=obj71)
obj72.Transform = 'Transform'
obj72.Transform.Translate = [0,18,-80]
obj72.Transform.Rotate = [0,0,0]

obj72Display = Show(obj72,renderView1)
obj72Display.Opacity = 0.1
obj72Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#
obj73 = Cylinder()

obj73.Resolution = 50
obj73.Height = 700
obj73.Radius = 3.2
obj73.Center = [0., -350,0.]

obj74 = Transform(Input=obj73)
obj74.Transform = 'Transform'
obj74.Transform.Rotate = [90.0, 0.0, 0.0]

obj73Display = Show(obj74,renderView1)
obj73Display.Opacity = 0.2
obj73Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj74)
obj75= Transform(Input=obj74)
obj75.Transform = 'Transform'
obj75.Transform.Translate = [350,0,0]
obj75.Transform.Rotate = [0,90,0]

obj75Display = Show(obj75,renderView1)
obj75Display.Opacity = 0.1
obj75Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj75)
obj76= Transform(Input=obj75)
obj76.Transform = 'Transform'
obj76.Transform.Translate = [0,30,-120]
obj76.Transform.Rotate = [0,0,0]

obj76Display = Show(obj76,renderView1)
obj76Display.Opacity = 0.1
obj76Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#

obj77= GroupDatasets(Input=[obj56,obj60,obj64,obj68,obj72,obj76]) 
obj77Display = Show(obj77, renderView1)
obj77Display.Opacity = 0.1
obj77Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()

Hide(obj77)
obj78= Transform(Input=obj77)
obj78.Transform = 'Transform'
obj78.Transform.Translate = [0,-105,0]
obj78.Transform.Rotate = [0,0,0]

obj78Display = Show(obj78,renderView1)
obj78Display.Opacity = 0.1
obj78Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#

obj79= GroupDatasets(Input=[obj27,obj52,obj78]) 
obj79Display = Show(obj79, renderView1)
obj79Display.Opacity = 0.1
obj79Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
#

obj80= GroupDatasets(Input=[obj1,obj79]) 
obj80Display = Show(obj80, renderView1)
obj80Display.Opacity = 0.1
obj80Display.DiffuseColor = [0., 0., 1.0]
renderView1.ResetCamera()
