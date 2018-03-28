import vtk

reader = vtk.vtkXMLPolyDataReader()
path = "example_5a.vtp" #Archive path
reader.SetFileName(path)
reader.Update()
 
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputData(reader.GetOutput())
 
coneActor = vtk.vtkActor()
coneActor.SetMapper( mapper )

ren1= vtk.vtkRenderer()
ren1.AddActor( coneActor )
ren1.SetBackground( 0.1, 0.2, 0.4 )

renWin = vtk.vtkRenderWindow()
renWin.AddRenderer( ren1 )
renWin.SetSize( 600, 600 )


iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)
    
iren.Initialize()
renWin.Render()
iren.Start()
