""" VTK and tkinter does not work """

# Import the modules for the code
import tkinter
import sys
import vtk
# from vtk.tk.vtkTkRenderWindowInteractor import vtkTkRenderWindowInteractor

# Setup for root window
root = tkinter.Tk()
root.title("Tkinter Test")
frame = tkinter.Frame(root)
frame.pack(fill = tkinter.BOTH, expand = 1, side = tkinter.TOP)

# Setup for renderer
render = vtk.vtkRenderer()
render.SetBackground(0.329412, 0.34902, 0.427451)
render.ResetCameraClippingRange()

# Setup for rendering window
renWindow = vtk.vtkRenderWindow()
renWindow.AddRenderer(render)
# Setup for rendering window interactor
renWinInteract = vtk.tk.vtkTkRenderWindowInteractor(root, rw = renWindow, width = 400, height = 400)
renWinInteract.Initialize()
renWinInteract.pack(side = 'top', fill = 'both', expand = 1)
renWinInteract.Start()

# Begin execution by updating the renderer and
# starting the Tkinter loop
renWindow.Render()
root.mainloop()
