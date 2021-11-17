import sys; sys.path.append("../python/modules/"); sys.path.append("../../CPlantBox/");
sys.path.append("../../CPlantBox/src/python_modules")

import vtk_plot as vp
import vtk_tools as vt
import rsml_writer

import tkinter as tk


def hello():
    print("hello")


frame = tk.Tk()
frame.geometry("200x200")
label = tk.Label(frame, text = "Hallo Welt!")
label.grid(row = 1, column = 1)
button = tk.Button(frame, text = "OK", command = hello)
button.grid(row = 2, column = 1)
frame.mainloop()
