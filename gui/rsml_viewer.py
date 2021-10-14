import sys; sys.path.append("../src/python_modules/"); sys.path.append("../")

import rsml_reader
import plantbox as pb
import vtk_plot as vp

import tkinter
from tkinter import ttk
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backend_bases import key_press_handler  # Implement the default Matplotlib key bindings
import matplotlib.pyplot as plt
import numpy as np
from numpy import place

""" gui """
root = None  # tk main frame
# info tab
label_general = None
label_prop = None
label_fun = None
ax = None
canvas = None

""" data """
fname = None
polylines = None
properties = None
functions = None
analyser = None


def open_rsml(fname):
    """ opens an rsml file """
    global polylines
    global properties
    global functions
    polylines, properties, functions = rsml_reader.read_rsml(fname)


def convert_to_analyser():
    """ converts the polylines to a SegmentAnalyzer object """
    global analyser
    if ("parent-node" in properties) and ("parent-poly" in properties):
        nodes, segs = rsml_reader.get_segments(polylines, properties)  # fetch nodes and segments
        radii, cts, types = rsml_reader.get_parameter(polylines, functions, properties)  # paramter per node
        segRadii = np.zeros((segs.shape[0], 1))  # convert to paramter per segment
        segCTs = np.zeros((segs.shape[0], 1))
        segTypes = np.zeros((segs.shape[0], 1))
        for i, s in enumerate(segs):
            segRadii[i] = radii[s[1] - 1]  # seg to node index
            segCTs[i] = cts[s[1] - 1]
            segTypes[i] = types[s[1] - 1]
        segs_ = [pb.Vector2i(s[0], s[1]) for s in segs]  # convert to CPlantBox types
        nodes_ = [pb.Vector3d(n[0], n[1], n[2]) for n in nodes]
        analyser = pb.SegmentAnalyser(nodes_, segs_, segCTs, segRadii)
        analyser.data["subType"] = segTypes
    else:
        tkinter.messagebox.showwarning("Warning", "'parent-node' and or 'parent-poly' property is missing\nCannot create SegmentAnalyzer")


def update_info():
    """ update info tab """
    c = 0
    for pl in polylines:
        c += 1
        for p in pl:
            c += 1
    gstr = "\nFilename \t\t" + fname + "\n"
    gstr += "Number of roots \t\t{:g}\n".format(len(polylines))
    gstr += "Number of nodes \t\t{:g}\n".format(c)
    gstr += "Bounding box \t\t{:s}-{:s}\n".format(str(analyser.getMinBounds()), str(analyser.getMaxBounds()))
    label_general.set(gstr)
    pstr = "\n"
    for k in properties.keys():
        v = np.array(properties[k])
        pstr += k.ljust(30)
        pstr += "\t[{:g}, {:g}]\n".format(np.min(v), np.max(v))
    label_prop.set(pstr)
    fstr = "\n"
    for k in functions.keys():
        v_ = functions[k]
        v = []
        for vv in v_:
            v.extend(vv)
        v = np.array(v)
        fstr += k.ljust(30)
        fstr += "\t[{:g}, {:g}]\n".format(np.min(v), np.max(v))
    label_fun.set(fstr)


def update_profile():
    """ depth profile """
    n = np.ceil(-analyser.getMinBounds().z)
    z_ = np.linspace(-0.5, -n + 0.5, n)
    # print(z_)
    d = analyser.distribution("length", 0., float(-n), int(n), True)
    ax.plot(d, z_, "-*")
    ax.set_ylabel("Depth (cm)")
    ax.set_xlabel("Length per 1 cm layer (cm)")
    canvas.draw()


def update_all():
    """ updates the view """
    if polylines is not None:
        convert_to_analyser()
        update_info()
        update_profile()


def file_open():
    """ file menu item """
    global fname
    fname = tkinter.filedialog.askopenfilename(title = 'Please select a RSML root system',
                                              filetypes = [('Image Files', ['.rsml', '.RSML', '.xml'])])
    if isinstance(fname, str):
        open_rsml(fname)
        update_all()


def view_vtk_plot():
    """ view menu item """
    if analyser is not None:
        vp.plot_roots(analyser, "creationTime")
    else:
        tkinter.messagebox.showwarning("Warning", "No SegmentAnalyser")


def view_vtk_anim():
    """ view menu item """
    print("na")


def view_about():
    """ view menu item """
    tkinter.messagebox.showinfo("About", "RSML Viewer \nby Daniel Leitner, 2021 \n\nPart of CPlantBox")

# def _quit():
#     root.quit()  # stops mainloop
#     root.destroy()  # this is necessary on Windows to prevent
#                     # Fatal Python Error: PyEval_RestoreThread: NULL tstate


root = tkinter.Tk()
root.wm_title("RSML Viewer")
root.geometry("850x800")
menu = tkinter.Menu(root)
menu_file = tkinter.Menu(menu, tearoff = 0)
menu_file.add_command(label = "Open...", command = file_open)
menu.add_cascade(label = "File", menu = menu_file)
menu_view = tkinter.Menu(menu, tearoff = 0)
menu_view.add_command(label = "Vtk plot...", command = view_vtk_plot)
menu_view.add_command(label = "Vtk animation...", command = view_vtk_anim)
menu_view.add_separator()
menu_view.add_command(label = "About...", command = view_about)
menu.add_cascade(label = "View", menu = menu_view)
root.config(menu = menu)

tabControl = ttk.Notebook(root)
tab_info = ttk.Frame(tabControl)
tab_depth = ttk.Frame(tabControl)
tab_development = ttk.Frame(tabControl)
tab_suf = ttk.Frame(tabControl)
tabControl.add(tab_info, text = 'Information')
tabControl.add(tab_depth, text = 'Depth  profile')
tabControl.add(tab_development, text = 'Development')
tabControl.add(tab_suf, text = 'Hydraulic properties')
tabControl.pack(expand = 1, fill = "both")

# tab_info
lf_general = ttk.LabelFrame(tab_info, text = 'General')
lf_general.grid(column = 0, row = 0, padx = 20, pady = 10)
lf_prop = ttk.LabelFrame(tab_info, text = 'Properties (values per root)')
lf_prop.grid(column = 0, row = 1, padx = 20, pady = 10)
lf_fun = ttk.LabelFrame(tab_info, text = 'Functions (values per node)')
lf_fun.grid(column = 0, row = 2, padx = 20, pady = 10)
label_general = tkinter.StringVar()
label1 = ttk.Label(lf_general, textvariable = label_general, anchor = "w", width = 100).pack()
label_prop = tkinter.StringVar()
label2 = ttk.Label(lf_prop, textvariable = label_prop, anchor = "w", width = 100).pack()
label_fun = tkinter.StringVar()
label3 = ttk.Label(lf_fun, textvariable = label_fun, anchor = "w", width = 100).pack()

# tab_profile
fig, ax = plt.subplots(1, 1, figsize = (15, 10))
canvas = FigureCanvasTkAgg(fig, master = tab_depth)  # A tk.DrawingArea.
canvas.draw()
canvas.get_tk_widget().pack(side = tkinter.TOP, fill = tkinter.BOTH, expand = 1)

tkinter.mainloop()
