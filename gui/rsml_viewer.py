import sys; sys.path.append("../src/python_modules/"); sys.path.append("../")

import rsml_reader
import plantbox as pb
import vtk_plot as vp
import xylem_flux
import viewer_conductivities

import tkinter
from tkinter import ttk
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backend_bases import key_press_handler  # Implement the default Matplotlib key bindings
import matplotlib.pyplot as plt
import numpy as np
from numpy import place

""" gui """
root = None  # tk main frame
label_general_l, label_general_r = None, None  # info tab
label_prop_l, label_prop_r = None, None
label_fun_l, label_fun_r = None, None
type_str = ["length", "surface", "volume"]
unit_str = ["(cm)", "(cm$^2$)", "(cm$^3$)"]
ax = None  # depth profile
canvas = None
ax2 = None  # development
canvas2 = None
ax3 = None  # hydraulics
canvas3 = None

""" data """
fname = None
polylines, properties, functions, metadata = None, None, None, None  # open_rsml
radii, cts, types = None, None, None
analyser = None  # convert_to_analyser
max_ct = 0.
mapped_segments = None


def open_rsml(fname):
    """ opens an rsml file """
    global polylines, properties, functions, metadata
    global radii, cts, types
    polylines, properties, functions, metadata = rsml_reader.read_rsml(fname)
    radii, cts, types = rsml_reader.get_parameter(polylines, functions, properties)  # paramter per node


def scale_polylines_(scale):
    """ scales coordinates """
    global polylines
    if polylines is not None:
        for i in range(0, len(polylines)):
            for j in range(0, len(polylines[i])):
                for k in range(0, 3):
                    polylines[i][j][k] *= scale
    else:
        tkinter.messagebox.showwarning("Warning", "Open RSML file first")


def convert_to_analyser():
    """ converts the polylines to a SegmentAnalyser object,
        converts the SegmentAnalyser to a mapped segments object """
    global analyser, max_ct, mapped_segments
    global cts
    global subTypes
    if "parent-poly" in properties:
        scale_polylines_(metadata.scale_to_cm)
        nodes, segs = rsml_reader.get_segments(polylines, properties)  # fetch nodes and segments
        segRadii = np.zeros((segs.shape[0], 1))  # convert to paramter per segment
        segCTs = np.zeros((segs.shape[0], 1))
        subTypes = np.zeros((segs.shape[0], 1))
        for i, s in enumerate(segs):
            segRadii[i] = radii[s[1] - 1]  # seg to node index
            segCTs[i] = cts[s[1] - 1]
            subTypes[i] = types[s[1] - 1]
        max_ct = np.max(segCTs)
        segs_ = [pb.Vector2i(s[0], s[1]) for s in segs]  # convert to CPlantBox types
        nodes_ = [pb.Vector3d(n[0], n[1], n[2]) for n in nodes]
        analyser = pb.SegmentAnalyser(nodes_, segs_, segCTs, segRadii)
        analyser.addData("subType", subTypes)
        if np.isnan(cts[0]):
            cts = np.zeros((len(nodes_),))
        if np.isnan(subTypes[0]):
            subTypes = np.ones((len(segs_),), dtype=np.int64)
        mapped_segments = pb.MappedSegments(nodes_, np.array(cts), segs_, np.array(segRadii), np.array(subTypes))
    else:
        print(properties.keys())
        tkinter.messagebox.showwarning("Warning", "'parent-node' and or 'parent-poly' property is missing, cannot create SegmentAnalyser")


def update_info():
    """ update info tab """
    c = 0
    for pl in polylines:
        c += 1
        for p in pl:
            c += 1
    lstr = "\nFilename \nNumber of roots \nNumber of nodes \nBounding box\n"
    rstr = "\n" + fname + "\n{:g}\n{:g}\n{:s}-{:s}\n".format(len(polylines), c, str(analyser.getMinBounds()), str(analyser.getMaxBounds()))
    label_general_l.set(lstr)
    label_general_r.set(rstr)
    lstr, rstr = "\n", "\n"
    for k in properties.keys():
        v = np.array(properties[k])
        lstr += k + "\n"
        rstr += "[{:g}, {:g}]\n".format(np.min(v), np.max(v))
    label_prop_l.set(lstr)
    label_prop_r.set(rstr)
    lstr, rstr = "\n", "\n"
    for k in functions.keys():
        v_ = functions[k]
        v = []
        for vv in v_:
            v.extend(vv)
        v = np.array(v)
        lstr += k + "\n"
        rstr += "[{:g}, {:g}]\n".format(np.min(v), np.max(v))
    label_fun_l.set(lstr)
    label_fun_r.set(rstr)


def update_profile():
    """ depth profile """
    j = 0  # pick length, surface, volume
    ax.clear()
    n = int(np.ceil(-analyser.getMinBounds().z))
    z_ = np.linspace(-0.5, -n + 0.5, n)
    d = analyser.distribution(type_str[j], 0., float(-n), int(n), True)
    ax.plot(d, z_, "-*", label="total")
    max_type = int(np.max(analyser.data["subType"]))
    for i in range(0, max_type + 1):
        ana = pb.SegmentAnalyser(analyser)  # copy
        ana.filter("subType", i)
        segn = len(ana.segments)
        if segn > 0:
            d = ana.distribution(type_str[j], 0., float(-n), int(n), True)
            ax.plot(d, z_, "-*", label="type {:g}".format(i))
    ax.set_ylabel("Depth (cm)")
    ax.set_xlabel("Root system " + type_str[j] + " per 1 cm layer " + unit_str[j])
    ax.legend()
    canvas.draw()


def update_development():
    """ root system development """
    j = 0
    ax2.clear()
    weights = [analyser.getSegmentLength(i) for i in range(0, len(analyser.segments))]
    cts = np.array(analyser.data["creationTime"])
    l_, t_ = np.histogram(cts, 100, weights=weights)
    ax2.plot(0.5 * (t_[1:] + t_[:-1]), np.cumsum(l_), "-", label="total")
    max_type = int(np.max(analyser.data["subType"]))
    for i in range(0, max_type + 1):
        ana = pb.SegmentAnalyser(analyser)  # copy
        ana.filter("subType", i)
        n = len(ana.segments)
        if n > 0:
            weights = [ana.getSegmentLength(ii) for ii in range(0, n)]
            cts = np.array(ana.data["creationTime"])
            l_, t_ = np.histogram(cts, 100, weights=weights)
            ax2.plot(0.5 * (t_[1:] + t_[:-1]), np.cumsum(l_), "-", label="type {:g}".format(i))
    ax2.set_xlabel("Time (days)")
    ax2.set_ylabel("Root system " + type_str[j] + " " + unit_str[j])
    ax2.legend()
    canvas2.draw()


def update_hydraulics():
    global analyser
    """ hydraulic properties """
    ax3.clear()
    r = xylem_flux.XylemFluxPython(mapped_segments)
    viewer_conductivities.init_conductivities_scenario_jan_const(r)  # TODO pick scenario
    krs, _ = r.get_krs(max_ct)  # TODO move to label
    suf = r.get_suf(max_ct)
    analyser.addData("SUF", suf)
    n = int(np.ceil(-analyser.getMinBounds().z))
    z_ = np.linspace(-0.5, -n + 0.5, n)
    d = analyser.distribution("SUF", 0., float(-n), int(n), True)
    ax3.plot(d, z_, "-*", label="total")
    max_type = int(np.max(analyser.data["subType"]))
    for i in range(0, max_type + 1):
        ana = pb.SegmentAnalyser(analyser)  # copy
        ana.filter("subType", i)
        segn = len(ana.segments)
        if segn > 0:
            d = ana.distribution("SUF", 0., float(-n), int(n), True)
            ax3.plot(d, z_, "-*", label="type {:g}".format(i))
    ax3.set_title("Root system krs {:g}".format(krs))
    ax3.set_ylabel("Depth (cm)")
    ax3.set_xlabel("Root system surface uptake fraction (SUF) (1)")
    ax3.legend()
    canvas3.draw()


def update_all():
    """ updates the view """
    if polylines is not None:
        convert_to_analyser()
        update_info()
        update_profile()
        update_development()
        update_hydraulics()


def file_open():
    """ file menu item """
    global fname
    fname = tkinter.filedialog.askopenfilename(title='Please select a RSML root system',
                                              filetypes=[('Image Files', ['.rsml', '.RSML', '.xml'])])
    if isinstance(fname, str):
        if fname:
            open_rsml(fname)
            update_all()


def view_vtk_plot(name):
    """ view menu item """
    if analyser is not None:
        vp.plot_roots(analyser, name)
    else:
        tkinter.messagebox.showwarning("Warning", "Open RSML file first")


def view_vtk_plot_subtype():
    view_vtk_plot(("subType"))


def view_vtk_plot_creationtime():
    view_vtk_plot("creationTime")


def view_vtk_plot_suf():
    view_vtk_plot("SUF")


def view_vtk_anim():
    """ view menu item """
    print("na")


def view_about():
    """ view menu item """
    tkinter.messagebox.showinfo("About", "RSML Viewer \nby Daniel Leitner, 2021 \n\nPart of CPlantBox")


def _quit():
    root.quit()  # stops mainloop
    root.destroy()  # this is necessary on Windows to prevent
                    # Fatal Python Error: PyEval_RestoreThread: NULL tstate


root = tkinter.Tk()
root.wm_title("RSML Viewer")
root.geometry("850x800")
menu = tkinter.Menu(root)
menu_file = tkinter.Menu(menu, tearoff=0)
menu_file.add_command(label="Open...", command=file_open)
menu_file.add_separator()
menu_file.add_command(label="Exit", command=_quit)
menu.add_cascade(label="File", menu=menu_file)

menu_view = tkinter.Menu(menu, tearoff=0)
menu_view.add_command(label="Type...", command=view_vtk_plot_subtype)
menu_view.add_command(label="Creation Time...", command=view_vtk_plot_creationtime)
menu_view.add_command(label="SUF...", command=view_vtk_plot_suf)
menu_view.add_command(label="Animation...", command=view_vtk_anim)
menu_view.add_separator()
menu_view.add_command(label="About...", command=view_about)
menu.add_cascade(label="View", menu=menu_view)
root.config(menu=menu)

tabControl = ttk.Notebook(root)
tab_info = ttk.Frame(tabControl)
tab_depth = ttk.Frame(tabControl)
tab_development = ttk.Frame(tabControl)
tab_suf = ttk.Frame(tabControl)
tabControl.add(tab_info, text='Information')
tabControl.add(tab_depth, text='Depth  profile')
tabControl.add(tab_development, text='Development')
tabControl.add(tab_suf, text='Hydraulic properties')
tabControl.pack(expand=1, fill="both")

# tab_info
lf_general = ttk.LabelFrame(tab_info, text='General')
lf_general.grid(column=0, row=0, padx=20, pady=10)
lf_prop = ttk.LabelFrame(tab_info, text='Properties (values per root)')
lf_prop.grid(column=0, row=1, padx=20, pady=10)
lf_fun = ttk.LabelFrame(tab_info, text='Functions (values per node)')
lf_fun.grid(column=0, row=2, padx=20, pady=10)
label_general_l = tkinter.StringVar()
label_general_r = tkinter.StringVar()
ttk.Label(lf_general, textvariable=label_general_l, anchor="w", width=30).grid(column=0, row=0)
ttk.Label(lf_general, textvariable=label_general_r, anchor="w", width=70).grid(column=1, row=0)
label_prop_l = tkinter.StringVar()
label_prop_r = tkinter.StringVar()
ttk.Label(lf_prop, textvariable=label_prop_l, anchor="w", width=30).grid(column=0, row=0)
ttk.Label(lf_prop, textvariable=label_prop_r, anchor="w", width=70).grid(column=1, row=0)
label_fun_l = tkinter.StringVar()
label_fun_r = tkinter.StringVar()
ttk.Label(lf_fun, textvariable=label_fun_l, anchor="w", width=30).grid(column=0, row=0)
ttk.Label(lf_fun, textvariable=label_fun_r, anchor="w", width=70).grid(column=1, row=0)

# tab_profile
fig, ax = plt.subplots(1, 1, figsize=(15, 10))
canvas = FigureCanvasTkAgg(fig, master=tab_depth)  # A tk.DrawingArea.
canvas.draw()
canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)

# tab_development
fig2, ax2 = plt.subplots(1, 1, figsize=(15, 10))
canvas2 = FigureCanvasTkAgg(fig2, master=tab_development)  # A tk.DrawingArea.
canvas2.draw()
canvas2.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)

# hydraulic properties
fig3, ax3 = plt.subplots(1, 1, figsize=(15, 10))
canvas3 = FigureCanvasTkAgg(fig3, master=tab_suf)  # A tk.DrawingArea.
canvas3.draw()
canvas3.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)

tkinter.mainloop()
