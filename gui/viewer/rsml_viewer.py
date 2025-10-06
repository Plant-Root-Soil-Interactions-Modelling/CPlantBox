import sys; sys.path.append("../.."); sys.path.append("../../src/")

import functional.xylem_flux as xylem_flux
import visualisation.vtk_plot as vp
import visualisation.vtk_tools as vt
from viewer_data import ViewerDataModel
import viewer_plots
import viewer_conductivities

import tkinter
from tkinter import ttk
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backend_bases import key_press_handler  # Implement the default Matplotlib key bindings
import matplotlib.pyplot as plt
import numpy as np

""" TODO animation, and optional remesh and RSML (re)write would be nice, writing artifical shoot(?!) """
""" log file """


class App:

    def __init__(self, root):
        self.data = ViewerDataModel()
        self.root = root
        self.root.wm_title("RSML Viewer")
        self.root.geometry("850x800")
        # Menu
        menu = tkinter.Menu(root)
        menu_file = tkinter.Menu(menu, tearoff = 0)
        menu_file.add_command(label = "Open (.rsml)...", command = self.file_open)
        menu_file.add_command(label = "Save (.rsml)...", command = self.file_save)
        menu_file.add_command(label = "Save (.vtp)...", command = self.file_save_vtp)
        menu_file.add_separator()
        menu_file.add_command(label = "Exit", command = self.file_quit)
        menu.add_cascade(label = "File", menu = menu_file)
        menu_edit = tkinter.Menu(menu, tearoff = 0)
        menu_edit.add_command(label = "Add shoot", command = self.edit_add_shoot)
        menu_edit.add_command(label = "Add creation time", command = self.edit_add_creation_times)
        menu.add_cascade(label = "Edit", menu = menu_edit)
        menu_view = tkinter.Menu(menu, tearoff = 0)
        menu_view.add_command(label = "Type...", command = self.view_vtk_plot_subtype)
        menu_view.add_command(label = "Segment length...", command = self.view_vtk_plot_length)
        menu_view.add_command(label = "Creation time...", command = self.view_vtk_plot_creationtime)
        menu_view.add_command(label = "SUF...", command = self.view_vtk_plot_suf)
        menu_view.add_command(label = "Multiple SUF...", command = self.view_vtk_plot_multiple_suf)
        menu_view.add_command(label = "Animation...", command = self.view_vtk_anim)
        menu_view.add_separator()
        menu_view.add_command(label = "About...", command = self.view_about)
        menu.add_cascade(label = "View", menu = menu_view)
        self.root.config(menu = menu)
        # Tabs
        tabControl = ttk.Notebook(self.root)
        tab_info = ttk.Frame(tabControl)
        tab_depth = ttk.Frame(tabControl)
        tab_development = ttk.Frame(tabControl)
        tab_suf = ttk.Frame(tabControl)
        tab_krs = ttk.Frame(tabControl)
        tabControl.add(tab_info, text = 'Information')
        tabControl.add(tab_depth, text = 'Root depth  profile')
        tabControl.add(tab_development, text = 'Root development')
        tabControl.add(tab_suf, text = 'Hydraulic properties')
        tabControl.add(tab_krs, text = 'Hydraulic development')
        tabControl.pack(expand = 1, fill = "both")
        # tab_info
        lf_general = ttk.LabelFrame(tab_info, text = 'General')
        lf_general.grid(column = 0, row = 0, padx = 20, pady = 10)
        lf_prop = ttk.LabelFrame(tab_info, text = 'Properties (values per root)')
        lf_prop.grid(column = 0, row = 1, padx = 20, pady = 10)
        lf_fun = ttk.LabelFrame(tab_info, text = 'Functions (values per node)')
        lf_fun.grid(column = 0, row = 2, padx = 20, pady = 10)
        lf_use = ttk.LabelFrame(tab_info, text = 'Using')
        lf_use.grid(column = 0, row = 3, padx = 20, pady = 10)
        self.label_general_l = tkinter.StringVar()
        self.label_general_r = tkinter.StringVar()
        ttk.Label(lf_general, textvariable = self.label_general_l, anchor = "w", width = 30).grid(column = 0, row = 0)
        ttk.Label(lf_general, textvariable = self.label_general_r, anchor = "w", width = 70).grid(column = 1, row = 0)
        self.label_prop_l = tkinter.StringVar()
        self.label_prop_r = tkinter.StringVar()
        ttk.Label(lf_prop, textvariable = self.label_prop_l, anchor = "w", width = 30).grid(column = 0, row = 0)
        ttk.Label(lf_prop, textvariable = self.label_prop_r, anchor = "w", width = 70).grid(column = 1, row = 0)
        self.label_fun_l = tkinter.StringVar()
        self.label_fun_r = tkinter.StringVar()
        ttk.Label(lf_fun, textvariable = self.label_fun_l, anchor = "w", width = 30).grid(column = 0, row = 0)
        ttk.Label(lf_fun, textvariable = self.label_fun_r, anchor = "w", width = 70).grid(column = 1, row = 0)
        self.label_use_l = tkinter.StringVar()
        self.label_use_r = tkinter.StringVar()
        ttk.Label(lf_use, textvariable = self.label_use_l, anchor = "w", width = 30).grid(column = 0, row = 0)
        ttk.Label(lf_use, textvariable = self.label_use_r, anchor = "w", width = 70).grid(column = 1, row = 0)
        # tab_profile
        self.combo1 = ttk.Combobox(tab_depth, values = [ "Length", "Surface", "Volume"])
        self.combo1.pack(pady = 10)
        self.combo1.current(0)
        self.combo1.bind("<<ComboboxSelected>>", self.update_profile)
        fig, self.ax = plt.subplots(1, 1, figsize = (7, 7))
        self.canvas = FigureCanvasTkAgg(fig, master = tab_depth)  # A tk.DrawingArea.
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill = tkinter.BOTH, expand = 1)
        # tab_development
        self.combo2 = ttk.Combobox(tab_development, values = [ "Length", "Surface", "Volume"])
        self.combo2.pack(pady = 10)
        self.combo2.current(0)
        self.combo2.bind("<<ComboboxSelected>>", self.update_development)
        fig2, self.ax2 = plt.subplots(1, 1, figsize = (15, 10))
        self.canvas2 = FigureCanvasTkAgg(fig2, master = tab_development)  # A tk.DrawingArea.
        self.canvas2.draw()
        self.canvas2.get_tk_widget().pack(side = tkinter.TOP, fill = tkinter.BOTH, expand = 1)
        # hydraulic properties
        tab_suf_frame = tkinter.Frame(tab_suf)
        self.combo3 = ttk.Combobox(tab_suf_frame, values = [ "Constant scenario 1", "Constant scenario 2", "Dynamic scenario 1", "Dynamic scenario 2", "Wine"])
        self.combo3.pack(side = tkinter.LEFT, pady = 5, padx = 10)
        self.combo3.current(0)
        self.combo3.bind("<<ComboboxSelected>>", self.update_hydraulics)
        button = tkinter.Button(tab_suf_frame, text = "plot conductivities", command = self.plot_conductivities)
        button.pack(side = tkinter.LEFT, pady = 5, padx = 10)
        tab_suf_frame.pack(side = tkinter.TOP)
        fig3, self.ax3 = plt.subplots(1, 1, figsize = (15, 10))
        self.canvas3 = FigureCanvasTkAgg(fig3, master = tab_suf)  # A tk.DrawingArea.
        self.canvas3.draw()
        self.canvas3.get_tk_widget().pack(side = tkinter.BOTTOM, fill = tkinter.BOTH, expand = 1)
        # hydraulic development
        tab_krs_frame = tkinter.Frame(tab_krs)
        self.combo4 = ttk.Combobox(tab_krs_frame, values = [ "Constant scenario 1", "Constant scenario 2", "Dynamic scenario 1", "Dynamic scenario 2", "Wine"])
        self.combo4.pack(side = tkinter.LEFT, pady = 5, padx = 10)
        self.combo4.current(0)
        self.combo4.bind("<<ComboboxSelected>>", self.update_krs)
        button2 = tkinter.Button(tab_krs_frame, text = "plot conductivities", command = self.plot_conductivities)
        button2.pack(side = tkinter.LEFT, pady = 5, padx = 10)
        tab_krs_frame.pack(side = tkinter.TOP)
        fig4, self.ax4 = plt.subplots(1, 1, figsize = (15, 10))
        self.canvas4 = FigureCanvasTkAgg(fig4, master = tab_krs)  # A tk.DrawingArea.
        self.canvas4.draw()
        self.canvas4.get_tk_widget().pack(side = tkinter.BOTTOM, fill = tkinter.BOTH, expand = 1)

    def update_info(self):
        """ update info tab """
        # label_general
        c = 0  # node counter
        for i, pl in enumerate(self.data.polylines):
            # c += 1
            for p in pl:
                c += 1
        lstr = "\nSoftware\nFilename \nNumber of plants (base nodes)\nNumber of base roots (segments)\nNumber of roots\nNumber of nodes\n"
        lstr += "Bounding box \nUnit (length scale)\nResolution\n"
        min_str = str(self.data.analyser.getMinBounds())
        max_str = str(self.data.analyser.getMaxBounds())
        metadata = self.data.metadata
        nop = len(self.data.base_nodes)
        nobr = len(self.data.base_segs)
        rstr = "\n{:s}\n{:s}\n{:g}\n{:g}\n{:g}\n{:g}\n".format(metadata.software, fname, nop, nobr, len(self.data.polylines), c)
        rstr += "{:s} - {:s}\n{:s}\n{:s} (dots per {:s})\n".format(min_str, max_str + " (cm)", metadata.unit, str(metadata.resolution), metadata.unit)
        self.label_general_l.set(lstr)
        self.label_general_r.set(rstr)
        # label_prop
        lstr, rstr = "\n", "\n"
        for k in self.data.properties.keys():
            v = np.array(self.data.properties[k])
            lstr += k + "\n"
            rstr += "[{:g}, {:g}]".format(np.min(v), np.max(v))
            if k in metadata.properties:
                rstr += " (" + metadata.properties[k].unit + ")\n"
            else:
                rstr += "\n"
        self.label_prop_l.set(lstr)
        self.label_prop_r.set(rstr)
        # label_fun
        lstr, rstr = "\n", "\n"
        for k in self.data.functions.keys():
            v_ = self.data.functions[k]
            v = []
            for vv in v_:
                v.extend(vv)
            v = np.array(v)
            lstr += k + "\n"
            rstr += "[{:g}, {:g}]".format(np.min(v), np.max(v))
            if k in metadata.properties:
                rstr += " (" + metadata.properties[k].unit + ")\n"
            else:
                rstr += "\n"
        self.label_fun_l.set(lstr)
        self.label_fun_r.set(rstr)
        # label_using
        tagnames = self.data.tagnames
        rstr = "\n"
        lstr = "\nRadius \nCreation time \nTypes \n"
        if tagnames[0]:
            rstr += "from tag '{:s}' within [{:g}, {:g}] cm\n".format(tagnames[0], np.min(self.data.radii), np.max(self.data.radii))
        else:
            rstr += "not found (set to 0.1 cm) \n"
        if tagnames[1]:
            rstr += "from tag '{:s}' within [{:g}, {:g}] days\n".format(tagnames[1], np.min(self.data.cts), np.max(self.data.cts))
        else:
            rstr += "not found\n"
        if tagnames[2]:
            rstr += "from tag '{:s}' within [{:g}, {:g}]\n".format(tagnames[2], np.min(self.data.types), np.max(self.data.types))
        else:
            rstr += "not found\n"  # , derived from root order
        self.label_use_l.set(lstr)
        self.label_use_r.set(rstr)

    def update_profile(self, event):
        """ updates depth profile plot """
        if self.data.exists():
            viewer_plots.plot_depth_profile(self.data.analyser, self.ax, self.combo1.current())
        self.canvas.draw()

    def update_development(self, event):
        viewer_plots.plot_rootsystem_development(self.data.analyser, self.ax2, self.combo2.current())
        self.canvas2.draw()

    def update_hydraulics(self, event):
        """ updates hydraulic properties plot """
        if self.data.exists():
            viewer_plots.plot_suf(self.data, self.ax3, self.combo3.current())
            self.canvas3.draw()

    def update_krs(self, event):
        """ updates hydraulic properties plot """
        if self.data.exists():
            viewer_plots.plot_krs(self.data, self.ax4, self.combo4.current())
            self.canvas4.draw()

    def view_vtk_plot(self, name):
        """ vtk plot coloring name """
        if self.data.exists():
            vp.plot_roots(self.data.analyser, name)
        else:
            tkinter.messagebox.showwarning("Warning", "Open RSML file first")

    def update_all(self):
        """ updates the view """
        if self.data.exists():
            self.update_info()
            self.update_profile(None)
            self.update_development(None)
            self.update_hydraulics(None)
            self.update_krs(None)

    def plot_conductivities(self):
        """ button """
        if self.data.exists():
            j = self.combo3.current()
            root = tkinter.Tk()
            root.wm_title("Root conductivities")
            r = self.data.xylem_flux  # rename
            if j == 0:
                viewer_conductivities.init_constant_scenario1(r)
            elif j == 1:
                viewer_conductivities.init_constant_scenario2(r)
            elif j == 2:
                viewer_conductivities.init_dynamic_scenario1(r)
            elif j == 3:
                viewer_conductivities.init_dynamic_scenario2(r)
            elif j == 4:
                viewer_conductivities.init_constant_scenario_wine(r)

            fig = r.plot_conductivities(monocot = False, plot_now = False, axes_ind = [0], lateral_ind = [1, 2, 3])
            canvas = FigureCanvasTkAgg(fig, master = root)  # A tk.DrawingArea.
            canvas.draw()
            canvas.get_tk_widget().pack(side = tkinter.TOP, fill = tkinter.BOTH, expand = 1)
            tkinter.mainloop()

    def file_open(self):
        """ menu item: open rsml file """
        global fname
        fname = tkinter.filedialog.askopenfilename(title = 'Please select a RSML root system',
                                                  filetypes = [('Image Files', ['.rsml', '.RSML', '.xml'])])
        if isinstance(fname, str):
            if fname:
                self.data.open_rsml(fname)
                self.update_all()

    def file_save(self):
        """ menu item: save rsml file (polylines)"""
        if self.data.exists():
            fname = tkinter.filedialog.asksaveasfilename(defaultextension = ".rsml")
            print(fname, type(fname))
            if isinstance(fname, str):
                if fname:
                    pd = vp.segs_to_polydata(self.data.analyser, zoom_factor = 1., param_names = ["subType", "radius", "creationTime"])
                    vt.write_rsml(fname, pd, 0, None, self.data.base_nodes)

    def file_save_vtp(self):
        """ menu item: save save vtp file containing segments (not polylines) """
        if self.data.exists():
            fname = tkinter.filedialog.asksaveasfilename(defaultextension = ".vtp")
            print(fname, type(fname))
            if isinstance(fname, str):
                if fname:
                    self.data.analyser.write(fname)  # segment analyser's writer

    def file_quit(self):
        """ menu item: quits application """
        self.root.quit()  # stops mainloop
        self.root.destroy()  # this is necessary on Windows to prevent

    def edit_add_shoot(self):
        """ adds an artifical shoot """
        if self.data.exists():
            bni = self.data.base_nodes
            if len(bni) > 1:
                self.data.add_artificial_shoot()
                self.update_all()
            else:
                tkinter.messagebox.showwarning("Warning", "Only a single base node (no roots to connect)")
        else:
            tkinter.messagebox.showwarning("Warning", "Open RSML file first")

    def edit_add_creation_times(self):
        """ adds creation times by interpolation """
        if self.data.exists():
            if not self.data.tagnames[1]:
                max_ct = tkinter.simpledialog.askstring("Creation time interpolation", "Final root system age \nin days?")
                # try:
                self.data.max_ct = float(max_ct)
                self.data.add_creation_times()
                self.update_all()
                # except:
                #     tkinter.messagebox.showerror("Error", "Something went wrong")
            else:
                tkinter.messagebox.showwarning("Warning", "Creation times are aready set by tag " + self.data.tagnames[1])
        else:
            tkinter.messagebox.showwarning("Warning", "Open RSML file first")
        pass

    def view_vtk_plot_subtype(self):
        """ menu item """
        self.view_vtk_plot("subType")

    def view_vtk_plot_length(self):
        """ menu item """
        self.view_vtk_plot("length")

    def view_vtk_plot_creationtime(self):
        """ menu item """
        self.view_vtk_plot("creationTime")

    def view_vtk_plot_multiple_suf(self):
        """ menu item """
        r = self.data.xylem_flux
        if len(self.data.base_segs) == len(self.data.base_nodes):
            n = len(self.data.base_segs)
            krs = [None] * n
            depth = [None] * n
            for i in range(0, n):
                r.dirichlet_ind = [ self.data.base_nodes[i] ]  # krs is based on dirichlet boundary condition
                krs[i], _ = r.get_krs(self.data.max_ct, [self.data.base_segs[i]])
                print("plot for base segment", self.data.base_segs[i], "krs", krs)  #
                r.neumann_ind = [self.data.base_nodes[i]]  # suf is based on neumann boundary condititon
                suf = r.get_suf(self.data.max_ct)
                self.data.analyser.addData("SUF", suf)
                self.view_vtk_plot("SUF")
        else:
            tkinter.messagebox.showwarning("Warning", "Only implemented for multiple dicots emerging from different nodes")

    def view_vtk_plot_suf(self):
        """ menu item """
        self.view_vtk_plot("SUF")

    def view_vtk_anim(self):
        """ shows animation (todo) """
        print("na")

    def view_about(self):
        """ menu item: view about dialog """
        tkinter.messagebox.showinfo("About", "RSML Viewer \nby Daniel Leitner, 2021 \n\nPart of CPlantBox")


if __name__ == '__main__':
    root = tkinter.Tk()
    app = App(root)
    root.mainloop()
