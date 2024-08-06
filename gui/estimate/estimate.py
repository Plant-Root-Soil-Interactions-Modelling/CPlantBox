import sys; sys.path.append("../.."); sys.path.append("../../src/")

import visualisation.vtk_plot as vp
from estimate_data import EstimateDataModel
import estimate_plots

import tkinter
from tkinter import simpledialog
from tkinter import ttk
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backend_bases import key_press_handler  # Implement the default Matplotlib key bindings
import matplotlib.pyplot as plt
import numpy as np

SMALL_SIZE = 16
MEDIUM_SIZE = 16
BIGGER_SIZE = 16
plt.rc('font', size = SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize = SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize = MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize = SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize = BIGGER_SIZE)  # fontsize of the figure title


class App:

    def __init__(self, root):
        self.data = None
        self.root = root
        self.root.wm_title("Estimate")
        self.root.geometry("1200x800")
        # Menu
        menu = tkinter.Menu(self.root)
        menu_file = tkinter.Menu(menu, tearoff = 0)
        menu_file.add_command(label = "Open folder...", command = self.file_open)
        menu_file.add_command(label = "Save parameters...", command = self.file_save)
        menu_file.add_separator()
        menu_file.add_command(label = "Exit", command = self.file_quit)
        menu.add_cascade(label = "File", menu = menu_file)
        menu_edit = tkinter.Menu(menu, tearoff = 0)
        menu_edit.add_command(label = "Shift time...", command = self.edit_shift_time)
        menu.add_cascade(label = "Edit", menu = menu_edit)
        menu_view = tkinter.Menu(menu, tearoff = 0)
        menu_view.add_command(label = "Estimated root age...", command = self.view_vtk_plot_age)
        menu_view.add_separator()
        menu_view.add_command(label = "About...", command = self.view_about)
        menu.add_cascade(label = "View", menu = menu_view)
        self.root.config(menu = menu)
        # top frame
        self.tab_base_top_frame = ttk.Frame(self.root)
        self.tab_base_top_frame.pack(side = tkinter.TOP)
        self.combo0 = ttk.Combobox(self.tab_base_top_frame, values = [ "Apical delay", "Apical length"], width = 15)  #
        self.combo0.pack(pady = 10, side = tkinter.LEFT, expand = 1)
        self.combo0.current(0)
        self.combo0.bind("<<ComboboxSelected>>", self.update_all)
        self.combo1 = ttk.Combobox(self.tab_base_top_frame, values = [ "Multiple dicots", "Multiple dicots (fixed lmax)",
                                                       "Monocot linear model", "Monocot linear model (fixed lmax)"], width = 30)
        self.combo1.pack(pady = 10, side = tkinter.LEFT, expand = 1)
        self.combo1.current(0)
        self.combo1.bind("<<ComboboxSelected>>", self.update_all)
        self.combo2 = ttk.Combobox(self.tab_base_top_frame, values = [ "per order", "per RSML type tag", "per clustering 2", "per clustering 3"], width = 15)
        self.combo2.pack(pady = 10, side = tkinter.LEFT, expand = 1)
        self.combo2.current(0)
        self.combo2.bind("<<ComboboxSelected>>", self.update_all)
        self.label1 = ttk.Label(self.tab_base_top_frame, text = "lmax 0 order roots")
        self.label1.pack(padx = 10, expand = 1, side = tkinter.LEFT)
        self.lmax0 = tkinter.StringVar()
        self.lmax0.set("100")
        self.entry1 = ttk.Entry(self.tab_base_top_frame, textvariable = self.lmax0, width = 7)
        self.entry1.pack(pady = 10, expand = 1, side = tkinter.LEFT)
        self.lmax1 = tkinter.StringVar()
        self.lmax1.set("50")
        self.label2 = ttk.Label(self.tab_base_top_frame, text = "1st order")
        self.label2.pack(padx = 10, expand = 1, side = tkinter.LEFT)
        self.entry2 = ttk.Entry(self.tab_base_top_frame, textvariable = self.lmax1, width = 7)
        self.entry2.pack(pady = 10, expand = 1, side = tkinter.LEFT)
        self.lmax2 = tkinter.StringVar()
        self.lmax2.set("10")
        self.label3 = ttk.Label(self.tab_base_top_frame, text = "2nd order")
        self.label3.pack(padx = 10, expand = 1, side = tkinter.LEFT)
        self.entry3 = ttk.Entry(self.tab_base_top_frame, textvariable = self.lmax2, width = 7)
        self.entry3.pack(pady = 10, expand = 1, side = tkinter.LEFT)
        # Tabs
        tabControl = ttk.Notebook(self.root)
        tab_info = ttk.Frame(tabControl)
        tab_base = ttk.Frame(tabControl)
        tab_laterals = ttk.Frame(tabControl)
        tab_base_params = ttk.Frame(tabControl)
        tab_plant_params = ttk.Frame(tabControl)
        tab_lateral_params = ttk.Frame(tabControl)
        tabControl.add(tab_info, text = 'Information')
        tabControl.add(tab_base, text = 'Base roots')
        tabControl.add(tab_laterals, text = 'Laterals')
        tabControl.add(tab_base_params, text = 'Parameters')
        tabControl.add(tab_plant_params, text = 'Plant Parameters')
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
        # tab_base
        fig, self.ax = plt.subplots(1, 1, figsize = (7, 7))
        self.canvas = FigureCanvasTkAgg(fig, master = tab_base)  # A tk.DrawingArea.
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill = tkinter.BOTH, expand = 1, side = tkinter.BOTTOM)
       # tab_laterals
        fig2, self.ax2 = plt.subplots(1, 1, figsize = (7, 7))
        self.canvas2 = FigureCanvasTkAgg(fig2, master = tab_laterals)  # A tk.DrawingArea.
        self.canvas2.draw()
        self.canvas2.get_tk_widget().pack(fill = tkinter.BOTH, expand = 1)
        # tab_base_parameters
        self.label_basal_params_l, self.label_basal_params_r = [], []
        names_ = ['Tap and basal root parameters [mean, sd]', 'First order [mean, sd]', 'Second order [mean, sd]']  # TODO clustering
        for i in range(0, 3):
            lf_basal_params = ttk.LabelFrame(tab_base_params, text = names_[i])
            lf_basal_params.grid(column = 0, row = i, padx = 20, pady = 10)
            self.label_basal_params_l.append(tkinter.StringVar())
            self.label_basal_params_r.append(tkinter.StringVar())
            ttk.Label(lf_basal_params, textvariable = self.label_basal_params_l[-1], anchor = "w", width = 50).grid(column = 0, row = 0)
            ttk.Label(lf_basal_params, textvariable = self.label_basal_params_r[-1], anchor = "w", width = 50).grid(column = 1, row = 0)
        # tab_plant_parameters
        self.label_plant_params_l, self.label_plant_params_r = [], []
        names_ = 'Plant parameters [mean, sd]'
        lf_plant_params = ttk.LabelFrame(tab_plant_params, text = names_)
        lf_plant_params.grid(column = 0, row = 1, padx = 20, pady = 10)
        self.label_plant_params_l = tkinter.StringVar()
        self.label_plant_params_r = tkinter.StringVar()
        ttk.Label(lf_plant_params, textvariable = self.label_plant_params_l, anchor = "w", width = 50).grid(column = 0, row = 0)
        ttk.Label(lf_plant_params, textvariable = self.label_plant_params_r, anchor = "w", width = 70).grid(column = 1, row = 0)

    def parse_gui(self):
        """ converts values of Entry fields into CPlantbox lmax parameter """
        self.data.parameters[0].lmax = float(self.lmax0.get())
        self.data.parameters[1].lmax = float(self.lmax1.get())
        self.data.parameters[2].lmax = float(self.lmax2.get())

    def update_info(self):
        """ update info tab """
        lstr = "\nFolder\nSoftware\nNumber of files\nNumber of base roots\nNumber of laterals\n"
        metadata = self.data.rsmls[0].metadata
        fname = self.data.folder_name
        nof = len(self.data.rsmls)
        brc = 0
        for root in self.data.base_root_indices:
            for _ in root:
                brc += 1
        lc = 0
        for root in self.data.lat_root_indices:
            for _ in root:
                lc += 1
        rstr = "\n{:s}\n{:s}\n{:g}\n{:g}\n{:g}\n".format(fname, metadata.software, nof, brc, lc)
        self.label_general_l.set(lstr)
        self.label_general_r.set(rstr)
        # label_prop
        lstr, rstr = "\n", "\n"
        for k in self.data.rsmls[0].properties.keys():
            v = np.array(self.data.rsmls[0].properties[k])
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
        for k in self.data.rsmls[0].functions.keys():
            v_ = self.data.rsmls[0].functions[k]
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
        tagnames = self.data.rsmls[0].tagnames
        rstr = "\n"
        lstr = "\nRadius \nCreation time \nTypes \n"
        if tagnames[0]:
            rstr += "from tag '{:s}' within [{:g}, {:g}] cm\n".format(tagnames[0], np.min(self.data.rsmls[0].radii), np.max(self.data.rsmls[0].radii))
        else:
            rstr += "not found\n"
        if tagnames[1]:
            rstr += "from tag '{:s}' within [{:g}, {:g}] days\n".format(tagnames[1], np.min(self.data.rsmls[0].cts), np.max(self.data.rsmls[0].cts))
        else:
            rstr += "not found (trying to parse filenames)\n"
        if tagnames[2]:
            rstr += "from tag '{:s}' within [{:g}, {:g}]\n".format(tagnames[2], np.min(self.data.rsmls[0].types), np.max(self.data.rsmls[0].types))
        else:
            rstr += "not found\n"  # , derived from root order
        self.label_use_l.set(lstr)
        self.label_use_r.set(rstr)

    def parameters_rstr_(self, index):
        p = self.data.parameters[index]
        rstr = "{:g} +-{:g}\n{:g} +-{:g}\n\n".format(p.r, p.rs, p.lmax, p.lmaxs)
        rstr += "{:g} +-{:g}\n{:g} +-{:g}\n{:g} +-{:g}\n\n".format(p.lb, p.lbs, p.la, p.las, p.ln, p.lns)
        rstr += "{:g} +-{:g}\n{:g} +-{:g}\n".format(p.a, p.a_s, p.theta, p.thetas)
        rstr += "\n"
        # int tropismT = 1;        ///< Root tropism parameter (Type)
        # double tropismN = 1.;   ///< Root tropism parameter (number of trials)
        # double tropismS = 0.2;  ///< Root tropism parameter (mean value of expected changeg) [1/cm]
        # double ldelay = 1.;     ///< Lateral root emergence delay [day], only used by RootDelay, @see RootDelay, RootSystem::initializeDB
        # double ldelays = 0.;     ///< Standard deviation of lateral root emergence delay [day]
        return rstr

    def pparameters_pstr_(self):
        srp = self.data.pparameters
        pstr = "\n\n{:g} +-{:g}\n{:g} +-{:g}\n{:g} +-{:g}\n\n".format(srp.delayB, srp.delayBs, srp.firstB, srp.firstBs, srp.maxB, srp.maxBs)
        pstr += "\n"
        return pstr

    def update_parameters_tap(self):
        """ update first parameter tap """
        for i in range(0, len(self.label_basal_params_l)):
        # for i in range(0, np.max(self.data.orders)+1):
            lstr = "\nInitial growth rate [cm day-1]\nMaximal root length [cm]\n\n"
            lstr += "Basal zone [cm]\nApical zone [cm]\nInter-lateral distance [cm]\n\n"
            lstr += "Root radius [cm]\nAngle between root and parent root [deg]\n"
            lstr += "Lateral subtypes (type, propability; ...)\n"
            self.label_basal_params_l[i].set(lstr)
            self.label_basal_params_r[i].set(self.parameters_rstr_(i))

    def update_plant_parameters(self):
        """ update first parameter tap """
        lstr = "\nTime delay between the basal roots [day]\nEmergence of first basal root [day]\nMaximal number of basal roots [1]\n\n"
        self.label_plant_params_l.set(lstr)
        self.label_plant_params_r.set(self.pparameters_pstr_())

    def update_all(self, event = None):
        """ updates the view """
        if self.data:
            self.parse_gui()  # fills self.data.parameters[:].lmax values
            self.update_info()
            self.data.create_params(self.combo0.current(), self.combo1.current(), self.combo2.current())  # does the fitting for the current settings
            self.update_parameters_tap()
            self.update_plant_parameters()
            estimate_plots.plot_baseroots(self.data, self.combo1.current(), self.ax)
            estimate_plots.plot_laterals(self.data, self.combo1.current(), self.combo2.current(), self.ax2)
            self.canvas.draw()
            self.canvas2.draw()

    def file_open(self):
        """ menu item: open rsml file """
        fname = tkinter.filedialog.askdirectory(mustexist = True)
        if isinstance(fname, str):
            if fname:
                self.data = EstimateDataModel()  # new data model
                self.data.open_folder(fname)
                self.update_all()

    def file_save(self):
        """ menu item: save cplantbox parameters """
        if self.data:
            fname = tkinter.filedialog.asksaveasfilename(defaultextension = ".xml")
            print(fname)
            print(type(fname))
            if isinstance(fname, str):
                if fname:
                    self.data.write_parameters(fname)
        else:
            tkinter.messagebox.showwarning("Warning", "Open a rsml folder first")

    def file_quit(self):
        """ menu item: quits application """
        self.root.quit()  # stops mainloop
        self.root.destroy()  # this is necessary on Windows to prevent

    def view_about(self):
        """ menu item: view about dialog """
        tkinter.messagebox.showinfo("About", "Estimate (parametrisation tool) \nby Daniel Leitner, 2023 \n\nPart of CPlantBox")

    def edit_shift_time(self):
        """ shifts the measurement time """
        if self.data:
            shift_str = simpledialog.askfloat("Input", "Shift measurement time for [day]", parent = self.root)
            shift = float(shift_str)
            for i in range(0, len(self.data.times)):
                self.data.times[i] += shift
            self.update_all()
        else:
            tkinter.messagebox.showwarning("Warning", "Open a rsml folder first")

    def view_vtk_plot_age(self):
        """ vtk plot with estimated root ages """
        if self.data:
            # shift_str = simpledialog.askfloat("Input", "Shift measurement time for [day]", parent = self.root)
            vp.plot_roots(self.data.analyser, name)
        else:
            tkinter.messagebox.showwarning("Warning", "Open a rsml folder first")


if __name__ == '__main__':
    root = tkinter.Tk()
    app = App(root)
    root.mainloop()
