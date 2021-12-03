import sys; sys.path.append("../../src/python_modules/"); sys.path.append("../../")

import vtk_plot as vp
from estimate_data import EstimateDataModel
import estimate_plots

import tkinter
from tkinter import ttk
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backend_bases import key_press_handler  # Implement the default Matplotlib key bindings
import matplotlib.pyplot as plt
import numpy as np


class App:

    def __init__(self, root):
        self.data = EstimateDataModel()
        self.root = root
        self.root.wm_title("Estimate")
        self.root.geometry("850x800")
        # Menu
        menu = tkinter.Menu(root)
        menu_file = tkinter.Menu(menu, tearoff = 0)
        menu_file.add_command(label = "Open folder...", command = self.file_open)
        menu_file.add_command(label = "Save parameters...", command = self.file_save)
        menu_file.add_separator()
        menu_file.add_command(label = "Exit", command = self.file_quit)
        menu.add_cascade(label = "File", menu = menu_file)
        menu_view = tkinter.Menu(menu, tearoff = 0)
        self.root.config(menu = menu)
        # Tabs
        tabControl = ttk.Notebook(self.root)
        tab_info = ttk.Frame(tabControl)
        tab_base = ttk.Frame(tabControl)
        tab_laterals = ttk.Frame(tabControl)
        tab_base_params = ttk.Frame(tabControl)
        tab_lateral_params = ttk.Frame(tabControl)
        tabControl.add(tab_info, text = 'Information')
        tabControl.add(tab_base, text = 'Base roots')
        tabControl.add(tab_base_params, text = 'Parameters')
        tabControl.add(tab_laterals, text = 'Laterals')
        tabControl.add(tab_lateral_params, text = 'Parameters')
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
        self.combo1 = ttk.Combobox(tab_base, values = [ "Length", "Surface", "Volume"])
        self.combo1.pack(pady = 10)
        self.combo1.current(0)
        self.combo1.bind("<<ComboboxSelected>>", self.update_baseroots)
        fig, self.ax = plt.subplots(1, 1, figsize = (7, 7))
        self.canvas = FigureCanvasTkAgg(fig, master = tab_base)  # A tk.DrawingArea.
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill = tkinter.BOTH, expand = 1)
        # tab_development
        # self.combo2 = ttk.Combobox(tab_development, values = [ "Length", "Surface", "Volume"])
        # self.combo2.pack(pady = 10)
        # self.combo2.current(0)
        # self.combo2.bind("<<ComboboxSelected>>", self.update_development)
        # fig2, self.ax2 = plt.subplots(1, 1, figsize = (15, 10))
        # self.canvas2 = FigureCanvasTkAgg(fig2, master = tab_development)  # A tk.DrawingArea.
        # self.canvas2.draw()
        # self.canvas2.get_tk_widget().pack(side = tkinter.TOP, fill = tkinter.BOTH, expand = 1)

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
        # tagnames = self.data.tagnames
        # rstr = "\n"
        # lstr = "\nRadius \nCreation time \nTypes \n"
        # if tagnames[0]:
        #     rstr += "from tag '{:s}' within [{:g}, {:g}] cm\n".format(tagnames[0], np.min(self.data.radii), np.max(self.data.radii))
        # else:
        #     rstr += "not found\n"
        # if tagnames[1]:
        #     rstr += "from tag '{:s}' within [{:g}, {:g}] days\n".format(tagnames[1], np.min(self.data.cts), np.max(self.data.cts))
        # else:
        #     rstr += "not found\n"
        # if tagnames[2]:
        #     rstr += "from tag '{:s}' within [{:g}, {:g}]\n".format(tagnames[2], np.min(self.data.types), np.max(self.data.types))
        # else:
        #     rstr += "not found\n"  # , derived from root order
        # self.label_use_l.set(lstr)
        # self.label_use_r.set(rstr)

    def update_baseroots(self, event):
        """ updates base roots scatter plot """
        estimate_plots.plot_baseroots(self.data, self.ax)
        self.canvas.draw()

    def update_all(self):
        """ updates the view """
        if self.data.exists():
            self.update_info()
            self.update_baseroots(None)

    def file_open(self):
        """ menu item: open rsml file """
        fname = tkinter.filedialog.askdirectory(mustexist = True)
        if isinstance(fname, str):
            if fname:
                self.data.open_folder(fname)
                self.update_all()

    def file_save(self):
        """ menu item: save cplantbox parameters """
        pass
        # fname = tkinter.filedialog.askopenfilename(title = 'Please select a RSML root system',
        #                                           filetypes = [('Image Files', ['.rsml', '.RSML', '.xml'])])
        # if isinstance(fname, str):
        #     if fname:
        #         self.data.open_rsml(fname)
        #         self.update_all()

    def file_quit(self):
        """ menu item: quits application """
        self.root.quit()  # stops mainloop
        self.root.destroy()  # this is necessary on Windows to prevent

    def view_about(self):
        """ menu item: view about dialog """
        tkinter.messagebox.showinfo("About", "RSML Viewer \nby Daniel Leitner, 2021 \n\nPart of CPlantBox")


if __name__ == '__main__':
    root = tkinter.Tk()
    app = App(root)
    root.mainloop()
