import sys; sys.path.append("../../src/python_modules/"); sys.path.append("../../")

from rsml_data import RsmlData
import plantbox as pb
import rsml_reader
from estimate_params import *

import numpy as np
import os


class EstimateDataModel:
    """ 
    Model in the sense of model view controller (MVC), stores most that is presented in the view    
    """

    def __init__(self):
        self.rsmls = []  # list of rsml data to analyse, rsml_data.RsmlData
        self.estimates = [{}]  # list of dictionary of parameters per root
        self.parameters = []  # list of CPlantBox parameters of type , index = subType, max = 10
        self.times = []
        self.base_root_indices = []
        self.tap_root_indices = []
        self.basal_root_indices = []
        self.root_indices = []
        self.folder_name = ""
        self.file_names = []  # for debugging

    def open_folder(self, folder_name):
        """ see RsmlData.open_rsml() in src/python_modules/rsml_data.py                
        """
        self.rsmls = []  # delete old data
        self.folder_name = folder_name
        for root, subdirs, files in os.walk(folder_name):
            for filename in files:
                if filename.endswith('.rsml'):
                    file_path = os.path.join(root, filename)
                    print()
                    print('file %s (full path: %s)\n' % (filename, file_path))
                    print()
                    self.file_names.append(filename)
                    file_data = RsmlData()
                    file_data.open_rsml(file_path)
                    if not "order" in file_data.properties:  # check if orders were created, if not create tag
                        file_data.properties["order"] = rsml_reader.get_root_orders(file_data.properties)
                    self.rsmls.append(file_data)
                    try:
                        s = filename.split(".")
                        str = ''.join([n for n in s[0][-3:] if n.isdigit()])
                        self.times.append(int(str))
                    except:
                        print("filename", filename)
                        self.times.append(0)
        self.estimates = [{}] * len(self.times)
        plant = pb.Plant()  # dummy parameter
        self.parameters = [pb.RootRandomParameter(plant) for _ in range(0, 10)]
        self.create_length()  # add a length tag to properties
        self.initialize_roots_()  # find base roots

    def exists(self):
        """ true if a rsml folder was set """
        return len(self.rsmls) > 0

    def initialize_roots_(self):
        """ sets base root indices, tap root indices, and basal root indices.
        base root indices are the sum of taps and basals (treating shootborne and basals as equal for now) """
        # 1. find base root indices and all indices
        self.base_root_indices = []
        self.root_indices = []
        for i, data in enumerate(self.rsmls):  # plants
            pp = data.properties["parent-poly"]
            self.base_root_indices.append([])
            self.root_indices.append([])
            for j, _ in enumerate(data.polylines):
                self.root_indices[-1].append(j)
                if pp[j] == -1:  # add if base root
                    self.base_root_indices[-1].append(j)
        # 2. find tap and basal indices
        self.tap_root_indices = []
        self.basal_root_indices = []
        for i in range(0, len(self.times)):  # plants
            bri = self.base_root_indices[i]
            self.tap_root_indices.append([])  # case for not found
            self.basal_root_indices.append([])
            software = self.rsmls[i].metadata.software
            if (software == "OpenSimRoot" or software == "archisimple") and (len(bri) > 1):  # some come with an artificial shoot
                self.tap_root_indices[-1].append(bri[1])
                startindex = 2
                print("speical case for ", software)
            else:
                self.tap_root_indices[-1].append(bri[0])
                startindex = 1
            for j in bri[startindex:]:  # other base roots
                self.basal_root_indices[-1].append(j)
        # print("taps")
        # print(self.tap_root_indices)
        # print("basals")
        # print(self.basal_root_indices)

    def create_length(self):
        """ calculates root lengths, and adds it to properties if not already present)"""
        for k, _ in enumerate(self.rsmls):  # measurements
            if not "length" in self.rsmls[k].properties:
                self.rsmls[k].properties["length"] = []
                for i, p in enumerate(self.rsmls[k].polylines):
                    l = 0.
                    for j, n1 in enumerate(p[:-1]):
                        n2 = p[j + 1]
                        l += np.linalg.norm(np.array(n2) - np.array(n1))
                    self.rsmls[k].properties["length"].append(l)
        else:
            print("EstimateDataModel.create_length: 'length' tag is already available")

    def create_params(self, base_method, calibration_method):
        """
        """
        self.estimate_zones_(self.root_indices)
        p = self.parameters[0]
        if base_method < 2:  # tap and basals are treated the same way (multiple dicots)
            self.aggregate_parameters_(self.base_root_indices, target_type = 0)

        for i, j_ in enumerate(self.base_root_indices):
            for j in j_:
                self.estimates[i][(j, "age")] = self.times[i]

        order = 0
        indices = self.base_root_indices  # starting with base roots indices, then is ascending topological order
        while len(indices) > 0:
            print("computing order ", order)
            self.add_r_(indices, calibration_method, order)
            self.compute_age(indices, order)

            # TODO lateral callibration ..., to obtain new r, k for laterals

            order += 1
            indices = self.pick_order(order)
            print("new length", len(indices))  # TODO change while criteria [[],[],[]] will pass

    def estimate_zones_(self, indices):
        """ creates lb, ln, la per root (if possible)
        todo theta, r
        see aggregate_parameters, to obtain mean and sd 
        """
        for i, j_ in enumerate(indices):
            for j in j_:
                baseroot_length = self.rsmls[i].properties["length"][j]
                kids = self.pick_kids(i, j)  # print("kids", kids[:, 0])
                # print(kids.shape)
                if kids.shape[0] > 0:
                    ii = [self.rsmls[i].properties["parent-node"][kids[m, 0]] for m in range (0, kids.shape[0])]
                    base_polyline = self.rsmls[i].polylines[j]
                    ii.sort()
                    n0 = 0
                    nl0 = ii[0]
                    lb = polyline_length(n0, nl0, base_polyline)  # length of basal zone
                    ne = ii[-1]
                    nle = len(base_polyline) - 1
                    la = polyline_length(ne, nle, base_polyline)  # length of apical zone
                    ln_ = []
                    for k in range(0, len(ii) - 1):  # laterals
                        ln_.append(polyline_length(ii[k], ii[k + 1], base_polyline))  # internodal distances
                    # print("found values for ", i, j)
                    self.estimates[i][(j, "lb")] = lb
                    self.estimates[i][(j, "ln")] = ln_
                    self.estimates[i][(j, "la")] = la

    def aggregate_parameters_(self, indices, target_type):
        """ aggregates the individual root parameters into mean and sd, 
        and adds it to the parameters (of type list of RootRandomParameters) at index target_type  
        """
        la_, lb_, ln_ = [], [], []
        # a_, theta
        for i, j_ in enumerate(indices):
            for j in j_:
                # print(i, j)
                if (j, "la") in self.estimates[i]:
                    la_.append(self.estimates[i][(j, "la")])
                if (j, "lb") in self.estimates[i]:
                    lb_.append(self.estimates[i][(j, "lb")])
                if (j, "ln") in self.estimates[i]:
                    ln_.extend(self.estimates[i][(j, "ln")])
        # print(la_)
        # print(lb_)
        # print(ln_)
        p = self.parameters[target_type]
        p.la = np.mean(la_)
        p.las = np.std(la_)
        p.lb = np.mean(lb_)
        p.lbs = np.std(lb_)
        p.ln = np.mean(ln_)
        p.lns = np.std(ln_)

    def compute_age(self, indices, target_type):
        p = self.parameters[target_type]
        for i, j_ in enumerate(indices):
            for j in j_:
                t = self.estimates[i][(j, "age")]
                kids = self.pick_kids(i, j)
                if len(kids) > 0:  # add estimated age
                    for k in kids[:, 0]:  # kids
                        # print("age", i, j)
                        r = self.estimates[i][(j, "r")]
                        la = self.estimates[i][(j, "la")]
                        il = self.rsmls[i].properties["parent-node"][k]
                        bl = polyline_length(0, il, self.rsmls[i].polylines[j])  # latearal base length
                        bl_la = min(bl + la, p.lmax * 0.99)
                        et = negexp_age(max(bl_la, 0.), r, p.lmax)  # time the lateral emerged
                        age = max(t - et, 0.)
                        age = min(age, t)
                        self.estimates[i][(k, "age")] = age  # tap root et =  0

    def add_r_(self, indices, calibration_method, target_type):
        # add r as drawn from the distribution
        p = self.parameters[target_type]
        for i, j_ in enumerate(indices):
            for j in j_:
                l = self.rsmls[i].properties["length"][j]
                t = self.times[i]
                r = negexp_rate(l, p.lmax, t)
                self.estimates[i][(j, "r")] = r

    def pick_order(self, root_order, order_tag = "order"):
        """
        """
        indices = []
        for i, data in enumerate(self.rsmls):  # plants
            indices.append([])
            for j, _ in enumerate(data.polylines):
                o = data.properties[order_tag][j]
                if o == root_order:  # add if of right topological order
                    indices[-1].append(j)
        return indices

    def pick_kids(self, measurement_index, base_root_index):
        """ returns the indices of all roots that are directly connected with 
        root @param index base_root in measurement @param measurement_index
        """
        ppi = np.array(self.rsmls[measurement_index].properties["parent-poly"])
        # print(base_root_index)
        # print("ppi", ppi)
        kids_indices = np.array(np.argwhere(ppi == base_root_index * np.ones(ppi.shape)), dtype = np.int64)
        return kids_indices

