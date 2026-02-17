import sys; sys.path.append(".."); sys.path.append("../..");

import numpy as np
import matplotlib.pyplot as plt
import json

import plantbox as pb
from plantbox import PlantHydraulicParameters as PlantHydraulicParametersCPP
import rsml.rsml_reader as rsml


class PlantHydraulicParameters(PlantHydraulicParametersCPP):
    """       
        This class handles axial (kx [cm3 day-1]) and radial (kr [1 day-1]) root hydraulic conductivities, using the callback functions                  
        - kr_f(int segment_index, double age, int subType, int organType) 
        - kx_f(int segment_index, double age, int subType, int organType)        
        see also PlantHydraulicParameters.h    
        
        Values are set in dependence on subType and organType, are given
        constant:                 setKrConst, setKxConst,
        age dependent:            setKrAgeDependent, setKxAgeDependent, 
        distance dependent:       setKrDistanceDependent, setKxDistanceDependent, or  
        per segment:              setKrValues
        
        Python setter are more flexible: 
            default: subType = -1 means all subTypes, 
            or single type (int), e.g. subType = 1 
            or multiple types (list), e.g. subType = [2,3]
        
        set_kr_const, set_kx_const -> setKrConst, setKxConst
        set_kr_age_dependent, set_kx_age_dependent -> setKrAgeDependent, setKxAgeDependent
        set_kr_distance_dependent, set_kx_distance_dependent -> setKrDistanceDependent, setKxDistanceDependent
        
        it is not possible to use different methods (const, age, or distance) for different subTypes or organTypes
        it is possible to use different methods for kr and kx 
        
        can write and read parameter json files 
        
        @author Daniel Leitner, 2025     
    """

    @staticmethod
    def sinusoidal(t):
        """ sinusoidal function (used for transpiration) (integral over one day is 1)"""
        return np.sin(2. * np.pi * np.array(t) - 0.5 * np.pi) + 1.

    @staticmethod
    def sinusoidal2(t, dt):
        """ sinusoidal function from 6:00 - 18:00, 0 otherwise (integral over one day is 1)"""
        return np.maximum(0., np.pi * (np.cos(2 * np.pi * (t - 0.5)) + np.cos(2 * np.pi * ((t + dt) - 0.5))) / 2)

    def set_flexible_(self, f, params, subType, organType):
        """ calls cpp functions @param f with parameters @param params,  
        enabling multiple subTypes and multiple organTypes """
        if subType == -1:
            subTypes = range(0, 10)
        else:
            if isinstance(subType, int):
                subTypes = [subType]
            else:
                subTypes = subType
        if organType == -1:
            organTypes = [pb.OrganTypes.root, pb.OrganTypes.stem, pb.OrganTypes.leaf]
        else:
            if isinstance(organType, (int, pb.OrganTypes)):
                organTypes = [organType]
            else:
                organTypes = organType
        for ot in organTypes:
            for st in subTypes:
                # print(f, st, ot)
                f(*params, st, ot)

    def set_kr_const(self, kr, subType = -1, organType = int(pb.OrganTypes.root)):
        """ sets a constant radial conductivity [1 day-1] for roots for a subTypes (default: all 0-9) and organTypes (default: roots)"""
        self.set_flexible_(self.setKrConst, (kr,), subType, organType)

    def set_kx_const(self, kx, subType = -1, organType = int(pb.OrganTypes.root)):
        """ sets a constant axial conductivity [cm3/day] for roots for a subType (-1 for all subTypes) """
        self.set_flexible_(self.setKxConst, (kx,), subType, organType)

    def set_kr_age_dependent(self, age, values, subType = -1, organType = int(pb.OrganTypes.root)):
        """ sets an age dependent radial conductivity [1 day-1] for roots for a subTypes (default: all 0-9) and organTypes (default: roots)"""
        self.set_flexible_(self.setKrAgeDependent, (age, values,), subType, organType)

    def set_kx_age_dependent(self, age, values, subType = -1, organType = int(pb.OrganTypes.root)):
        """ sets an age dependent axial conductivity [cm3/day] for roots for a subType (-1 for all subTypes) """
        self.set_flexible_(self.setKxAgeDependent, (age, values), subType, organType)

    def set_kr_distance_dependent(self, age, values, subType = -1, organType = int(pb.OrganTypes.root)):
        """ sets a distance dependent radial conductivity [1 day-1] for roots for a subTypes (default: all 0-9) and organTypes (default: roots)"""
        self.set_flexible_(self.setKrDistanceDependent, (age, values,), subType, organType)

    def set_kx_distance_dependent(self, age, values, subType = -1, organType = int(pb.OrganTypes.root)):
        """ sets a distant dependent axial conductivity [cm3/day] for roots for a subType (-1 for all subTypes) """
        self.set_flexible_(self.setKxDistanceDependent, (age, values), subType, organType)

    def write_parameters(self, filename):
        """ writes the parameters into a json file """
        if self.krMode == "perSegment" or self.kxMode == "perSegment":
            print("PlantHydraulicParameters.write_parameters(): Warning! perSegment works only if you reload the exact same root geometry (with the exact same segments)")
        mode = {"krMode": self.krMode, "kxMode": self.kxMode }
        json_dict = {"mode": mode, "kx_ages": self.kx_ages, "kx_values": self.kx_values, "kr_ages": self.kr_ages, "kr_values": self.kr_values,
                     "krPerSegment": self.krValues, "kxPerSegment": self.kxValues }
        with open(filename + ".json", "w+") as f:
            json.dump(json_dict, f)

    def read_parameters(self, filename):
        """ writes the parameters into a json file """
        with open(filename + ".json", "r") as f:
            json_dict = json.load(f)
        self.setKrValues(json_dict["krPerSegment"])  # can be empty
        self.setKxValues(json_dict["kxPerSegment"])
        for ot in [int(pb.OrganTypes.root), int(pb.OrganTypes.stem), int(pb.OrganTypes.leaf)]:
            for st in range(0, self.maxSubTypes):
                age_ = json_dict["kx_ages"][str(ot)][st]
                value_ = json_dict["kx_values"][str(ot)][st]
                self.setKxAgeDependent(age_, value_, st, ot)
                age_ = json_dict["kr_ages"][str(ot)][st]
                value_ = json_dict["kr_values"][str(ot)][st]
                self.setKrAgeDependent(age_, value_, st, ot)
        mode = json_dict["mode"]
        self.setMode(mode["krMode"], mode["kxMode"])

    def plot_conductivities(self, monocot = True, plot_now = True, axes_ind = [], lateral_ind = []):
        """ plots conductivity  
        @param monocot      indicates if monocot (True) or dicot (False)
        @param plot_now     indicates if the figure is shown, or just retruned
        @param axes_ind     for monocots a list of three root types for "tap root", "basal", "shoot borne" [1, 4, 5];
                            for dicots a list of one root type representing the "tap root" [1]
        @param lateral_ind  for monocots a list of two root types for "1st order laterals", "2nd order laterals" [2, 3]
                            for dicots a list of three root types for "1st order laterals", "2nd order laterals", "3rd order laterals" [2, 3, 4]          
        """
        ot_root = pb.OrganTypes.root
        axes_age = naxes_age = np.linspace(0, 50, 500)
        lateral_age = np.linspace(0, 25, 125)
        lateral_cols = ["r", "g:", "m--", "b--"]
        axes_cols = ["r", "g:", "m--", "b--"]
        if monocot:
            axes_str = ["tap root", "basal", "shoot borne"]
            lateral_str = ["1st order laterals", "2nd order laterals"]
            if axes_ind == []:
                axes_ind = [1, 4, 5]
            if lateral_ind == []:
                lateral_ind = [2, 3]
        else:  # dicot
            axes_str = ["tap root"]
            lateral_str = ["1st order laterals", "2nd order laterals", "3rd order laterals"]
            if axes_ind == []:
                axes_ind = [1]
            if lateral_ind == []:
                lateral_ind = [2, 3, 4]

        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize = (16, 10))
        for j, st in enumerate(axes_ind):
            kx_ = [ self.kx_f(0, axes_age[i], st, ot_root) for i in range(0, len(axes_age)) ]
            ax1.plot(axes_age, kx_, axes_cols[j])
        kx_max = np.max(kx_)
        ax1.legend(axes_str)
        ax1.set_title("Axis")
        ax1.set_xlabel("age [day]")
        ax1.set_ylabel("axial conductance [cm$^3$ day$^{-1}$]")
        for j, st in enumerate(lateral_ind):
            kx_ = [ self.kx_f(0, lateral_age[i], st, ot_root) for i in range(0, len(lateral_age)) ]
            ax2.plot(lateral_age, kx_, axes_cols[j])
        kx_max = max(kx_max, np.max(kx_))
        ax2.legend(lateral_str)
        ax2.set_title("Laterals")
        ax2.set_xlabel("age [day]")
        ax2.set_ylabel("axial conductance [cm$^3$ day$^{-1}$]")
        for j, st in enumerate(axes_ind):
            kr_ = [ self.kr_f(0, axes_age[i], st, ot_root) for i in range(0, len(axes_age)) ]
            ax3.plot(axes_age, kr_, axes_cols[j])
        kr_max = np.max(kr_)
        ax3.legend(axes_str)
        ax3.set_title("Axis")
        ax3.set_xlabel("age [day]")
        ax3.set_ylabel("radial conductance [day$^{-1}$]")
        for j, st in enumerate(lateral_ind):
            kr_ = [ self.kr_f(0, lateral_age[i], st, ot_root) for i in range(0, len(lateral_age)) ]
            ax4.plot(lateral_age, kr_, axes_cols[j])
        kr_max = max(kr_max, np.max(kr_))
        ax4.legend(lateral_str)
        ax4.set_title("Laterals")
        ax4.set_xlabel("age [day]")
        ax4.set_ylabel("radial conductance [day$^{-1}$]")
        print(kx_max)
        print(kr_max)
        # ax1.set_ylim([0, kx_max * 1.1])
        # # ax2.set_ylim([0, kx_max * 1.1])
        # ax3.set_ylim([0, kr_max * 1.1])
        # ax4.set_ylim([0, kr_max * 1.1])
        print()
        for st in range(0, 5):
            print("SubType {:g} for negative age: kx = {:g}, kr = {:g}".format(st, self.kx_f(0, -1, st, ot_root), self.kr_f(0, -1, st, ot_root)))
            print("SubType {:g} for old root age: kx = {:g}, kr = {:g}".format(st, self.kx_f(0, 100, st, ot_root), self.kr_f(0, 100, st, ot_root)))
        print()
        if plot_now:
            plt.tight_layout()
            plt.show()
        return fig


if __name__ == "__main__":

    params = PlantHydraulicParameters()

    # Maize, from Couvreur et al. (2012), originally from Doussan et al. (1998)
    kr00 = np.array([[0., 0.]])  # artificial shoot
    kr0 = np.array([[-1e4, 0.], [-0.1, 0.], [0., 0.000181], [8., 0.000181], [10, 0.0000648], [18, 0.0000648], [25, 0.0000173], [300, 0.0000173]])  # time, value; [day], [1 day-1]
    kr1 = np.array([[-1e4, 0.], [-0.1, 0.], [0., 0.000181], [10., 0.000181], [16, 0.0000173], [300, 0.0000173]])  #   time, value; [1 day-1]
    values = [kr00[:, 1], kr0[:, 1], kr1[:, 1], kr1[:, 1], kr0[:, 1], kr0[:, 1]]
    ages = [kr00[:, 0], kr0[:, 0], kr1[:, 0], kr1[:, 0], kr0[:, 0], kr0[:, 0]]
    [params.setKrAgeDependent(ages[i], values[i], i) for i in range(0, 6)]
    kx00 = np.array([[0., 1.e3]])  # artificial shoot
    kx0 = np.array([[0., 0.000864], [5., 0.00173], [12., 0.0295], [15., 0.0295], [20., 0.432], [300., 0.432]])  #  time, value; [cm3 day-1]
    kx1 = np.array([[0., 0.0000864], [5., 0.0000864], [10., 0.0000864], [12., 0.0006048], [20., 0.0006048], [23., 0.00173], [300., 0.00173]])  #  time, value; [cm3 day-1]
    ages = [kx00[:, 0], kx0[:, 0], kx1[:, 0], kx1[:, 0], kx0[:, 0], kx0[:, 0]]
    values = [kx00[:, 1], kx0[:, 1], kx1[:, 1], kx1[:, 1], kx0[:, 1], kx0[:, 1]]
    [params.setKxAgeDependent(ages[i], values[i], i) for i in range(0, 6)]
    params.plot_conductivities(True)

    params.write_parameters("couvreur2012")

    params2 = PlantHydraulicParameters()
    params2.read_parameters("couvreur2012")
    params2.plot_conductivities(True)
