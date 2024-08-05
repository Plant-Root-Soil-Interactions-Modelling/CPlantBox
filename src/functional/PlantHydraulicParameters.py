import sys; sys.path.append(".."); sys.path.append("../..");

import numpy as np
import matplotlib.pyplot as plt

import plantbox as pb
from plantbox import PlantHydraulicParameters as PlantHydraulicParametersCPP
import rsml.rsml_reader as rsml


class PlantHydraulicParameters(PlantHydraulicParametersCPP):
    """       
        This class handles axial (kx [cm3 day-1]) and radial (kr [1 day-1]) root hydraulic conductivities  
        
        setKr            sets constant values per subType and organType
        setKrTable       sets age depent values (linearly interpolated) per subType and organType 
        setKxValues      sets values per segment
        (same for Kx)
        
        kr_f(int segment_index, double age, int subType, int organType) 
        kx_f(int segment_index, double age, int subType, int organType)
        
        see also PlantHydraulicParameters.h    
        
        @author Daniel Leitner, 2023        
    """

    @staticmethod
    def sinusoidal(t):
        """ sinusoidal function (used for transpiration) (integral over one day is 1)"""
        return np.sin(2. * np.pi * np.array(t) - 0.5 * np.pi) + 1.

    @staticmethod
    def sinusoidal2(t, dt):
        """ sinusoidal function from 6:00 - 18:00, 0 otherwise (integral over one day is 1)"""
        return np.maximum(0., np.pi * (np.cos(2 * np.pi * (t - 0.5)) + np.cos(2 * np.pi * ((t + dt) - 0.5))) / 2)

    def kr_f(self, age, st, ot = 2 , seg_ind = 0, cells = False):
        """ root radial conductivity [1 day-1] for backwards compatibility """
        return self.kr_f_cpp(seg_ind, age, st, ot)  # kr_f_cpp is XylemFlux::kr_f

    def kx_f(self, age, st, ot = 2, seg_ind = 0):
        """ root axial conductivity [cm3 day-1]  for backwards compatibility """
        return self.kx_f_cpp(seg_ind, age, st, ot)  # kx_f_cpp is XylemFlux::kx_f

    def plot_conductivities(self, monocot = True, plot_now = True, axes_ind = [], lateral_ind = []):
        """ plots conductivity  
        @param monocot      indicates if monocot (True) or dicot (False)
        @param plot_now     indicates if the figure is shown, or just retruned
        @param axes_ind     for monocots a list of three root types for "tap root", "basal", "shoot borne" [1, 4, 5];
                            for dicots a list of one root type representing the "tap root" [1]
        @param lateral_ind  for monocots a list of two root types for "1st order laterals", "2nd order laterals" [2, 3]
                            for dicots a list of three root types for "1st order laterals", "2nd order laterals", "3rd order laterals" [2, 3, 4]          
        """
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
            kx_ = [ self.kx_f(0, axes_age[i], int(st), 2) for i in range(0, len(axes_age)) ]
            ax1.plot(axes_age, kx_, axes_cols[j])
        kx_max = np.max(kx_)
        ax1.legend(axes_str)
        ax1.set_title("Axis")
        ax1.set_xlabel("age [day]")
        ax1.set_ylabel("axial conductance [cm$^3$ day$^{-1}$]")
        for j, st in enumerate(lateral_ind):
            kx_ = [ self.kx_f(0, lateral_age[i], st, 2) for i in range(0, len(lateral_age)) ]
            ax2.plot(lateral_age, kx_, axes_cols[j])
        kx_max = max(kx_max, np.max(kx_))
        ax2.legend(lateral_str)
        ax2.set_title("Laterals")
        ax2.set_xlabel("age [day]")
        ax2.set_ylabel("axial conductance [cm$^3$ day$^{-1}$]")
        for j, st in enumerate(axes_ind):
            kr_ = [ self.kr_f(0, axes_age[i], st, 2) for i in range(0, len(axes_age)) ]
            ax3.plot(axes_age, kr_, axes_cols[j])
        kr_max = np.max(kr_)
        ax3.legend(axes_str)
        ax3.set_title("Axis")
        ax3.set_xlabel("age [day]")
        ax3.set_ylabel("radial conductance [day$^{-1}$]")
        for j, st in enumerate(lateral_ind):
            kr_ = [ self.kr_f(0, lateral_age[i], st, 2) for i in range(0, len(lateral_age)) ]
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
            print("SubType {:g} for negative age: kx = {:g}, kr = {:g}".format(st, self.kx_f(0, -1, st, 2), self.kr_f(0, -1, st, 2)))
            print("SubType {:g} for old root age: kx = {:g}, kr = {:g}".format(st, self.kx_f(0, 100, st, 2), self.kr_f(0, 100, st, 2)))
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
    params.setKrTables([values], [ages])  # values, ages per organType per subType
    kx00 = np.array([[0., 1.e3]])  # artificial shoot
    kx0 = np.array([[0., 0.000864], [5., 0.00173], [12., 0.0295], [15., 0.0295], [20., 0.432], [300., 0.432]])  #  time, value; [cm3 day-1]
    kx1 = np.array([[0., 0.0000864], [5., 0.0000864], [10., 0.0000864], [12., 0.0006048], [20., 0.0006048], [23., 0.00173], [300., 0.00173]])  #  time, value; [cm3 day-1]
    ages = [kx00[:, 0], kx0[:, 0], kx1[:, 0], kx1[:, 0], kx0[:, 0], kx0[:, 0]]
    values = [kx00[:, 1], kx0[:, 1], kx1[:, 1], kx1[:, 1], kx0[:, 1], kx0[:, 1]]
    params.setKxTables([values], [ages])
    params.plot_conductivities(True)

