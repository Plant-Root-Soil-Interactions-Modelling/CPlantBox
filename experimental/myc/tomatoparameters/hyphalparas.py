import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
import numpy as np

mycp = pb.MycorrhizalPlant()

path = "../tomatoparameters/"
name = "TomatoJohanna_WildTypeTwoHyphaeTypes"

# mycp.readParameters(path + name + ".xml", fromFile = True, verbose = True)

# hyphae = mycp.getOrganRandomParameter(pb.hyphae)

# v_typical_hour = 200 # micrometer per hour
# v_typical_day = v_typical_hour*24/10000 # convert to cm per day

# branching_freq_hour = 0.04 #branching events per hour
# branching_freq_day = 24*branching_freq_hour # branching events per day

# ln_typical_day = v_typical_day*branching_freq_day # typical interlaterateral distance in cm branching frequency times speed


# print("Typical speed: ",v_typical_day)
# print("Typical branching frequency: ", branching_freq_day)
# print("Typical lateral distance: ", ln_typical_day)

# for hp in hyphae:
#     hp.a = 0.01
#     hp.ln = ln_typical_day
#     hp.b_prob = 0.
#     hp.b = branching_freq_day
#     hp.lb = hp.ln
#     hp.dx = 0.01
#     hp.v = v_typical_day
#     hp.tropismS = 1.0
#     hp.distTH = 0.01   # distance for anastomosis
#     mycp.setOrganRandomParameter(hp)


new_name = "TwoHyphaePlusBAS"
# mycp.writeParameters(new_name + ".xml")


mycp.readParameters(new_name + ".xml", fromFile = True, verbose = True)

hyphae = mycp.getOrganRandomParameter(pb.hyphae)

for hp in hyphae:
    hp.tropismS = 1.0
    hp.distTH = 0.01   # distance for anastomosis
    mycp.setOrganRandomParameter(hp)


new_name = "TwoHyphaePlusBAS"
mycp.writeParameters(new_name + ".xml")