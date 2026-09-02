"""Carrot"""

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
from plantbox.visualisation.vtk_animate import AnimateRoots

plant = pb.Plant()

# Open plant and root parameter from a file
path = "../../modelparameter/structural/rootsystem/"
name = "Daucus_carota"
plant.readParameters(path + name + ".xml")

plant.initialize()

sim_time = 40
dt = 1  # simulation step [day], animation shows one frame per step

# Animate growth in an interactive vtk window (close window, or press 'e', to continue)
anim = AnimateRoots(plant)
anim.root_name = "age"
anim.start(axis="v")
for i in range(int(sim_time / dt)):
    plant.simulate(dt)
    anim.simtime = (i + 1) * dt
    anim.update()
anim.run()

# # Export final result (as vtp)
# plant.write("results/example_plant.vtp")

# ana = pb.SegmentAnalyser(plant)
# ana.write("results/example_plant_segs.vtp")

# # Interactive plot, using vtk
# vp.plot_plant(plant, "age")
