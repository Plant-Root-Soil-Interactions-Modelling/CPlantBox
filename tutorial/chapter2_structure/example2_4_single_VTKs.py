"""increase axial resolution (e.g. for animation)"""

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp  # |\label{3f:importvtk}|

# plant  |\label{3f:plantStart}|
plant = pb.MappedPlant(0)
path = "../../modelparameter/structural/plant/"
filename = "fspm2023"
plant.readParameters(path + filename + ".xml")
sim_time = 60.0  # days

# Parameters for animation
sim_time = 15  # days
fps = 24  # frames per second
anim_time = 5  # seconds
n_steps = fps * anim_time
dt = sim_time / n_steps  # days
filename = "animate"

# Modify axial resolution of all relevant organs
for organ_type in [pb.root, pb.stem, pb.leaf]:
    for p in plant.getOrganRandomParameter(organ_type):
        p.dxMin = 0.05  # cm
        p.dx = 0.1  # adjust resolution (cm)

# Simulate
plant.initialize()
for i in range(0, n_steps):
    plant.simulate(dt)
    vp.write_plant(f"results/{filename}{i:04d}", plant)
    print(i + 1, "/", n_steps)
