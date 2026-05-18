"""
GrassLeaf animation example
============================
Demonstrates how to animate a growing GrassLeaf plant step-by-step in an
interactive VTK window using :class:`AnimateRoots` from ``vtk_animate.py``.

The simulation runs for ``total_days`` with a time step of ``dt`` days.
After each step the render window is refreshed so growth can be watched live.
If ``save_frames = True`` the frames are also written to JPEG files that can
later be assembled into a video with, e.g.::

    ffmpeg -r 10 -i grassleaf_%d.jpg -c:v libx264 grassleaf.mp4

Run from this directory:
    python grassleaf_animate.py
"""

import plantbox as pb
from plantbox.visualisation.vtk_animate import AnimateRoots

# --------------------------------------------------------------------------- #
#  Subclass: use GrassLeaf instead of the default Leaf organ
# --------------------------------------------------------------------------- #


class Poaceae(pb.Plant):
    """Plant subclass that creates :class:`pb.GrassLeaf` organs."""

    def createLeaf(self, subType, delay, parent, pni):
        return pb.GrassLeaf(self, subType, delay, parent, pni)


# --------------------------------------------------------------------------- #
#  1.  Build the plant and set random parameters
# --------------------------------------------------------------------------- #

plant = Poaceae()

# -- Seed --
seed_rp = pb.SeedRandomParameter(plant)
seed_rp.subType = 0
plant.setOrganRandomParameter(seed_rp)

# -- Tap root (short anchor only) --
root_rp = pb.RootRandomParameter(plant)
root_rp.subType = 1
root_rp.lmax = 0.1
root_rp.r = 1.0
root_rp.theta = 0.0
plant.setOrganRandomParameter(root_rp)

# -- Stem that bears the GrassLeaf laterals --
stem_rp = pb.StemRandomParameter(plant)
stem_rp.subType = 1
stem_rp.lmax = 15.0
stem_rp.r = 2.0  # cm day⁻¹ elongation rate
stem_rp.la = 1.0
stem_rp.lb = 2.0  # first leaf at 2 cm
stem_rp.ln = 5.0  # internode spacing 5 cm
stem_rp.lnf = 0  # homogeneous distances
stem_rp.theta = 0.0
stem_rp.successor = [[1]]
stem_rp.successorP = [[1.0]]
stem_rp.successorOT = [[pb.leaf]]
stem_rp.successorNo = [1]
plant.setOrganRandomParameter(stem_rp)

# -- GrassLeaf random parameters (deterministic: std = 0) --
gl_rp = pb.GrassLeafRandomParameter(plant)
gl_rp.subType = 1
gl_rp.a = 0.02
gl_rp.bladeAngle = 0.4
gl_rp.bladeAngles = 0.0
gl_rp.bladeLength = 12.0  # cm
gl_rp.bladeLengths = 0.0
gl_rp.bladeWidth = 0.8  # cm
gl_rp.bladeWidths = 0.0
gl_rp.sheathLength = 6.0  # cm
gl_rp.sheathLengths = 0.0
gl_rp.sheathDuration = 8.0  # days
gl_rp.sheathDurations = 0.0
gl_rp.bladeDelay = 1.0  # days after sheath complete
gl_rp.bladeDelays = 0.0
gl_rp.bladeDuration = 15.0  # days to full blade
gl_rp.bladeDurations = 0.0
gl_rp.f_gf = pb.LinearGrowth()
plant.setOrganRandomParameter(gl_rp)

# --------------------------------------------------------------------------- #
#  2.  Simulation settings
# --------------------------------------------------------------------------- #

total_days = 25.0
dt = 0.5  # days per step
steps = int(total_days / dt)

# Refresh the VTK window every N simulation steps (1 = every step).
# Increase this value to speed up the simulation at the cost of fewer frames.
render_every = 1

# Set to a base filename (e.g. "grassleaf_") to export JPEG frames, or
# leave as None to skip frame export.
save_frames = "grassleaf_"  # e.g. "grassleaf_"

# --------------------------------------------------------------------------- #
#  3.  Warm-start: simulate until the plant has at least one segment so that
#      AnimateRoots can determine scene bounds for the initial camera setup.
# --------------------------------------------------------------------------- #

plant.initialize(verbose=False)

# Simulate just past the point where the stem starts growing (lb = 2 cm,
# r = 2 cm/day → about 1 day).  Any positive segment count is sufficient.
warmup_days = 1.5
warmup_steps = int(warmup_days / dt)
print(f"Warm-up: simulating {warmup_days} days …")
for _ in range(warmup_steps):
    plant.simulate(dt, verbose=False)
elapsed = warmup_steps * dt

# --------------------------------------------------------------------------- #
#  4.  Set up the animation
# --------------------------------------------------------------------------- #

anim = AnimateRoots(plant)
anim.root_name = "age"
anim.plant = True  # include leaf surface polygons via plot_plant
if save_frames:
    anim.avi_name = save_frames

anim.start(axis="v")  # open window with oblique camera

# --------------------------------------------------------------------------- #
#  5.  Main simulation + animation loop
# --------------------------------------------------------------------------- #

remaining_steps = steps - warmup_steps
print(f"Animating {remaining_steps} steps ({remaining_steps * dt:.1f} days) …")
for i in range(remaining_steps):
    plant.simulate(dt, verbose=False)
    elapsed += dt
    if i % render_every == 0:
        anim.simtime = elapsed
        anim.update()

print(f"\nSimulation finished at t = {elapsed:.1f} days.")
print("Close the VTK window or press 'e' to exit.")

# --------------------------------------------------------------------------- #
#  6.  Enter the interactive event loop (keeps the window open)
# --------------------------------------------------------------------------- #

anim.run()
