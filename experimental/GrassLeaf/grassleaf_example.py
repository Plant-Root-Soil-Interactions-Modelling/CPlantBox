"""
GrassLeaf example
=================
Demonstrates how to create and simulate a single GrassLeaf organ in isolation,
inspect its parameter classes, and visualise the resulting polyline with
matplotlib.

Run from this directory:
    python grassleaf_example.py
"""

import math

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp


class Poaceae(pb.Plant):
    """Plant subclass that creates GrassLeaf organs instead of the default Leaf."""

    def createLeaf(self, subType, delay, parent, pni):
        print("wild things are going on", flush=True)
        return pb.GrassLeaf(self, subType, delay, parent, pni)


# --------------------------------------------------------------------------- #
#  1.  Build a minimal Plant with one stem that hosts GrassLeaf laterals
# --------------------------------------------------------------------------- #

plant = Poaceae()  #  pb.Plant()

# -- Seed (required boilerplate) --
seed_rp = pb.SeedRandomParameter(plant)
seed_rp.subType = 0
plant.setOrganRandomParameter(seed_rp)

# -- Tap root (minimal, very short) --
root_rp = pb.RootRandomParameter(plant)
root_rp.subType = 1
root_rp.lmax = 0.1
root_rp.r = 1.0
root_rp.theta = 0.0
plant.setOrganRandomParameter(root_rp)

# -- Stem: carries the GrassLeaf laterals --
stem_rp = pb.StemRandomParameter(plant)
stem_rp.subType = 1
stem_rp.lmax = 15.0
stem_rp.r = 2.0  # 2 cm/day elongation rate
stem_rp.la = 1.0
stem_rp.lb = 2.0  # first leaf at 2 cm
stem_rp.ln = 5.0  # internode 5 cm
stem_rp.lnf = 0  # homogeneous distances
stem_rp.theta = 0.0
stem_rp.successor = [[1]]  # GrassLeaf subType 1
stem_rp.successorP = [[1.0]]
stem_rp.successorOT = [[pb.leaf]]
stem_rp.successorNo = [1]  # one leaf per stem
plant.setOrganRandomParameter(stem_rp)

# print(stem_rp.successorNo, flush=True)
# print("Stem", int(pb.stem), "Leaf", int(pb.leaf), flush=True)


# -- GrassLeaf random parameters --
gl_rp = pb.GrassLeafRandomParameter(plant)
gl_rp.subType = 1
gl_rp.a = 0.1
gl_rp.bladeAngle = 0.4  # ~23 deg bend at ligule
gl_rp.bladeAngles = 0.05
gl_rp.bladeLength = 12.0  # cm
gl_rp.bladeLengths = 1.0
gl_rp.bladeWidth = 0.8  # cm
gl_rp.bladeWidths = 0.05
gl_rp.sheathLength = 6.0  # cm
gl_rp.sheathLengths = 2
gl_rp.sheathDuration = 8.0  # days
gl_rp.sheathDurations = 0.0
gl_rp.bladeDelay = 1.0  # days after sheath complete
gl_rp.bladeDelays = 0.0
gl_rp.bladeDuration = 15.0  # days to full blade
gl_rp.bladeDurations = 1.0
gl_rp.f_gf = pb.LinearGrowth()  # for other organs this is set in initCallbacks() from parameters

# print(gl_rp.successorNo, flush=True)
# gl_rp = pb.LeafRandomParameter(plant)
# gl_rp.lmax = 12.0
# gl_rp.subType = 1
plant.setOrganRandomParameter(gl_rp)


# print(gl_rp.f_gf)
# print(gl_rp.f_gf.getLength(5.0, 1.0, 1.0, None))  # example call to growth function
# print("done")
# ss

# --------------------------------------------------------------------------- #
#  2.  Initialise and run the simulation
# --------------------------------------------------------------------------- #

plant.initialize(verbose=False)

total_days = 100.0
dt = 0.5  # days per step
steps = int(total_days / dt)

print(f"Simulating {total_days} days in {steps} steps …")
for i in range(steps):
    print(i, end=" ", flush=True)
    plant.simulate(dt, verbose=False)


# --------------------------------------------------------------------------- #
#  3.  Inspect the GrassLeaf organs
# --------------------------------------------------------------------------- #

# organs = plant.getOrgans(pb.OrganTypes.leaf, True)
# print(f"\nFound {len(organs)} leaf organ(s)")

# for gl in organs:

#     p = gl.param()
#     print("\n--- GrassLeaf ---")
#     print(f"  age            = {gl.getAge():.2f} days")
#     print(f"  total length   = {gl.getLength():.2f} cm")
#     print(f"  sheathLength   = {gl.getSheathLength():.2f} / {p.sheathLength:.2f} cm")
#     print(f"  bladeLengthGrown = {gl.getBladeLengthGrown():.2f} / {p.bladeLength:.2f} cm")
#     print(f"  nodes          = {gl.getNumberOfNodes()}")
#     print(gl.getParent().getNode(gl.parentNI))  # parent node where leaf is attached
#     for i in range(0, gl.getNumberOfNodes()):
#         print(gl.getNode(i), end="; ")


ana = pb.SegmentAnalyser(plant)
ana.addAge(total_days)
# vp.plot_roots(ana, "age")  # plot roots only |\label{l13:plot_roots}|
vp.plot_roots(ana, "organType")

# --------------------------------------------------------------------------- #
#  4.  Step-by-step time series for ONE leaf (first one found)
# --------------------------------------------------------------------------- #

# Re-run a fresh plant, record state every 2 days
plant2 = Poaceae()
plant2.setOrganRandomParameter(pb.SeedRandomParameter(plant2))
rp2 = pb.RootRandomParameter(plant2)
rp2.subType = 1
rp2.lmax = 0.1
rp2.r = 1.0
rp2.theta = 0.0
plant2.setOrganRandomParameter(rp2)
s2 = pb.StemRandomParameter(plant2)
s2.subType = 1
s2.lmax = 30.0
s2.r = 2.0
s2.la = 0.0
s2.lb = 2.0
s2.ln = 5.0
s2.lnf = 1
s2.theta = 0.0
# s2.successor = [1]
# s2.successorP = [1.0]
plant2.setOrganRandomParameter(s2)
g2 = pb.GrassLeafRandomParameter(plant2)
g2.subType = 1
g2.bladeAngle = 0.4
g2.bladeAngles = 0.0
g2.bladeLength = 12.0
g2.bladeLengths = 0.0
g2.bladeWidth = 0.8
g2.bladeWidths = 0.0
g2.sheathLength = 6.0
g2.sheathLengths = 0.0
g2.sheathDuration = 8.0
g2.sheathDurations = 0.0
g2.bladeDelay = 1.0
g2.bladeDelays = 0.0
g2.bladeDuration = 15.0
g2.bladeDurations = 0.0
plant2.setOrganRandomParameter(g2)
plant2.initialize(verbose=False)

times, sheath_vals, blade_vals = [], [], []
record_dt = 0.5
for step in range(int(total_days / record_dt)):
    plant2.simulate(record_dt, verbose=False)
    leaves = plant2.getOrgans(pb.OrganTypes.leaf)
    if leaves:
        gl = leaves[0]
        times.append((step + 1) * record_dt)
        sheath_vals.append(gl.getSheathLength())
        blade_vals.append(gl.getBladeLengthGrown())

# --------------------------------------------------------------------------- #
#  5.  Visualise: growth curves
# --------------------------------------------------------------------------- #

fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(times, sheath_vals, label="Sheath length grown [cm]")
ax.plot(times, blade_vals, label="Blade length grown [cm]")
ax.set_xlabel("Time [days]")
ax.set_ylabel("Length [cm]")
ax.set_title("GrassLeaf growth over time (first leaf)")
ax.legend()
ax.grid(True)
plt.tight_layout()
plt.savefig("grassleaf_growth.png", dpi=150)
print("\nGrowth curve saved to grassleaf_growth.png")

# --------------------------------------------------------------------------- #
#  6.  Visualise: 3-D polyline of the first fully grown leaf
# --------------------------------------------------------------------------- #

leaves_final = plant.getOrgans(pb.OrganTypes.leaf)
if leaves_final:
    gl = leaves_final[0]
    n = gl.getNumberOfNodes()
    nodes = np.array([[gl.getNode(i).x, gl.getNode(i).y, gl.getNode(i).z] for i in range(n)])

    fig3d = plt.figure(figsize=(6, 8))
    ax3d = fig3d.add_subplot(111, projection="3d")
    ax3d.plot(nodes[:, 0], nodes[:, 1], nodes[:, 2], "g-o", markersize=4, linewidth=2)
    ax3d.scatter(*nodes[0], color="blue", s=60, zorder=5, label="base (node 0)")
    ax3d.scatter(*nodes[-1], color="red", s=60, zorder=5, label="tip")
    # mark sheath / blade boundary
    n_sheath = int(gl.getSheathLength() / max(gl.length, 1e-9) * (n - 1)) + 1
    if 0 < n_sheath < n:
        ax3d.scatter(*nodes[n_sheath - 1], color="orange", s=80, zorder=5, label="ligule (approx.)")
    ax3d.set_xlabel("x [cm]")
    ax3d.set_ylabel("y [cm]")
    ax3d.set_zlabel("z [cm]")
    ax3d.set_title("GrassLeaf polyline (3D)")
    ax3d.legend()
    plt.tight_layout()
    plt.savefig("grassleaf_3d.png", dpi=150)
    print("3-D polyline saved to  grassleaf_3d.png")

# --------------------------------------------------------------------------- #
#  7.  Direct use of Turtle3D and Meristem (unit-level test)
# --------------------------------------------------------------------------- #

print("\n--- Turtle3D / Meristem standalone ---")
t = pb.Turtle3D()
t.forward(5.0)
t.pitchDown(math.radians(30))
t.forward(8.0)
print("Turtle position after forward(5)+pitchDown(30°)+forward(8):", t.getPosition())

m = pb.Meristem()
m.addNodeBack(5.0)  # straight segment
m.addNodeBack(8.0, 0.0, math.radians(30))  # pitch down 30°
print(f"Meristem size: {m.size()} nodes")
poly = m.getPolyline()
print("Polyline tip:", poly[-1])

plt.show()
print("\nDone.")
