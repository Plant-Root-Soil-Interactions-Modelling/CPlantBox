"""
Turtle3D demo – 3-D polyline visualisation
==========================================
Demonstrates the ``pb.Turtle3D`` class from ``mymath.h`` / PyPlantBox by
tracing several named paths, each built with a different combination of turtle
commands, and then plotting them in a shared 3-D matplotlib figure.

Turtle3D coordinate frame (default orientation)
------------------------------------------------
  H (heading)  = +x   – the direction the turtle moves on ``forward()``
  L (left)     = +y   – used as the rotation axis for pitch
  U (up)       = +z   – used as the rotation axis for yaw

Available commands
------------------
  t.forward(dist)        – move dist along H
  t.turnLeft(angle)      – yaw left  (rotate H/L around U)
  t.turnRight(angle)     – yaw right
  t.pitchUp(angle)       – pitch up  (rotate H/U around L)
  t.pitchDown(angle)     – pitch down
  t.rollLeft(angle)      – roll left  (rotate L/U around H)
  t.rollRight(angle)     – roll right
  t.getPosition()        – returns current Vector3d position
  t.heading() / .left() / .up() – current frame axes as Vector3d
  t.setPosition(p)       – teleport without changing orientation
  t.setFrame(m)          – replace the full 3×3 frame matrix
  str(t)                 – human-readable state string

Run from this directory:
    PYTHONPATH=../../src python polyline_3d_example.py
"""

import math
import sys

sys.path.insert(0, "../../src")

import matplotlib.pyplot as plt
import numpy as np
import plantbox as pb
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

# --------------------------------------------------------------------------- #
#  Helper
# --------------------------------------------------------------------------- #


def pos_to_tuple(v):
    """Convert a plantbox Vector3d to a plain Python (x, y, z) tuple."""
    return (v.x, v.y, v.z)


# --------------------------------------------------------------------------- #
#  1.  Straight line  –  just forward()
# --------------------------------------------------------------------------- #

anchorframe = pb.Matrix3d.ons(pb.Vector3d(0, 0, 1))
anchropos = pb.Vector3d(0, 0, 0)

print("Anchor frame")
print(anchorframe)
print("Pos")
print(anchropos)
t = pb.Turtle3D(anchropos, anchorframe)
# t.turnLeft(math.pi / 2)
t.pitchUp(math.pi / 2)
straight = [pos_to_tuple(t.getPosition())]
for _ in range(10):
    # print("Frame")
    # print(t.getFrame())
    t.forward(1.0)
    straight.append(pos_to_tuple(t.getPosition()))

print("Straight line:")
for p in straight:
    print(f"  {p}")


# --------------------------------------------------------------------------- #
#  2.  Planar arc  –  turnLeft() + forward()  (yaw only)
# --------------------------------------------------------------------------- #

t = pb.Turtle3D()
arc = [pos_to_tuple(t.getPosition())]
for _ in range(12):
    t.turnLeft(math.pi / 6)  # 30° yaw each step
    t.forward(2.0)
    arc.append(pos_to_tuple(t.getPosition()))

print("\nPlanar arc (turnLeft + forward):")
for p in arc:
    print(f"  {p}")


# # --------------------------------------------------------------------------- #
# #  3.  3-D helix  –  turnLeft() + pitchUp() + forward()
# # --------------------------------------------------------------------------- #

# t = pb.Turtle3D()
# helix = [pos_to_tuple(t.getPosition())]
# for _ in range(24):
#     t.turnLeft(math.pi / 6)  # 30° yaw
#     t.pitchUp(math.pi / 18)  # 10° pitch per step
#     t.forward(1.5)
#     helix.append(pos_to_tuple(t.getPosition()))

# print("\nHelix (turnLeft + pitchUp + forward):")
# for p in helix:
#     print(f"  {p}")


# # --------------------------------------------------------------------------- #
# #  4.  Grass-leaf shape  –  straight sheath, then pitchDown() into blade
# #
# #  We first reorient the turtle so it grows upward (H = +z) via pitchDown(π/2).
# # --------------------------------------------------------------------------- #

# t = pb.Turtle3D()
# t.pitchDown(math.pi / 2)  # H now points to +z (upward)

# grassleaf = [pos_to_tuple(t.getPosition())]

# # Sheath: grow straight up, 6 steps × 1 cm
# for _ in range(6):
#     t.forward(1.0)
#     grassleaf.append(pos_to_tuple(t.getPosition()))

# # Ligule: one pitchDown at the sheath/blade boundary
# t.pitchDown(math.pi / 5)  # ~36° bend at ligule

# # Blade: continue forward, 6 more steps
# for _ in range(6):
#     t.forward(1.0)
#     grassleaf.append(pos_to_tuple(t.getPosition()))

# print("\nGrass-leaf shape (pitchDown at ligule):")
# for p in grassleaf:
#     print(f"  {p}")


# # --------------------------------------------------------------------------- #
# #  5.  Root path  –  pitchDown() + forward() (growing downward)
# # --------------------------------------------------------------------------- #

# t = pb.Turtle3D()
# t.pitchUp(math.pi / 2)  # H now points to -z (downward)

# root = [pos_to_tuple(t.getPosition())]
# for i in range(8):
#     t.turnLeft(math.pi / 8)  # 22.5° yaw per step
#     t.pitchDown(math.pi / 20)  # gentle downward drift
#     t.forward(1.2)
#     root.append(pos_to_tuple(t.getPosition()))

# print("\nRoot path (pitchUp to start downward, then turnLeft + pitchDown):")
# for p in root:
#     print(f"  {p}")


# # --------------------------------------------------------------------------- #
# #  Plot all paths
# # --------------------------------------------------------------------------- #

POLYLINES = {
    "straight (forward)": straight,
    "planar arc (turnLeft)": arc,
    #     "helix (turnLeft+pitchUp)": helix,
    #     "grass leaf (pitchDown)": grassleaf,
    #     "root (pitchUp+turnLeft)": root,
}

COLOURS = ["steelblue", "darkorange", "purple", "green", "sienna"]

fig = plt.figure(figsize=(9, 7))
ax = fig.add_subplot(111, projection="3d")

for colour, (name, pts) in zip(COLOURS, POLYLINES.items()):
    nodes = np.array(pts, dtype=float)
    ax.plot(nodes[:, 0], nodes[:, 1], nodes[:, 2], "-o", color=colour, markersize=4, linewidth=2, label=name)
    ax.scatter(*nodes[0], color=colour, s=70, marker="^", zorder=5)  # base
    ax.scatter(*nodes[-1], color=colour, s=70, marker="*", zorder=5)  # tip

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.set_title("Turtle3D demo")
ax.legend(fontsize=8)
ax.set_aspect("equal")
plt.tight_layout()
plt.savefig("polyline_3d_example.png", dpi=150)
print("\nSaved polyline_3d_example.png")
plt.show()
