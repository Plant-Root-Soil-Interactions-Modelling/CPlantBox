"""Dump the plant state to a json file, restore it into a new plant, and compare the continued growth"""

import json

import matplotlib.pyplot as plt
import numpy as np

import plantbox as pb

path = pb.data_path() + "/structural/plant/"
name = "fspm2023"
param_file = path + name + ".xml"

ptypes = {pb.seed: pb.SeedSpecificParameter, pb.root: pb.RootSpecificParameter, pb.stem: pb.StemSpecificParameter, pb.leaf: pb.LeafSpecificParameter}
otypes = {pb.root: pb.Root, pb.stem: pb.Stem, pb.leaf: pb.Leaf}


def all_organs(plant):
    """The seed and all organs below it (including organs that have not started to grow yet)"""
    result, stack = [], [plant.getSeed()]
    while stack:
        o = stack.pop()
        result.append(o)
        stack.extend(o.getChild(i) for i in range(o.getNumberOfChildren()))
    return sorted(result, key=lambda o: o.getId())


def dump_plant(plant, file_name):
    """Writes the state of the organism and of all organs (linked by organ ids) into a json file"""
    data = {"organism": plant.getState(), "organs": {str(o.getId()): o.getState() for o in all_organs(plant)}}
    with open(file_name, "w") as f:
        json.dump(data, f)


def load_plant(file_name):
    """Creates a new plant from a json file written by dump_plant"""
    with open(file_name) as f:
        data = json.load(f)
    plant = pb.Plant()
    plant.readParameters(param_file)  # the random parameters are not part of the state
    plant.initialize(False)  # sets up tropisms, growth functions, and the parameters the seed derives; the organs it creates are replaced below
    organs = {}
    for key, s in data["organs"].items():  # 1. create empty organs
        ot, oid = int(s["organType"][0]), int(s["id"][0])
        if ot == pb.seed:
            o = pb.Seed(oid, ptypes[ot](), True, True, 0.0, 0.0, False, 0)
        else:
            o = otypes[ot](oid, ptypes[ot](), True, True, 0.0, 0.0, pb.Vector3d(0, 0, 1), 0, False, 0)
        o.setOrganism(plant)
        organs[oid] = o
    for oid, o in organs.items():  # 2. set their state, and link them using the ids
        o.setState(data["organs"][str(oid)], organs)
    plant.setState(data["organism"], organs)  # 3. counters, random number generators, base organs
    return plant


def plot_plant(ax, plant, title, color):
    nodes = np.array([[n.x, n.y, n.z] for n in plant.getNodes()])
    for s in plant.getSegments():
        p = nodes[[s.x, s.y]]
        ax.plot(p[:, 0], p[:, 1], p[:, 2], color=color, lw=0.7)
    ax.set_title(title)
    ax.set_xlabel("x (cm)")
    ax.set_ylabel("y (cm)")
    ax.set_zlabel("z (cm)")



t1, t2 = 20, 20  # days before the dump, and days after the dump

# Plant 1: grow, dump, continue growing
plant1 = pb.Plant(1)
plant1.readParameters(param_file)
plant1.initialize(False)
plant1.simulate(t1, False)
dump_plant(plant1, "results/example_state.json")
plant1.simulate(t2, False)

# Plant 2: new plant from the dumped state, continue growing
plant2 = load_plant("results/example_state.json")
plant2.simulate(t2, False)

nodes1 = np.array([[n.x, n.y, n.z] for n in plant1.getNodes()])
nodes2 = np.array([[n.x, n.y, n.z] for n in plant2.getNodes()])
print("plant 1: {} organs, {} nodes".format(plant1.getNumberOfOrgans(), len(nodes1)))
print("plant 2: {} organs, {} nodes".format(plant2.getNumberOfOrgans(), len(nodes2)))
print("maximal node difference: {:g} cm".format(np.abs(nodes1 - nodes2).max()))

# Plot 1: the two plants side by side
fig = plt.figure(figsize=(10, 6))
ax1 = fig.add_subplot(1, 2, 1, projection="3d")
ax2 = fig.add_subplot(1, 2, 2, projection="3d")
plot_plant(ax1, plant1, "Plant 1 (continuous growth)", "tab:green")
plot_plant(ax2, plant2, "Plant 2 (restored from state)", "tab:blue")
lim = np.array([ax.get_xlim() for ax in (ax1, ax2)] + [ax.get_ylim() for ax in (ax1, ax2)])
for ax in (ax1, ax2):  # same view
    ax.set_xlim(lim.min(), lim.max())
    ax.set_ylim(lim.min(), lim.max())
    ax.set_zlim(*ax1.get_zlim())
plt.tight_layout()
plt.savefig("results/example_state_plants.png", dpi=150)

# Plot 2: node coordinates of plant 2 against plant 1 (all points on the diagonal if both plants are the same)
fig, ax = plt.subplots(figsize=(5, 5))
for i, c in enumerate("xyz"):
    ax.plot(nodes1[:, i], nodes2[:, i], ".", ms=3, label=c)
ax.plot(ax.get_xlim(), ax.get_xlim(), "k--", lw=0.8)
ax.set_xlabel("plant 1 (cm)")
ax.set_ylabel("plant 2 (cm)")
ax.set_title("Node coordinates, max diff {:g} cm".format(np.abs(nodes1 - nodes2).max()))
ax.legend()
plt.tight_layout()
plt.savefig("results/example_state_compare.png", dpi=150)
plt.show()
