"""
Plant organ overview
====================
Simulates the same Poaceae plant as grassleaf_example.py and then prints a
table with one row per organ showing:

  - organ index
  - organ type name  (seed / root / stem / leaf)
  - organ sub-type
  - organ ID         (global id assigned by CPlantBox)
  - number of nodes
  - number of segments
  - list of global node IDs

Run from this directory:
    python plant_organ_overview.py
"""

import plantbox as pb


# --------------------------------------------------------------------------- #
#  Organ-type helpers
# --------------------------------------------------------------------------- #

_OT_NAMES = {
    int(pb.OrganTypes.organ): "organ",
    int(pb.OrganTypes.seed):  "seed",
    int(pb.OrganTypes.root):  "root",
    int(pb.OrganTypes.stem):  "stem",
    int(pb.OrganTypes.leaf):  "leaf",
}


def ot_name(organ):
    """Return the human-readable organ-type string for *organ*."""
    ot = int(organ.getParameter("organType"))
    return _OT_NAMES.get(ot, f"ot{ot}")


# --------------------------------------------------------------------------- #
#  1.  Rebuild the same plant as in grassleaf_example.py
# --------------------------------------------------------------------------- #

class Poaceae(pb.Plant):
    def createLeaf(self, subType, delay, parent, pni):
        return pb.GrassLeaf(self, subType, delay, parent, pni)


plant = Poaceae()

seed_rp = pb.SeedRandomParameter(plant)
seed_rp.subType = 0
plant.setOrganRandomParameter(seed_rp)

root_rp = pb.RootRandomParameter(plant)
root_rp.subType = 1
root_rp.lmax = 0.1
root_rp.r = 1.0
root_rp.theta = 0.0
plant.setOrganRandomParameter(root_rp)

stem_rp = pb.StemRandomParameter(plant)
stem_rp.subType = 1
stem_rp.lmax = 15.0
stem_rp.r = 2.0
stem_rp.la = 1.0
stem_rp.lb = 2.0
stem_rp.ln = 5.0
stem_rp.lnf = 0
stem_rp.theta = 0.0
stem_rp.successor = [[1]]
stem_rp.successorP = [[1.0]]
stem_rp.successorOT = [[pb.leaf]]
stem_rp.successorNo = [1]
plant.setOrganRandomParameter(stem_rp)

gl_rp = pb.GrassLeafRandomParameter(plant)
gl_rp.subType = 1
gl_rp.a = 0.1
gl_rp.bladeAngle = 0.4
gl_rp.bladeAngles = 0.05
gl_rp.bladeLength = 12.0
gl_rp.bladeLengths = 1.0
gl_rp.bladeWidth = 0.8
gl_rp.bladeWidths = 0.05
gl_rp.sheathLength = 6.0
gl_rp.sheathLengths = 2.0
gl_rp.leafGrowthDuration = 20.0  # days to full leaf (sheath + blade)
gl_rp.leafGrowthDurations = 1.0
gl_rp.f_gf = pb.LinearGrowth()
plant.setOrganRandomParameter(gl_rp)

plant.initialize(verbose=False)

total_days = 100.0
dt = 0.5
steps = int(total_days / dt)

print(f"Simulating {total_days} days …", flush=True)
for _ in range(steps):
    plant.simulate(dt, verbose=False)

# --------------------------------------------------------------------------- #
#  2.  Collect all organs and print the overview table
# --------------------------------------------------------------------------- #

organs = plant.getOrgans()   # all organ types, all organs

col_widths = (5, 8, 8, 9, 8, 12, 0)   # last col is unbounded
header = (
    f"{'#':>5}  "
    f"{'type':<8}  "
    f"{'subType':>8}  "
    f"{'organID':>9}  "
    f"{'nodes':>8}  "
    f"{'segments':>12}  "
    f"node IDs"
)
separator = "-" * max(80, len(header))

print(f"\nPlant organ overview  ({len(organs)} organs total)\n")
print(header)
print(separator)

for idx, organ in enumerate(organs):
    n_nodes = organ.getNumberOfNodes()
    n_segs  = organ.getNumberOfSegments()
    node_ids = list(organ.getNodeIds())

    print(
        f"{idx:>5}  "
        f"{ot_name(organ):<8}  "
        f"{int(organ.getParameter('subType')):>8}  "
        f"{organ.getId():>9}  "
        f"{n_nodes:>8}  "
        f"{n_segs:>12}  "
        f"{node_ids}"
    )

print(separator)
print(
    f"Total nodes across all organs : "
    f"{sum(o.getNumberOfNodes() for o in organs)}"
)
print(
    f"Total segments across all organs: "
    f"{sum(o.getNumberOfSegments() for o in organs)}"
)
