"""
Checks parameter sets for plausibility
"""

import os
import sys

import numpy as np

import plantbox as pb

OrganParameterRanges = {
    "a": [0.01, 2.5],
    "dx": [0.1, 1],
    "dxMin": [0.0, 0.1],
}

OrganParameterValues = {"organType": [pb.seed, pb.root, pb.stem, pb.leaf], "subType": range(0, 100)}

RootParameterRanges = {
    "la": [0.0, 50.0],  # apical zone [cm]
    "lb": [0.0, 50.0],  # basal zone [cm]
    "ln": [0.0, 30.0],  # inter-lateral distance [cm]
    "lmax": [0.0, 500.0],  # maximal root length [cm]
    "r": [0.001, 50.0],  # initial growth rate [cm/day]
    "rlt": [1.0, 1e12],  # root life time [day]
    "theta": [0.0, np.pi],  # insertion angle [rad]
    "tropismN": [0.0, 30.0],  # number of tropism trials
    "tropismS": [0.0, 2 * np.pi],  # tropism sensitivity [1/cm]
}
RootParameterValues = {
    "tropismT": [0, 1, 2, 3, 4, 5],
    "gf": [1, 2],
}

StemParameterRanges = {
    "la": [0.0, 50.0],
    "lb": [0.0, 50.0],
    "ln": [0.0, 50.0],
    "lmax": [0.0, 1000.0],
    "r": [0.001, 50.0],
    "rlt": [1.0, 1e12],
    "theta": [0.0, np.pi / 2],  # angle from parent stem [rad]
    "tropismN": [0.0, 30.0],
    "tropismS": [0.0, 2 * np.pi],
}
StemParameterValues = {
    "tropismT": [0, 1, 2, 3, 4, 5, 6],
    "gf": [1, 2],
    "lnf": [0, 1, 2, 3, 4],
    "nodalGrowth": [0, 1],
}

LeafParameterRanges = {
    "la": [0.0, 100.0],
    "lb": [0.0, 100.0],
    "lmax": [0.0, 200.0],
    "r": [0.001, 50.0],
    "rlt": [1.0, 1e12],
    "theta": [0.0, np.pi],
    "areaMax": [0.001, 2000.0],  # leaf area [cm²]
    "tropismN": [0.0, 30.0],
    "tropismS": [0.0, 2 * np.pi],
}
LeafParameterValues = {
    "tropismT": [0, 1, 2, 3, 4, 5, 6],
    "gf": [1, 2],
    "shapeType": [0, 1, 2],
    "parametrisationType": [0, 1],
}

SeedParameterRanges = {
    "maxB": [0.0, 500.0],  # max basal roots
    "firstB": [0.0, 365.0],  # first basal root [day]
    "delayB": [0.0, 365.0],  # delay between basal roots [day]
}

# std::vector<std::vector<int> > successorWhere = std::vector<std::vector<int>>(0, std::vector<int> (0, 0)); ///< Lateral types [1]
# std::vector<std::vector<int> > successorOT = std::vector<std::vector<int>>(0, std::vector<int> (0, 0)); ///< Lateral types [1]
# std::vector<std::vector<int> > successorST = std::vector<std::vector<int>>(0, std::vector<int> (0, 0)); ///< Lateral types [1]
# std::vector<std::vector<double>> successorP = std::vector<std::vector<double>>(0, std::vector<double> (0, 0)); ///< Probabilities of lateral type to emerge (sum of values == 1) [1]
# std::vector<int>  successorNo = std::vector<int>(0); ///< Lateral types [1]


def checkFolder(folder_path):
    """Runs checkFile() on every file in folder_path, printing the path before each attempt."""
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        if os.path.isfile(file_path):
            print("try to open:", file_path)
            checkFile(file_path)


def checkFile(file):
    """
    Load a CPlantBox XML parameter file, initialise the plant, and run plausibility
    checks on all organ types (seed, root, stem, leaf).

    Returns a nested dict  results[organType][subType] -> bool (True = passed).
    Raises an Exception if the file cannot be opened or an organType mismatch is found.
    """
    results = {}

    try:
        plant = pb.Plant()
        plant.readParameters(file)
        plant.initialize(False)
    except Exception as e:
        raise Exception("checkFile(): could not open file {:s}: {:s}".format(file, str(e)))

    organ_names = {pb.seed: "Seed", pb.root: "Root", pb.stem: "Stem", pb.leaf: "Leaf"}
    checkers = {pb.root: check_root, pb.stem: check_stem, pb.leaf: check_leaf, pb.seed: check_seed}

    for ot in [pb.seed, pb.root, pb.stem, pb.leaf]:
        results[ot] = {}
        for p in plant.getOrganRandomParameter(ot):
            if p.organType != ot:
                raise Exception("checkFile(): organType mismatch {:g} (expected {:s}={:g})".format(p.organType, organ_names[ot], ot))
            st = p.subType
            if st in results[ot]:
                raise Exception("checkFile(): {:s} subType {:d} is defined multiple times".format(organ_names[ot], st))
            results[ot][st] = checkers[ot](p)

    return results


def check_root(p):
    """
    Check plausibility of a RootRandomParameter set p.

    In addition to the shared organ checks, verifies:
    - la + lb <= lmax  (apical + basal zones must fit within total length)
    - lmaxs <= lmax    (std dev must not allow near-zero or negative lengths)
    - rs <= r          (std dev must not allow negative growth rates)

    Returns True if all checks pass, False otherwise.
    """
    name = "Root '{:s}' (subType {:d})".format(p.name, p.subType)
    ok = check_organ(p, name)
    ok = check_range_(p, RootParameterRanges, name) and ok
    ok = check_values_(p, RootParameterValues, name) and ok
    if p.la + p.lb > p.lmax:
        print(name + ": la ({:g}) + lb ({:g}) = {:g} exceeds lmax ({:g})".format(p.la, p.lb, p.la + p.lb, p.lmax))
        ok = False
    if p.lmaxs > p.lmax > 0:
        print(name + ": lmax std dev ({:g}) exceeds lmax ({:g}), may produce near-zero lengths".format(p.lmaxs, p.lmax))
        ok = False
    if p.rs > p.r > 0:
        print(name + ": r std dev ({:g}) exceeds r ({:g}), may produce negative growth rates".format(p.rs, p.r))
        ok = False
    return ok


def check_seed(p):
    """
    Check plausibility of a SeedRandomParameter set p.

    In addition to the shared organ checks, verifies:
    - seedPos.z <= 0   (seed must be at or below the soil surface)
    - firstB <= simtime when maxB > 0  (at least one basal root must emerge within the simulation)

    Returns True if all checks pass, False otherwise.
    """
    name = "Seed '{:s}' (subType {:d})".format(p.name, p.subType)
    ok = check_organ(p, name)
    ok = check_range_(p, SeedParameterRanges, name) and ok
    if p.seedPos.z > 0:
        print(name + ": seedPos.z ({:g}) is above the soil surface (expected <= 0)".format(p.seedPos.z))
        ok = False
    if p.maxB > 0 and p.firstB > p.simtime:
        print(name + ": firstB ({:g}) exceeds simtime ({:g}), no basal roots will emerge".format(p.firstB, p.simtime))
        ok = False
    return ok


def check_stem(p):
    """
    Check plausibility of a StemRandomParameter set p.

    In addition to the shared organ checks, verifies:
    - la + lb <= lmax  (apical + basal zones must fit within total length)
    - lmaxs <= lmax    (std dev must not allow near-zero or negative lengths)
    - rs <= r          (std dev must not allow negative growth rates)

    Returns True if all checks pass, False otherwise.
    """
    name = "Stem '{:s}' (subType {:d})".format(p.name, p.subType)
    ok = check_organ(p, name)
    ok = check_range_(p, StemParameterRanges, name) and ok
    ok = check_values_(p, StemParameterValues, name) and ok
    if p.la + p.lb > p.lmax:
        print(name + ": la ({:g}) + lb ({:g}) = {:g} exceeds lmax ({:g})".format(p.la, p.lb, p.la + p.lb, p.lmax))
        ok = False
    if p.lmaxs > p.lmax > 0:
        print(name + ": lmax std dev ({:g}) exceeds lmax ({:g}), may produce near-zero lengths".format(p.lmaxs, p.lmax))
        ok = False
    if p.rs > p.r > 0:
        print(name + ": r std dev ({:g}) exceeds r ({:g}), may produce negative growth rates".format(p.rs, p.r))
        ok = False
    return ok


def check_leaf(p):
    """
    Check plausibility of a LeafRandomParameter set p.

    In addition to the shared organ checks, verifies:
    - la + lb <= lmax      (petiole + apical zones must fit within total length)
    - lmaxs <= lmax        (std dev must not allow near-zero or negative lengths)
    - leafGeometryX present when shapeType == 2  (2D shape requires geometry data)

    Returns True if all checks pass, False otherwise.
    """
    name = "Leaf '{:s}' (subType {:d})".format(p.name, p.subType)
    ok = check_organ(p, name)
    ok = check_range_(p, LeafParameterRanges, name) and ok
    ok = check_values_(p, LeafParameterValues, name) and ok
    if p.la + p.lb > p.lmax:
        print(name + ": la ({:g}) + lb ({:g}) = {:g} exceeds lmax ({:g}), no blade zone remains".format(p.la, p.lb, p.la + p.lb, p.lmax))
        ok = False
    if p.lmaxs > p.lmax > 0:
        print(name + ": lmax std dev ({:g}) exceeds lmax ({:g}), may produce near-zero lengths".format(p.lmaxs, p.lmax))
        ok = False
    if p.shapeType == 2 and len(p.leafGeometryX) == 0:
        print(name + ": shapeType=2 (2D user-defined) but no leafGeometry is defined")
        ok = False
    return ok


def check_organ(p, msg=""):
    """
    Run checks that apply to every organ type (root, stem, leaf, seed).

    Checks OrganParameterRanges (a, dx, dxMin) and OrganParameterValues
    (organType, subType), then validates:
    - name length between 2 and 20 characters
    - dxMin >= dx/10  (minimum segment size should not be excessively small relative to dx)

    Returns True if all checks pass, False otherwise.
    """
    pass_ = True
    pass_ = check_range_(p, OrganParameterRanges, msg) and pass_
    pass_ = check_values_(p, OrganParameterValues, msg) and pass_
    if len(p.name) < 2:
        print(msg, "Parameters for {:s}: parameter name is short".format(p.name))
        pass_ = False
    if len(p.name) > 20:
        print(msg, "Parameters for {:s}: parameter name is long".format(p.name))
        pass_ = False
    if p.dxMin <= p.dx / 10:
        print(msg, "Parameters for {:s}: dxMin ({:g}) should be at least 1/10 of dx ({:g})".format(p.name, p.dxMin, p.dx))
        pass_ = False
    return pass_


def check_range_(p, ranges, msg=""):
    """
    Verify that each parameter listed in *ranges* falls within [min, max].

    ranges: dict mapping parameter name -> [min, max]
    Returns True if every parameter is within bounds, False on first violation.
    """
    for name, r in ranges.items():
        v = p.getParameter(name)
        if v < r[0]:
            print(msg, "Parameters for {:s}: parameter '{:s}' has value {:g} which seems too small (limit is {:g})".format(p.name, name, v, r[0]))
            return False
        if v > r[1]:
            print(msg, "Parameters for {:s}: parameter '{:s}' has value {:g} which seems too large (limit is {:g})".format(p.name, name, v, r[1]))
            return False
    return True


def check_values_(p, ranges, msg=""):
    """
    Verify that each parameter listed in *ranges* holds one of the allowed values.

    ranges: dict mapping parameter name -> list/range of valid values
    Returns True if every parameter is valid, False on first violation.
    """
    for name, r in ranges.items():
        v = p.getParameter(name)
        if not v in r:
            print(msg, "Parameters for {:s}: parameter '{:s}' has invalid value {:g}, allowed values are ".format(p.name, name, v), r)
            return False
    return True


if __name__ == "__main__":

    # print(os.getcwd())
    checkFolder("../../modelparameter/structural/rootsystem")
