"""
Checks parameter sets for plausibility 
"""

import sys; sys.path.append(".."); sys.path.append("../..");

import os
import numpy as np

import plantbox as pb

OrganParameterRanges = {
    "a": [0.01, 2.5], 
    "dx" : [0.1, 1],
    "dxMin" : [0., 0.1],
}

OrganParameterValues = {
    "organType": [pb.seed, pb.root, pb.stem, pb.leaf],
    "subType" : range(0, 100)
}

# std::vector<std::vector<int> > successorWhere = std::vector<std::vector<int>>(0, std::vector<int> (0, 0)); ///< Lateral types [1]
# std::vector<std::vector<int> > successorOT = std::vector<std::vector<int>>(0, std::vector<int> (0, 0)); ///< Lateral types [1]
# std::vector<std::vector<int> > successorST = std::vector<std::vector<int>>(0, std::vector<int> (0, 0)); ///< Lateral types [1]
# std::vector<std::vector<double>> successorP = std::vector<std::vector<double>>(0, std::vector<double> (0, 0)); ///< Probabilities of lateral type to emerge (sum of values == 1) [1]
# std::vector<int>  successorNo = std::vector<int>(0); ///< Lateral types [1]

def checkFolder(folder_path):

    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        if os.path.isfile(file_path):
            print("try to open:", file_path)
            checkFile(file_path)
            
def checkFile(file):

    results = {}

    plant = pb.Plant()
    plant.readParameters(file)
    plant.initialize()     

    try: 
        plant = pb.Plant()
        plant.readParameters(file)
        plant.initialize(False)        
    except:
        raise Exception("checkFile(): something went wrong opening file "+ file)

    results[pb.root] = {}
    orp = plant.getOrganRandomParameter(pb.root)
    for p in orp:
        if not (p.organType == pb.root):
             raise "checkFile(): organType is of wrong type {:g} (should be root={:g})".format(p.organType,pb.root)
        st = p.subType
        if st in results[pb.root]:
             raise "checkFile(): Root     subType {:g} is defined multiple times".dormat(p.subType)
        results[pb.root][st] = check_root(p)
    
    return results    
    
    

def check_root(p):
    check_organ(p)
    pass

def check_seed(p):
    check_organ(p)
    pass

def check_stem(p):
    check_organ(p)
    pass

def check_leaf(p):
    check_organ(p)
    pass

def check_organ(p, msg = ""):    
    pass_ = True
    pass_ = check_range_(p, OrganParameterRanges, msg) and pass_
    pass_ = check_values_(p,OrganParameterValues, msg) and pass_
    if len(p.name) < 2 : 
        print(msg, "Parameters for {:s}: parameter name is short".format(p.name)) 
        pass_ = False
    if len(p.name) > 20 : 
        print(msg, "Parameters for {:s}: parameter name is long".format(p.name)) 
        pass_ = False
    if p.dxMin <= p.dx/10 :
        print(msg, "Parameters for {:s}: dxMin ({:g}) should be at least 1/10 of dx ({:g})".format(p.name, p.dxMin, p.dx)) 
        pass_ = False
        
def check_range_(p, ranges, msg = ""):
    for name, r in ranges.items():
        v = p.getParameter(name)
        if v<r[0]:
            print(msg, "Parameters for {:s}: parameter '{:s}' has value {:g} which seems too small (limit is {:g})"
                  .format(p.name, name, v, r[0])) 
            return False
        if v>r[1]:
            print(msg, "Parameters for {:s}: parameter '{:s}' has value {:g} which seems too large (limit is {:g})"
                  .format(p.name, name, v, r[1])) 
            return False
    return True
 
def check_values_(p, ranges, msg = ""):
    for name, r in ranges.items():
        v = p.getParameter(name)
        if not v in r:
            print(msg, "Parameters for {:s}: parameter '{:s}' has invalid value {:g}, allowed values are "
                  .format(p.name, name, v), r) 
            return False
    return True
    
if __name__ == "__main__":
    
    #print(os.getcwd())
    checkFolder("../../modelparameter/structural/rootsystem")
    
    