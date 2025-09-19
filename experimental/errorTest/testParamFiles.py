"""converts old xml files based on the RootSystem Class to be based on the new(er) Plant Class"""
import os
import sys; sys.path.append("../.."); sys.path.append("../../src/")
import plantbox as pb


plant_path = "../../modelparameter/structural/plant"
root_path = "../../modelparameter/structural/rootsystem"
errorPlant = []
errorRoot = []
errorPlantbis=[]
print("testing files in plant folder")

for filename in os.listdir(plant_path): # ["Juncus_squarrosus_Clausnitzer_199420.xml"]:#
    if 'xml' in filename:
        print('\t',filename)
        try:
            p = pb.MappedPlant()
            p.readParameters(os.path.join(plant_path, filename))
            p.initialize()
            p.simulate(1)
        except Exception as e:
            errorPlant.append(filename)
            errorPlantbis.append(e)
            #raise Exception
            
for filename in os.listdir(root_path):
    if 'xml' in filename:
        print('\t',filename)
        try:
            p = pb.MappedPlant()
            p.readParameters(os.path.join(root_path, filename))
            p.initialize()
            p.simulate(10)
        except:
            errorRoot.append(filename)
        
print("error with plant parameter files:",errorPlant)
if len(errorPlant)> 0:
    print("error messages:",errorPlantbis)
print("error with root parameter files:",errorRoot)