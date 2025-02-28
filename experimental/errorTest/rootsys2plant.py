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
for filename in os.listdir(plant_path):
    if 'xml' in filename:
        print('\t',filename)
        try:
            p = pb.MappedPlant()
            p.readParameters(os.path.join(plant_path, filename))
            p.initialize()
            p.simulate(10)
        except Exception as e:
            errorPlant.append(filename)
            errorPlantbis.append(e)
            #raise Exception
        
print("testing files in rot folder")
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
print(errorPlantbis)
print("error with root parameter files:",errorRoot)

# plant_directory = os.fsencode(plant_path)
# root_directory = os.fsencode(root_path)

# fileswerrors = []

# for plant_file in os.listdir(plant_directory):
#     filename = os.fsdecode(plant_file)
#     if filename.endswith(".xml"):
#         try:
#             p = pb.Plant()
#             p.readParameters(os.path.join(plant_path, filename))
#             p.initialize()
#             p.writeParameters(filename)
#         except:
#             print("error occured for " + filename + "\n")
#             fileswerrors.append(filename)



# for root_file in os.listdir(root_directory):
#     filename = os.fsdecode(root_file)
#     if filename.endswith(".xml") and not os.path.isfile(filename):
#          print(filename)
#          p = pb.Plant()
#          p.readParameters(os.path.join(root_path, filename))
#          p.writeParameters(filename)


# with open("parameterfileswerrors.txt",'w') as f:
#     f.write('\n'.join(fileswerrors))
