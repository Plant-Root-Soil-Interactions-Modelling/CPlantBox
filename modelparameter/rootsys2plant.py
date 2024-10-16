"""converts old xml files based on the RootSystem Class to be based on the new(er) Plant Class"""
import os
import sys; sys.path.append(".."); sys.path.append("../src/")
import plantbox as pb


plant_path = "structural/plant"
root_path = "structural/rootsystem"
filename = "Bench_lupin"

plant = pb.Plant()
plant.readParameters(root_path +"/"+ filename +".xml",verbose = True)
plant.writeParameters("testOfWrite.xml")


plant_directory = os.fsencode(plant_path)
# root_directory = os.fsencode(root_path)

# for plant_file in os.listdir(plant_directory):
#     filename = os.fsdecode(plant_file)
#     if filename.endswith(".xml") and not os.path.isfile(filename):
#          print(filename)
#          p = pb.Plant()
#          p.readParameters(os.path.join(plant_path, filename))
#          p.writeParameters(filename)


# for root_file in os.listdir(root_directory):
#     filename = os.fsdecode(root_file)
#     if filename.endswith(".xml") and not os.path.isfile(filename):
#          print(filename)
#          p = pb.Plant()
#          p.readParameters(os.path.join(root_path, filename))
#          p.writeParameters(filename)

