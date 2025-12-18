"""converts old xml files based on the RootSystem Class to be based on the new(er) Plant Class"""
import os
from os import walk
import sys; sys.path.append("../.."); sys.path.append("../../src/")
import plantbox as pb


plant_path = "plant"
root_path = "rootsystem"
# filename = "Anagallis_femina_Leitner_2010"

# plant = pb.Plant()
# plant.readParameters(root_path +"/"+ filename +".xml",verbose = True)
# plant.writeParameters(filename +".xml")


plant_directory = os.fsencode(plant_path)
root_directory = os.fsencode(root_path)

for plant_file in os.listdir(plant_directory):
    filename = os.fsdecode(plant_file)
    if filename.endswith(".xml") and not os.path.isfile(filename):
         print(filename)
         p = pb.Plant()
         p.readParameters(os.path.join(plant_path, filename))
         p.writeParameters(os.path.join(plant_path, filename))


for root_file in os.listdir(root_directory):
    filename = os.fsdecode(root_file)
    if filename.endswith(".xml") and not os.path.isfile(filename):
         print(filename)
         p = pb.Plant()
         p.readParameters(os.path.join(root_path, filename))
         p.writeParameters(os.path.join(plant_path, filename))