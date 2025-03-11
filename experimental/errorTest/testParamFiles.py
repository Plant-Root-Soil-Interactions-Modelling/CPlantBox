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
#['small_2020_low_xyl.xml', 'small_2020_abs.xml', 'small_2020.xml', 'small_2020_low_st.xml', 'new_tree.xml', 'small_2020_more_lateral.xml', 'small_2020_l.xml', 'morning_glory.xml', 'small_2020_high_xyl_x_10.xml', 'leaf_opposite_decussate.xml', 'testLatDef.xml', 'x.xml', 'carbon2020.xml','morning_glory_3m_d.xml', 'small_2020_low_xyl_10.xml', 'Anew_tree_3.xml']

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