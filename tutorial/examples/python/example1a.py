"""small example"""
import py_plantbox as rb

# print("hallo")

v = rb.Vector3d(1, 2, 3)
print(v)

# sdf = rb.SignedDistanceFunction();
# print(sdf)
# print("hello")

# import py_rootbox as rb
#
# rootsystem = rb.RootSystem()
#
# # Open plant and root parameter from a file
# name = "Zea_mays_1_Leitner_2010"  # "Anagallis_femina_Leitner_2010"
# rootsystem.readParameters("modelparameter/" + name + ".xml")
#
# # Initialize
# rootsystem.initialize()
#
# # Simulate
# rootsystem.simulate(30, True)
#
# # Export final result (as vtp)
# rootsystem.write("../results/example_1a.vtp")
#
# print("done.")
