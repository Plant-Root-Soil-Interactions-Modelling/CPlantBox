import os
import sys
from CPlantBox_PiafMunch import *
import sys; sys.path.append("../../.."); sys.path.append("../../../src/python_modules")
name = "morning_glory_7.xml"#"PMA2018.xml" # parameter name
# here are some optional parameter files to be tested
# name = "PMA2018" # Simulate a small plant with 3 leaves and two lateral root, you can comment the heliantus line and uncomment this line to see what happend.
time = 50# how many days the plant need to grow, make it smaller, for example 15 to see if the plant becomes smaller
plant1 = CPlantBox_PiafMunch(name, time, name) # make a plant object in python
# Visualization
