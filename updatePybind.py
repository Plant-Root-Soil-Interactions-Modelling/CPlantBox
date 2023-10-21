#!/usr/bin/env python3

"""
One click install script for CPlantBox
adapted from installdumux.py of the Dumux developers
"""
import os
import sys
import shutil
import subprocess


def show_message(message):
    print("*" * 120) 
    print(message)
    print("*" * 120)




if os.path.exists("./src/external/pybind11"):
    subprocess.run(['rm', '-rf', 'src/external/pybind11'])#delete folder
subprocess.run(['git', 'rm', '-r','--cached', 'src/external/pybind11'])#take out git cache for pybind11
subprocess.run(['git', 'submodule', 'add',  '--force', '-b', 'stable', 'https://github.com/pybind/pybind11.git', './src/external/pybind11'])
subprocess.run(['cmake', '.']) 
subprocess.run(['make'])  
os.chdir("..")


show_message("(3/3) Step completed. Succesfully configured and built CPlantBox.")


show_message("to test installation, run n\ cd CPlantBox/tutorial/examples/ \n python3 example1a_small.py")
