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




show_message("do not forget to run \n sudo apt update \n sudo apt upgrade \n\n only works for ubuntu >= 20.04")

#################################################################
#################################################################
## (1/3) Check some prerequistes
#################################################################
#################################################################

programs = ['wget', 'git', 'gcc', 'g++', 'cmake', 'pkg-config','clang', 'gfortran']#,'python3'] 
show_message("(1/3) Checking ubuntu prerequistes: " + " ".join(programs) + "...")

# check some prerequistes
if (sys.version_info.major < 3) or (sys.version_info.minor < 7):
    print("detected pyhton version:",sys.version_info)
    raise Exception("update python to at least version 3.7")
    
error = []
for program in programs:
    try:
        output2 = subprocess.run([program, "--version"], capture_output=True)
    except FileNotFoundError:
        error.append(program)


#is the script running on (agro)cluster?
#tried to make evaluation automatic but not sure it holds on all machines
isCluster = ('ENV' in os.environ.keys())
        
programs = ['default-jre', 'python3-pip','libeigen3-dev'] 
if not isCluster:
    programs.append('libboost-all-dev')
    
for program in programs:
    output = subprocess.run(["dpkg", "-l", program], capture_output=True)
    if ('no packages found' in str(output)):        
        error.append(program)
        
if len(error) > 0:
    print("Program(s) {0} has/have not been found. try running sudo apt-get install {0}".format(" ".join(error)))
    raise Exception('import modules')

import pip

# check some prerequistes
modules = ['numpy', 'scipy', 'matplotlib', 'vtk', 'mpi4py', 'astropy', 'pandas', 'pybind11[global]'] 
show_message("(2/3) Checking python prerequistes: " + " ".join(modules) + "...")

for mymodule in modules:
    #subprocess.run(["pip3", "install", mymodule]) 
    if ((mymodule =='vtk') and (sys.version_info.minor == 10)):
        subprocess.run(["pip3", "install", mymodule]) 	#conda install not working for vtk with py3.10 (?)
    else:
        try:
            subprocess.run(["conda", "install", mymodule]) 
        except:
            subprocess.run(["pip3", "install", mymodule]) 
      
show_message("(2/3) Step completed. All prerequistes found.")


#################################################################
#################################################################
## (3/3) Clone modules
#################################################################
#################################################################
# make a new folder containing everything
#os.makedirs("./CPB", exist_ok=True)
#os.chdir("CPB")

# CPlantBox
if not os.path.exists("CPlantBox"):
    subprocess.run(['git', 'clone', '--depth','1','-b', 'stable_v2.1', 'https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox.git'])
else:
    print("-- Skip cloning CPlantBox because the folder already exists.")
os.chdir("CPlantBox")

if os.path.exists("./src/external/pybind11"):
    subprocess.run(['rm', '-rf', 'src/external/pybind11'])#delete folder
subprocess.run(['git', 'rm', '-r','--cached', 'src/external/pybind11'])#take out git cache for pybind11
subprocess.run(['git', 'submodule', 'add',  '--force', '-b', 'stable', 'https://github.com/pybind/pybind11.git', './src/external/pybind11'])
subprocess.run(['cmake', '.']) 
subprocess.run(['make'])  
os.chdir("..")


show_message("(3/3) Step completed. Succesfully configured and built CPlantBox.")


show_message("to test installation, run n\ cd CPB/CPlantBox/tutorial/examples/ \n python3 example1a_small.py")
