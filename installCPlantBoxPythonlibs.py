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

programs = ['wget', 'git', 'gcc', 'g++', 'cmake', 'pkg-config','clang', 'gfortran', 'pip']#,'python3'] 
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
modules = ['numpy', 'scipy', 'matplotlib', 'vtk', 'pandas', 'pybind11[global]', 'ipython']  # 'mpi4py',
show_message("(1/2) (b) Checking python prerequistes: " + " ".join(modules) + "...")

for mymodule in modules:
    try:
        subprocess.run(["python3", "-m", "pip", "install", mymodule], check = True, text = True, capture_output = True)  # ,"-y"
    except subprocess.CalledProcessError as e:
        # Check for the "externally-managed-environment" error
        if "externally-managed-environment" in e.stderr:
            print(f"\n\n\n\n\033[31mError\033[0m: {mymodule} cannot be installed in an externally managed environment.")
            print("run:")
            print("\033[35m[ ! -d 'cpbenv' ] && python3 -m venv cpbenv &&  source cpbenv/bin/activate ||  source cpbenv/bin/activate\033[0m")
            print("\n\n\n\n")
        else:
            print(f"An error occurred: {e.stderr}")
        raise Exception

# pip3 install --upgrade pip setuptools wheel # necessary?
subprocess.run(["python3", "-m", "pip", "install", "--no-cache-dir", "mpi4py", "--verbose"])  # mpi4py can take a lot of time to install

show_message("(1/2) Step completed. All prerequistes found.")