#!/usr/bin/env python3

""" lkiugvofo
One click install script for CPlantBox
adapted from installdumux.py of the Dumux developers
"""
import os
import sys
import subprocess


def show_message(message):
    print("*" * 120)
    print(message)
    print("*" * 120)


# clear the log file
open('installCPlantBox.log', 'w').close()

show_message("Do not forget to run \nsudo apt update \nsudo apt upgrade")
show_message("We recommend you to setup CPlantBox in a virtual environment by running\npython3 -m venv cpbenv\nsource cpbenv/bin/activate")

#################################################################
#################################################################
# # (1/2) Check some prerequistes
#################################################################
#################################################################
programs = ['wget', 'git', 'gcc', 'g++', 'cmake', 'pkg-config', 'clang', 'gfortran']  # ,'python3']
show_message("(1/2) (a) Checking ubuntu prerequistes: " + " ".join(programs) + "...")

# check some prerequistes
if (sys.version_info.major < 3) or (sys.version_info.minor < 7):
    print("detected python version:", sys.version_info)
    raise Exception("update python to at least version 3.7")

error = []
for program in programs:
    try:
        output2 = subprocess.run([program, "--version"], capture_output = True)
    except FileNotFoundError:
        error.append(program)

# is the script running on (agro)cluster?
# tried to make evaluation automatic but not sure it holds on all machines
isCluster = ('MODULESHOME' in os.environ.keys())

programs = ['default-jre', 'libeigen3-dev' , 'python3-pip', 'openmpi-bin', 'libopenmpi-dev', 'python3-tk', 'libqt5x11extras5', 'libx11-dev']

if not isCluster:
    programs.append('libboost-all-dev')

    for program in programs:
        output = subprocess.run(["dpkg", "-l", program], capture_output = True)
        if ('no packages found' in str(output)):
            error.append(program)

    if len(error) > 0:
        print("Program(s) {0} has/have not been found. try running sudo apt-get install {0}".format(" ".join(error)))
        raise Exception('import modules')
        
# check some prerequistes
modules = ['numpy', 'scipy', 'matplotlib', 'vtk','cmake', 'pandas', 'pybind11[global]', 'ipython']  # 'mpi4py',
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

subprocess.run(["python3", "-m", "pip", "install", "--no-cache-dir", "mpi4py", "--verbose"])  # mpi4py can take a lot of time to install

show_message("(1/2) Step completed. All prerequistes found.")

#################################################################
#################################################################
# # (2/2) Clone CPlantBox
#################################################################
#################################################################

if not os.path.exists("CPlantBox"):
    if len(sys.argv) < 2:
        subprocess.run(['git', 'clone', '--depth', '1', '-b', 'master', 'https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox.git'])
    else:
        if sys.argv[1] == "full":
            subprocess.run(['git', 'clone', '-b', 'master', 'https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox.git'])
        else:
            show_message(f"Warning: Unknown argument: {sys.argv[1]}")
            subprocess.run(['git', 'clone', '--depth', '1', '-b', 'master', 'https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox.git'])
else:
    print("-- Skip cloning CPlantBox because the folder already exists.")

os.chdir("CPlantBox")
subprocess.run(['git', 'submodule', 'update', '--recursive', '--init'])
subprocess.run(['cmake', '.'])
subprocess.run(['make', 'install'])
os.chdir("..")

show_message("(2/2) Step completed. Succesfully configured and built CPlantBox.")

show_message("To test installation, run \n cd CPlantBox/tutorial/examples/ \n python3 example1a_small.py")

show_message("CPlantBox was installed in your python environment (master branch), use 'cmake . & make install' to recompile")




