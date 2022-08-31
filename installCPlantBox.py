#!/usr/bin/env python3

"""
One click install script for dumux-rosi
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


def run_command(command):
    with open("../installCPlantBox.log", "a") as log:
        popen = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        for line in popen.stdout:
            log.write(line)
            print(line, end='')
        for line in popen.stderr:
            log.write(line)
            print(line, end='')
        popen.stdout.close()
        popen.stderr.close()
        return_code = popen.wait()


def git_clone(url, branch=None):
    clone = ["git", "clone"]
    if branch:
        clone += ["-b", branch]
    result = run_command(command=[*clone, url])


# clear the log file
open('installdumux.log', 'w').close()

show_message("do not forget to run \n sudo apt update \n sudo apt upgrade \n\n only works for ubuntu >= 20.04")

#################################################################
#################################################################
## (1/3) Check some prerequistes
#################################################################
#################################################################
programs = ['wget', 'git', 'gcc', 'g++', 'cmake', 'pkg-config','clang', 'gfortran']#,'python3'] 
show_message("(1/3) (a) Checking ubuntu prerequistes: " + " ".join(programs) + "...")

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

   
programs = ['default-jre', 'libboost-all-dev', 'python3-pip'] 

for program in programs:
    output = subprocess.run(["dpkg", "-l", program], capture_output=True)
    if ('no packages found' in str(output)):
        error.append(program)
        
if len(error) > 0:
    print("Program(s) {0} has/have not been found. try running sudo apt-get install {0}".format(" ".join(error)))
    raise Exception('import modules')

import pip

# check some prerequistes
modules = ['numpy', 'scipy', 'matplotlib', 'vtk', 'mpi4py', 'astropy', 'pandas'] 
show_message("(1/3) (b) Checking python prerequistes: " + " ".join(modules) + "...")

for mymodule in modules:
	#subprocess.run(["pip3", "install", mymodule]) 
	subprocess.run(["conda", "install", mymodule]) 
      
show_message("(1/3) Step completed. All prerequistes found.")


#################################################################
#################################################################
## (2/3) Clone modules
#################################################################
#################################################################
# make a new folder containing everything
os.makedirs("./CPB", exist_ok=True)
os.chdir("CPB")

# CPlantBox
if not os.path.exists("CPlantBox"):
    git_clone('https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox.git', branch = 'master')
    os.chdir("CPlantBox")
    subprocess.run(['cmake', '.']) 
    subprocess.run(['make']) 
    os.chdir("..")
else:
    print("-- Skip cloning CPlantBox because the folder already exists.")
    os.chdir("CPlantBox")
    subprocess.run(['cmake', '.']) 
    subprocess.run(['make']) 
    os.chdir("..")


show_message("(3/3) Step completed. Succesfully configured and built CPlantBox, dune and dumux.")


show_message("to test installation, run n\ cd ~/CPB/CPlantBox/tutorial/examples/python \n python3 example1a_small.py")