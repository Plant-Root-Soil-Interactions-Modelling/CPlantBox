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
    with open("../installdumux.log", "a") as log:
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
        if return_code:
            print("\n")
            message = "\n    (Error) The command {} returned with non-zero exit code\n".format(command)
            message += "\n    If you can't fix the problem yourself consider reporting your issue\n"
            message += "    on the mailing list (dumux@listserv.uni-stuttgart.de) and attach the file 'installdumux.log'\n"
            show_message(message)
            sys.exit(1)


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
programs = ['wget', 'git', 'gcc', 'g++', 'cmake', 'pkg-config','clang', 'gfortran','python3'] 
show_message("(1/3) (a) Checking ubuntu prerequistes: " + " ".join(programs) + "...")

# check some prerequistes
error = []
for program in programs:
    try:
        output2 = subprocess.run([program, "--version"], capture_output=True)
        if program == 'python3':
            start_str = len("b'Python ")
            end_str = start_str+3
            pythonVersion = float(str(output2.stdout)[start_str:end_str])
            if pythonVersion < 3.7:
                error.append(program)
                print("update python to at least version 3.7")
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
	subprocess.run(["pip3", "install", mymodule]) 
      
show_message("(1/3) Step completed. All prerequistes found.")


#################################################################
#################################################################
## (2/3) Clone modules
#################################################################
#################################################################
# make a new folder containing everything
os.makedirs("./DUMUX", exist_ok=True)
os.chdir("DUMUX")

show_message("(2/3) Cloning repositories. This may take a while. Make sure to be connected to the internet...")

dune_version=2.6
dumux_version=3.0
# the core modules
for module in ['common', 'geometry', 'grid', 'localfunctions', 'istl' ]:
    if not os.path.exists("dune-{}".format(module)):
        git_clone('https://gitlab.dune-project.org/core/dune-{}.git'.format(module), "releases/{}".format(dune_version))
    else:
        print("-- Skip cloning dune-{} because the folder already exists.".format(module))
        
'''
#need a mpicommunication.hh file in dune-common for dune-grid-glue >=2.7:
https://dune-project.org/releases/2.7.1/
Deprecated header dune/common/parallel/collectivecommunication.hh which will be removed after Dune 2.7. 
Use dune/common/parallel/communication.hh instead!
Deprecated header dune/common/parallel/mpicollectivecommunication.hh which will be removed after Dune 2.7. 
Use dune/common/parallel/mpicommunication.hh instead!
'''

shutil.copyfile("dune-common/dune/common/parallel/collectivecommunication.hh", "dune-common/dune/common/parallel/communication.hh")
shutil.copyfile("dune-common/dune/common/parallel/mpicollectivecommunication.hh", "dune-common/dune/common/parallel/mpicommunication.hh")

# the extensions
for module in ['foamgrid','grid-glue', 'alugrid', 'spgrid' ]:
    if not os.path.exists("dune-{}".format(module)):
        if module == 'grid-glue' :
            git_clone('https://gitlab.dune-project.org/extensions/dune-{}.git'.format(module), "releases/{}".format(2.7)) #no 2.6 release
        else:
            git_clone('https://gitlab.dune-project.org/extensions/dune-{}.git'.format(module), "releases/{}".format(dune_version))
    else:
        print("-- Skip cloning dune-{} because the folder already exists.".format(module))
    
# staging    
for module in ['uggrid' ]:
    if not os.path.exists("dune-{}".format(module)):
        git_clone('https://gitlab.dune-project.org/staging/dune-{}.git'.format(module), "releases/{}".format(dune_version))
    else:
        print("-- Skip cloning dune-{} because the folder already exists.".format(module))

# pybind
for module in ['pybindxi' ]:
    if not os.path.exists("dune-{}".format(module)):
        git_clone('https://github.com/dune-community/dune-{}.git'.format(module), branch = 'master')
    else:
        print("-- Skip cloning dune-{} because the folder already exists.".format(module))

# dumux and course
if not os.path.exists("dumux"):
    git_clone('https://git.iws.uni-stuttgart.de/dumux-repositories/dumux.git', "releases/{}".format(dumux_version))
else:
    print("-- Skip cloning dumux because the folder already exists.")

#if not os.path.exists("dumux-course"):
 #   git_clone('https://git.iws.uni-stuttgart.de/dumux-repositories/dumux-course.git', "releases/{}".format(dumux_version))
#else:
 #   print("-- Skip cloning dumux-course because the folder already exists.")


# dumux-rosi
if not os.path.exists("dumux-rosi"):
    git_clone('https://github.com/Plant-Root-Soil-Interactions-Modelling/dumux-rosi.git', branch = 'master')
else:
    print("-- Skip cloning dumux-rosi because the folder already exists.")

# CPlantBox
if not os.path.exists("CPlantBox"):
    git_clone('https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox.git', branch = 'master')
    os.chdir("CPlantBox")
    subprocess.run(['cmake', '.']) 
    subprocess.run(['make']) 
    os.chdir("..")
else:
    print("-- Skip cloning CPlantBox because the folder already exists.")


show_message("(2/3) Step completed. All repositories have been cloned into a containing folder.")

#################################################################
#################################################################
## (3/3) Configure and build
#################################################################
#################################################################
show_message("(3/3) Configure and build dune modules and dumux using dunecontrol. This may take several minutes...")


#copy cplantbox so file to dumux-rosi folder
for f_name in os.listdir("CPlantBox"):
    if f_name.startswith('plantbox.cpython') and f_name.endswith('.so'):
        shutil.copyfile("CPlantBox/{0}".format(f_name), "dumux-rosi/{0}".format(f_name))
    

# run dunecontrol
if not os.path.isfile("cmake.opts"):
    shutil.move("dumux/cmake.opts", "cmake.opts")
else:
    print("-- The file cmake.opts already exists. The existing file will be used to configure dumux.")


subprocess.check_output(["./dune-common/bin/dunecontrol", "--opts=cmake.opts", "all"])

show_message("(3/3) Step completed. Succesfully configured and built CPlantBox, dune and dumux.")


show_message("to test installation, run n\ cd DUMUX/dumux-rosi/python/coupled \n python3 example7b_coupling.py")
 
