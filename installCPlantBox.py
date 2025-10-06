#!/usr/bin/env python3

"""
One click install script for CPlantBox
adapted from installdumux.py of the Dumux developers
"""
import os
import sys
import shutil
import subprocess
import fileinput


def show_message(message):
    print("*" * 120)
    print(message)
    print("*" * 120)


def run_command(command):
    with open("../installDumuxRosi.log", "a") as log:
        popen = subprocess.Popen(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE, universal_newlines = True)
        for line in popen.stdout:
            log.write(line)
            print(line, end = '')
        for line in popen.stderr:
            log.write(line)
            print(line, end = '')
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


def git_clone(url, branch = None):
    clone = ["git", "clone"]
    if branch:
        clone += ["-b", branch]
    result = run_command(command = [*clone, url])


# clear the log file
open('installDumuxRosi.log', 'w').close()

show_message("do not forget to run \nsudo apt update \nsudo apt upgrade")
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
isCluster = ('ENV' in os.environ.keys())

programs = ['default-jre', 'libeigen3-dev' , 'python3-pip', 'openmpi-bin', 'libopenmpi-dev', 'python3-tk', 'libqt5x11extras5', 'libx11-dev']  #
# sudo apt install openmpi-bin libopenmpi-dev
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

#################################################################
#################################################################
# # (2/2) Clone modules
#################################################################
#################################################################
# make a new folder containing everything
# os.makedirs("./CPB", exist_ok=True)
# os.chdir("CPB")

# CPlantBox
if not os.path.exists("CPlantBox"):
    subprocess.run(['git', 'clone', '--depth', '1', '-b', 'master', 'https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox.git'])
else:
    print("-- Skip cloning CPlantBox because the folder already exists.")
os.chdir("CPlantBox")

subprocess.run(['git', 'submodule', 'update', '--recursive', '--init'])
subprocess.run(['cmake', '.'])
subprocess.run(['make'])
os.chdir("..")

show_message("(2/2) Step completed. Succesfully configured and built CPlantBox.")

show_message("to test installation, run \n cd CPlantBox/tutorial/examples/ \n python3 example1a_small.py")

show_message("CPlantBox is currently at stable branch, use \n $git switch master \n to obtain the latest version, use cmake . & make to recompile")


