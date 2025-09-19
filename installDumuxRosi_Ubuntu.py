#!/usr/bin/env python3

"""
One click install script for dumux-rosi
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
show_message("We recommend you to setup CPlantBox in a virtual environment by running\n\npython3 -m venv cpbenv\nsource cpbenv/bin/activate")

#################################################################
#################################################################
# # (1/3) Check some prerequistes
#################################################################
#################################################################
programs = ['wget', 'git', 'gcc', 'g++', 'cmake', 'pkg-config', 'clang', 'gfortran']  # ,'python3']
show_message("(1/3) (a) Checking ubuntu prerequistes: " + " ".join(programs) + "...")

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
show_message("(1/3) (b) Checking python prerequistes: " + " ".join(modules) + "...")

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
subprocess.run(["python3", "-m", "pip","install", "--no-cache-dir", "mpi4py", "--verbose"])  # mpi4py can take a lot of time to install

show_message("(1/3) Step completed. All prerequistes found.")

#################################################################
#################################################################
# # (2/3) Clone modules
#################################################################
#################################################################

url = 'https://raw.githubusercontent.com/Plant-Root-Soil-Interactions-Modelling/dumux-rosi/master/installdumux.py'
subprocess.run(["wget", url, "-P", "."], check = True)

subprocess.run(["python3", "installdumux.py"], check = True)

os.chdir("./dumux")

subprocess.run(["python3", "dumux/bin/installexternal.py", "ug", "spgrid", "foamgrid"], check = True)

# dumux-rosi
if not os.path.exists("dumux-rosi"):
    subprocess.run(['git', 'clone', '--depth', '1', '-b', 'master', 'https://github.com/Plant-Root-Soil-Interactions-Modelling/dumux-rosi.git'])
else:
    print("-- Skip cloning dumux-rosi because the folder already exists.")

# CPlantBox
if not os.path.exists("CPlantBox"):
    subprocess.run(['git', 'clone', '--depth', '1', '-b', 'master', 'https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox.git'])
else:
    print("-- Skip cloning CPlantBox because the folder already exists.")

#################################################################
#################################################################
# # (3/3) Configure and build
#################################################################
#################################################################
show_message("(3/3) Configure and build dune modules and dumux using dunecontrol. This may take several minutes...")

os.chdir("CPlantBox")

subprocess.run(['git', 'submodule', 'update', '--recursive', '--init'])
subprocess.run(['cmake', '.'])
subprocess.run(['make'])
os.chdir("..")

# run dunecontrol
subprocess.run(["./dune-common/bin/dunecontrol", "--opts=dumux-rosi/cmake.opts", "all"])

print("(3/3) Step completed. Succesfully configured and built CPlantBox, dune and dumux.")

print("to test installation, run \n cd dumux/dumux-rosi/python/coupled \n python3 example7b_coupling.py")
