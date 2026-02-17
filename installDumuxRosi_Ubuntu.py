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

#################################################################
#################################################################
# # (2/3) Clone modules
#################################################################
#################################################
# CPlantBox
if not os.path.exists("CPlantBox"):
    subprocess.run(['git', 'clone', '--depth', '1', '-b', 'uclouvain3', 'https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox.git'])
else:
    print("-- Skip cloning CPlantBox because the folder already exists.")

#################################################################
#################################################################
# # (3/3) Configure and build
#################################################################
#################################################################
show_message("(3/3) Configure and build cpb. This may take several minutes...")

os.chdir("CPlantBox")

subprocess.run(['git', 'submodule', 'update', '--recursive', '--init'])
subprocess.run(['cmake', '.'])
subprocess.run(['make'])
os.chdir("..")

print("(3/3) Step completed. Succesfully configured and built CPlantBox.")

print("to test installation, run \n cd dumux/dumux-rosi/python/coupled \n python3 example7b_coupling.py")


