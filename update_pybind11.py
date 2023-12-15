#!/usr/bin/env python3

"""
clone pybind11
"""
import os
import sys
import shutil
import subprocess

if os.path.exists("./src/external/pybind11"):
    subprocess.run(['rm', '-rf', 'src/external/pybind11'])  # delete folder
subprocess.run(['git', 'rm', '-r', '--cached', 'src/external/pybind11'])  # take out git cache for pybind11
subprocess.run(['git', 'submodule', 'add', '--force', '-b', 'stable', 'https://github.com/pybind/pybind11.git', './src/external/pybind11'])

subprocess.run(['cmake', '.'])
subprocess.run(['make'])
os.chdir("..")

