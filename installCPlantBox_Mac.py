#!/usr/bin/env python3

"""
One click install script for CPlantBox (macOS version)
"""

import os
import sys
import subprocess


def show_message(message):
    print("*" * 120)
    print(message)
    print("*" * 120)


def run_command(command):
    with open("../installDumuxRosi.log", "a") as log:
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
            show_message(f"(Error) Command failed: {' '.join(command)}")
            sys.exit(1)


def git_clone(url, branch=None):
    clone = ["git", "clone"]
    if branch:
        clone += ["-b", branch]
    run_command([*clone, url])


# clear log
open('installCPlantBox.log', 'w').close()

show_message("macOS setup:\tMake sure Homebrew is installed\n\tMake sure python3 >= 3.10, < 3.14")
show_message("Recommended: use a virtual environment:\npython3 -m venv cpbenv\nsource cpbenv/bin/activate")

#################################################################
# (1/3) Check prerequisites
#################################################################

programs = ['git', 'cmake', 'pkg-config', 'clang', 'gfortran']
show_message("(1/3) Checking required tools: " + " ".join(programs))

error = []
for program in programs:
    try:
        subprocess.run([program, "--version"], capture_output=True)
    except FileNotFoundError:
        error.append(program)

if error:
    print("\nMissing tools:", error)
    print("Install with Homebrew:")
    print(f"brew install {' '.join(error)}")
    sys.exit(1)

# Check brew packages
brew_packages = [
    'cmake',
    'pkg-config',
    'gcc',
    'wget',
    'eigen',
    'boost',
    'open-mpi',
    'python',
    'qt'
]

missing = []
for pkg in brew_packages:
    result = subprocess.run(["brew", "list", pkg], capture_output=True)
    if result.returncode != 0:
        missing.append(pkg)

if missing:
    print("\nMissing Homebrew packages:", missing)
    print("Run:")
    print(f"brew install {' '.join(missing)}")
    sys.exit(1)

#################################################################
# Python dependencies
#################################################################

modules = ['numpy', 'scipy', 'matplotlib', 'vtk', 'pandas', 'pybind11', 'ipython']
show_message("(1/3) Installing Python modules: " + " ".join(modules))

for m in modules:
    subprocess.run(["python3", "-m", "pip", "install", m])

subprocess.run(["python3", "-m", "pip", "install", "--no-cache-dir", "mpi4py"])

show_message("(1/3) Step completed.")

#################################################################
# (2/3) Clone modules
#################################################################

# CPlantBox
if not os.path.exists("CPlantBox"):
    subprocess.run(['git', 'clone', '--depth', '1', '-b', 'master',
                    'https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox.git'])

#################################################################
# (3/3) Build
#################################################################

show_message("(3/3) Building modules...")

os.chdir("CPlantBox")

subprocess.run(['git', 'submodule', 'update', '--recursive', '--init'])
subprocess.run(['cmake', '-DPIAFMUNCH=OFF','.'])
subprocess.run(['make', 'install'])


if "VIRTUAL_ENV" in os.environ:
    pyver = f"{sys.version_info.major}.{sys.version_info.minor}"
    plantbox_path = os.path.join(
        os.environ["VIRTUAL_ENV"],
        f"lib/python{pyver}/site-packages/plantbox"
    )
    python_path = os.path.join(
        os.environ["VIRTUAL_ENV"],
        f"lib/python{pyver}/site-packages"
    )
else:
    raise Exception('virtual environment not activated')



from pathlib import Path

zshrc_path = Path.home() / ".zshrc"

lines = [f'export DYLD_LIBRARY_PATH="{plantbox_path}:$DYLD_LIBRARY_PATH"\n',
	f'export PYTHONPATH="{python_path}:$PYTHONPATH"\n'
		]
for line in lines:
    # Append only if it's not already there
    if zshrc_path.exists():
        content = zshrc_path.read_text()
    else:
        content = ""

    if line.strip() not in content:
        with open(zshrc_path, "a") as f:
            f.write("\n# Added by my Python script\n")
            f.write(line)
            
final_message = """Setup finished. Implement the changes made to your kernel by running\nsource ~/.zshrc
Test it by running\ncd CPlantBox/tutorial/chapter1_introduction\npython3 example1_3_helloplant.py
ATT: when running for the first time, 'import vtk' may take a long time."""
show_message(final_message)
