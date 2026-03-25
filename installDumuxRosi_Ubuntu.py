#!/usr/bin/env python

replace dunecontrol by the one in dumux-rosi
in \dune-spgrid\dune\grid\spgrid
replace partitionlist and linkage by the one in dumux-rosi
install MSYS2
and install MSYS MinGW inside.

cd ..

Right-click MSYS2 MinGW64 Shell => Run as Administrator
Then:

pacman -S mingw-w64-x86_64-tinyxml2
pacman -S mingw-w64-x86_64-eigen3
pacman -S mingw-w64-x86_64-pybind11

cd C:/Users/m.giraud/Documents/CPlantBox-windows
source CPlantBox-windows/.venv/Scripts/activate
cd CPlantBox
cmake -S . -B build \
  -DBUILD_CPLANTBOX_SHARED=OFF \
  -DPython_EXECUTABLE=C:/Users/m.giraud/Documents/CPlantBox-windows/CPlantBox-windows/Python312/Python312/python.exe \
  -DPython_ROOT_DIR=C:/Users/m.giraud/Documents/CPlantBox-windows/CPlantBox-windows/Python312/Python312 \
  -DPython_FIND_STRATEGY=LOCATION
cmake --build build --config Release


# not needed anymore?
# in powershell?:
# cd vcpkg 
# ./vcpkg install pkgconf:x64-windows

in git bash: <== to redo everytime
export PATH="/c/Users/m.giraud/Documents/CPlantBox-windows/vcpkg/installed/x64-windows/tools/pkgconf:$PATH"
export PATH="/c/Program Files/Git/cmd:$PATH"

# not needed anymore?
# ln -s pkgconf.exe \
#/c/Users/m.giraud/Documents/CPlantBox-windows/vcpkg/installed/x64-windows/tools/pkgconf/pkg-config.exe


git config --global --add safe.directory C:/Users/m.giraud/Documents/CPlantBox-windows/dune-common
git config --global --add safe.directory C:/Users/m.giraud/Documents/CPlantBox-windows/dune-geometry
git config --global --add safe.directory C:/Users/m.giraud/Documents/CPlantBox-windows/dune-foamgrid
git config --global --add safe.directory C:/Users/m.giraud/Documents/CPlantBox-windows/dune-grid
git config --global --add safe.directory C:/Users/m.giraud/Documents/CPlantBox-windows/dune-istl
git config --global --add safe.directory C:/Users/m.giraud/Documents/CPlantBox-windows/dune-localfunctions
git config --global --add safe.directory C:/Users/m.giraud/Documents/CPlantBox-windows/dune-spgrid
git config --global --add safe.directory C:/Users/m.giraud/Documents/CPlantBox-windows/dune-uggrid
git config --global --add safe.directory C:/Users/m.giraud/Documents/CPlantBox-windows/dune-alugrid
git config --global --add safe.directory C:/Users/m.giraud/Documents/CPlantBox-windows/dumux


# In Git Bash  <== to redo everytime?
export MSYS2_ARG_CONV_EXCL="*"
export MSYS_NO_PATHCONV=1

put posix_compat in the root dumux folder

if needed:
 ./dune-common/bin/dunecontrol bexec rm -r CMakeFiles CMakeCache.txt

./dune-common/bin/dunecontrol --verbose --opts=dumux-rosi/cmake_nompi.opts all -j $(nproc)

to just re-build specific stuff:
./dune-common/bin/dunecontrol --opts=dumux-rosi/cmake_nompi.opts --only=dumux-rosi all -j $(nproc)



