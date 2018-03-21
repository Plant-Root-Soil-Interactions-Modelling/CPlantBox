# CPlantBox

The fastest way to try CPlantBox is to read the examples. Just uncomment the example in the main.cpp file to try, compile and run it. 

The code should compile with any c++11 compiler, e.g. for g++:
$ g++ *.cpp -std=c++11
$ ./a.out


# Folder sructure

/				CPlantBox C++ codes

# Documentation

Create the documentation by running doxygen in the folder 
$ doxygen doxy_config

The documentation should now be located in the folder /doc

To build the shared library py_rootbox for coupling with Python pleaser refer to 'python building guide.txt'

# Example

Visualized in R, thanks  [@guillaumelobet](https://github.com/guillaumelobet)

![alt text](https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox/blob/master/results/plant.gif "Tree with leafs")


# Boost.Python Binding Guide in windows
Requirements
Microsoft Visual Studio 2017 (Community Edition)
boost_1_66_0
Python 3.6
 
1. Install Visual studio
and Install python3 from official website

2. Install  boost
Download boost from http://www.boost.org/users/history/version_1_66_0.html
Unzip to a PATH
After normal install http://www.boost.org/doc/libs/1_66_0/more/getting_started/windows.html
 ( .\bootstrap.bat then .\b2 ), project-config.jam in root of boost to
```
import option ; 
using msvc ; 
option.set keep-going : false ; 
using python : 3.6 ; 
```
Then run in cmd or powershell
```
.\b2 --with-python --toolset=msvc --build-type=complete
```
then test it follow http://www.discoversdk.com/blog/python-cpp-interop

If it is working then do the same to visual studio project of CPlantBox
 
put 
```
#define BOOST_LIB_NAME boost_python3
#define BOOST_PYTHON_STATIC_LIB
```
in PythonPlant.cpp if it is not there. Then run it if you see
```
Creating library C:\Users\xxxx\source\CPlantBox\x64\Debug\CPlantBoxVS.lib and object C:\Users\xxxx\source\CPlantBox\x64\Debug\CPlantBoxVS.exp

1>CPlantBoxVS.vcxproj -> C:\Users\xxxx\source\CPlantBox\x64\Debug\CPlantBoxVS.dll

1>Done building project "CPlantBoxVS.vcxproj".

========== Build: 1 succeeded, 0 failed, 0 up-to-date, 0 skipped ==========
```
Then copy "CPlantBoxVS.dll" to the root of CPlantBox and rename it to "py_plantbox.pyd"

then run the example1a.py to test
