 examples/Makefile                                   |   30 [32m+[m[31m-[m
 examples/benchmarks.h                               |  228 [32m+++++[m
 examples/cmake_install.cmake                        |    2 [32m+[m[31m-[m
 examples/example1.h                                 |  236 [32m++++[m[31m-[m
 examples/example2.h                                 |  234 [32m+++++[m
 examples/example3.h                                 |  159 [32m+++[m
 examples/example4.h                                 |  288 [32m++++++[m
 examples/example5.h                                 |  214 [32m++++[m
 examples/example_scenegraph.h                       |  225 [32m+++++[m
 examples/main.cpp                                   |   11 [32m+[m[31m-[m
 examples/modelparameter                             |    2 [32m+[m[31m-[m
 examples/shehan_RhizoTubes.h                        |  603 [32m+++++++++++[m
 examples/shehan_SoilCore.h                          |    5 [32m+[m
 examples/shehan_Trenches.h                          |  573 [32m+++++++++++[m
 examples/test_cplantbox                             |  Bin [31m628480[m -> [32m652488[m bytes
 modelparameter/AMAIZE_new.xml                       |  244 [31m-----[m
 modelparameter/morning_glory_7m.xml                 |   67 [32m+[m[31m-[m
 modelparameter/morning_glory_7m_new.xml             |   68 [32m+[m[31m-[m
 modelparameter/sympodial_monochasium.xml            |  411 [32m+[m[31m-------[m
 pyscript/tubePlot.py                                |   64 [32m+[m[31m-[m
 python/CPlantBox.ipynb                              |    2 [32m+[m[31m-[m
 python/Makefile                                     |   20 [32m+[m[31m-[m
 python/PiafMunch2_python_test_mg.ini                | 3251 [32m++++++++++++++++++++++++++++++++++++++++++++[m[31m---------------[m
 python/Test PMA1.ipynb                              | 2013 [32m++++++++++++++++++++++++++[m[31m----------[m
 python/Untitled1.ipynb                              |    6 [32m+[m[31m-[m
 python/__pycache__/rb_tools.cpython-36.pyc          |  Bin [31m3032[m -> [32m3012[m bytes
 python/cmake_install.cmake                          |    2 [32m+[m[31m-[m
 python/ee.py                                        |    6 [32m+[m[31m-[m
 python/modelparameter                               |    2 [32m+[m[31m-[m
 python/morning glory.ipynb                          | 1702 [32m++++++++[m[31m-----------------------[m
 python/morning glory_defolie 2.5m.ipynb             |    2 [32m+[m[31m-[m
 python/py_plantbox.so                               |  Bin [31m3571080[m -> [32m3718936[m bytes
 src/CMakeFiles/CPlantBox.dir/Leaf.cpp.o             |  Bin [31m57992[m -> [32m58032[m bytes
 src/CMakeFiles/CPlantBox.dir/LeafTropism.cpp.o      |  Bin [31m33736[m -> [32m32288[m bytes
 src/CMakeFiles/CPlantBox.dir/ModelParameter.cpp.o   |  Bin [31m249448[m -> [32m245832[m bytes
 src/CMakeFiles/CPlantBox.dir/Root.cpp.o             |  Bin [31m55496[m -> [32m53736[m bytes
 src/CMakeFiles/CPlantBox.dir/RootTropism.cpp.o      |  Bin [31m40064[m -> [32m38608[m bytes
 src/CMakeFiles/CPlantBox.dir/Stem.cpp.o             |  Bin [31m64032[m -> [32m64040[m bytes
 src/CMakeFiles/CPlantBox.dir/StemTropism.cpp.o      |  Bin [31m33736[m -> [32m32288[m bytes
 src/CMakeFiles/py_plantbox.dir/Leaf.cpp.o           |  Bin [31m53584[m -> [32m53624[m bytes
 src/CMakeFiles/py_plantbox.dir/LeafTropism.cpp.o    |  Bin [31m32752[m -> [32m32088[m bytes
 src/CMakeFiles/py_plantbox.dir/ModelParameter.cpp.o |  Bin [31m239312[m -> [32m238720[m bytes
 src/CMakeFiles/py_plantbox.dir/PythonPlant.cpp.o    |  Bin [31m4743864[m -> [32m4741832[m bytes
 src/CMakeFiles/py_plantbox.dir/Root.cpp.o           |  Bin [31m54016[m -> [32m51608[m bytes
 src/CMakeFiles/py_plantbox.dir/RootTropism.cpp.o    |  Bin [31m39048[m -> [32m38384[m bytes
 src/CMakeFiles/py_plantbox.dir/Stem.cpp.o           |  Bin [31m57464[m -> [32m57504[m bytes
 src/CMakeFiles/py_plantbox.dir/StemTropism.cpp.o    |  Bin [31m32752[m -> [32m32088[m bytes
 src/Leaf.cpp                                        |    4 [32m+[m[31m-[m
 src/LeafTropism.cpp                                 |   17 [32m+[m[31m-[m
 src/LeafTropism.h                                   |   18 [32m+[m[31m-[m
 src/Root.cpp                                        |   71 [32m+[m[31m-[m
 src/Root.h                                          |    1 [32m+[m
 src/RootTropism.cpp                                 |   13 [32m+[m[31m-[m
 src/RootTropism.h                                   |   12 [32m+[m[31m-[m
 src/Stem.cpp                                        |   26 [32m+[m[31m-[m
 src/StemTropism.cpp                                 |   17 [32m+[m[31m-[m
 src/StemTropism.h                                   |   18 [32m+[m[31m-[m
 src/libCPlantBox.a                                  |  Bin [31m1165198[m -> [32m1155532[m bytes
 src/mymath.h                                        |    5 [32m+[m[31m-[m
 59 files changed, 7347 insertions(+), 3525 deletions(-)
