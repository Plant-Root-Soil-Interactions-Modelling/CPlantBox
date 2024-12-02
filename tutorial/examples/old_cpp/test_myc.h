#include <vector>
#include <string>
#include <iostream>


#include "MycorrhizalPlant.h"
#include "MycorrhizalPlant.cpp"
#include "MycorrhizalRoot.h"
#include "MycorrhizalRoot.cpp"
#include "Mycorrhizalrootparameter.h"
#include "Mycorrhizalrootparameter.cpp"

int testmyc() {
    auto mycplant = std::make_shared<CPlantBox::MycorrhizalPlant>();
    mycplant -> readParameters("fspm2023.xml");
    mycplant ->initialize();

    auto mycroot = std::make_shared<CPlantBox::MycorrhizalRoot>();
    return 0;
}