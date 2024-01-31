#include "Plant.h"
#include <stdlib.h>
#include <memory>
#include <iostream>
#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <fstream>

namespace CPlantBox {


void cpb(std::string para, double time)
{
    auto plant = std::make_shared<Plant>();

 
    plant->openXML(para);

    plant->initialize();

    /*
     * Simulate
     */
    
    plant->simulate(time);
    /*
     * Export final result (as vtp)
     */
    plant->write("output.vtp");


}

} // end namespace CRootBox