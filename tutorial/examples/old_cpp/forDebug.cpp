#include "../../../src/MappedOrganism.h"
#include <iostream>


int main(int argc, char **argv) 
{
    auto rs = std::make_shared<CPlantBox::MappedPlant>(2);

    std::string path = "../../../modelparameter/plant/";
    std::string name = "manyleaves.xml"; // "maize_p1_zero_std", "maize_p2_zero_std", "maize_p3_zero_std"
    rs->readParameters(path+name);
	double dt =1.;
    rs->initialize(false, false);
	rs->simulate(dt);
    return 0;
}
