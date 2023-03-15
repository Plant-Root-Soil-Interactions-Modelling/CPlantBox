#include "../../../src/structural/MappedOrganism.h"
#include <iostream>
#include <string> 

int main(int argc, char **argv) 
{
    auto rs = std::make_shared<CPlantBox::MappedPlant>(2);

    std::string path = "../../../modelparameter/structural/plant/";
    std::string name ="test_relcoord.xml";
	//"manyleaves.xml";//"Triticum_aestivum_adapted_2023";
	// "test_relcoord.xml"; //"manyleaves.xml";
	rs->readParameters(path+name);
	double dt =1.;
    rs->initialize();//false, false);
	int steps = 100;
	for(int i=0;i<steps;i++)
	{	

		rs->simulate(dt, true);
		rs->write("test"+std::to_string(i)+".vtp");
		rs->plant()->write("testPlant"+std::to_string(i)+".vtp");
	}
	
    return 0;
}