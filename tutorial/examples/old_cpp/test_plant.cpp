#include "../../../src/structural/MappedOrganism.h"
#include <iostream>
#include <string> 

int main(int argc, char **argv) 
{
	std::string path = "../../../modelparameter/structural/plant/";
    std::string name = "Heliantus_Pagès_2013";
    auto rs = std::make_shared<CPlantBox::Plant>(2);//pb.RootSystem()  # the original
	rs->openXML(path + name + ".xml");
	
	rs->initializeDB(true);
	rs->simulate(76,true);//db
	std::cout <<path + name + ".xml"<<std::endl;
    return 0;
}