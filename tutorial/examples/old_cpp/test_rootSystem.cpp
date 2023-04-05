#include "../../../src/structural/RootSystem.h"
#include <iostream>
#include <string> 

int main(int argc, char **argv) 
{
	int seed = 110 ;// random seed
    std::string name = "Brassica_oleracea_Vansteenkiste_2014";
    auto rs = std::make_shared<CPlantBox::RootSystem>();//pb.RootSystem()  # the original
	rs->readParameters("../../../modelparameter/structural/rootsystem/" + name + ".xml");
	rs->setSeed(seed);
	rs->initialize(false);
	auto rs2 = rs->copy();// copy root system
	auto allORP = rs->getOrganRandomParameter(2) ;
	std::cout<<"num getOrganRandomParameter "<<allORP.size()<<std::endl;
	for(int i=0;i< allORP.size();i++){allORP.at(i)->toString();}
	rs->simulate(10);
    return 0;
}