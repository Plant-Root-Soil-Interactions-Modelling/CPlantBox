#include "../../../src/structural/MappedOrganism.h"
#include <iostream>
#include <string> 

int main(int argc, char **argv) 
{
    auto rs = std::make_shared<CPlantBox::MappedPlant>(2);

    std::string path = "../../../modelparameter/structural/plant/";
    std::string name ="testLatDef";
	//"Triticum_aestivum_adapted_2023";
	//"test_relcoord";
	//"manyleaves";//
			
	std::cout<<(path+name+".xml")<<std::endl;
	rs->readParameters(path+name+".xml");
    rs->initialize(false);//false, false);
	int steps = 50.;
	double dt =1.;
	for(int i=0;i<steps;i++)
	{	

		rs->simulate(dt, true);
		auto oneRoot = rs->getOrgans(2).at(0);
		auto orrp = oneRoot->getOrganRandomParameter();
		std::cout<< orrp->toString();
		auto oneStem = rs->getOrgans(3).at(0);
		auto osrp = oneStem->getOrganRandomParameter();
		std::cout<< osrp->toString();
		auto oneLeaf = rs->getOrgans(4).at(0);
		auto olrp = oneLeaf->getOrganRandomParameter();
		std::cout<< olrp->toString();
		rs->write("test_"+std::to_string(i)+".vtp");
		auto allLeaves = rs->getOrgans(4);
		std::cout<<"get leaves pni "<<std::endl;
		for(size_t lnum = 0; lnum < allLeaves.size(); lnum ++)
		{
			std::cout<<allLeaves.at(lnum)->parentNI<<" ";
		}std::cout<<std::endl;
		
	}
	
    return 0;
}