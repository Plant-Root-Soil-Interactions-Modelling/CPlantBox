#include "../../../src/structural/MappedOrganism.h"
#include "../../../src/functional/Photosynthesis.h"
#include <iostream>
#include <string> 

int main(int argc, char **argv) 
{
	std::string path = "../../../modelparameter/structural/plant/";
    std::string name = "Triticum_aestivum_adapted_2023";
    auto rs = std::make_shared<CPlantBox::MappedPlant>(2);//pb.RootSystem()  # the original
	rs->openXML(path + name + ".xml");
	
	
	auto sdf = std::make_shared<CPlantBox::SDF_PlantBox>(100., 100., 60. );

	rs->setGeometry(sdf) ;// creates soil space to stop roots from growing out of the soil

	
	rs->initialize(false);
	rs->simulate(10,false);
	

		
	//	constant total potential (hydraulic equilibrium)
	
	auto photo = std::make_shared<CPlantBox::Photosynthesis>(rs);
	std::vector<std::vector<double>> ages(0, std::vector<double> (0, 0));
	std::vector<double> krRoot{ 1.73e-5 };
	std::vector<double> krStem{0.};
	std::vector<double> krLeaf{ 3.83e-4};
	std::vector<std::vector<double>> krs{ krRoot,krStem,krLeaf};
	photo->setKr(krs,ages);
	std::vector<double> kxRoot{4.32e-2  };
	std::vector<double> kxStem{5.13};
	std::vector<double> kxLeaf{ 1.27};
	std::vector<std::vector<double>> kxs{ kxRoot,kxStem,kxLeaf};
	photo->setKx(kxs,ages);
	double RH_ = 0.5;
	double TairC_ = 25;
	double es = 6.112 * std::exp((17.67 * TairC_) / (TairC_ + 243.5));
    double ea = es * RH_ ;
	
	// prepare soil matric potentials per segment
	double p_top = -500;  // top soil pressure [cm]
	auto segs = rs->segments; // MappedRootSystem has access to segments and nodes 
	auto nodes = rs->nodes;
	auto organTypes = rs->organTypes;
	auto p_s = std::vector<double>(segs.size(),0.);// soil total potentials around each root segment
	for(int i=0;i<segs.size();i++)
	{
		int ot = organTypes.at(i);
		if(ot ==2)
		{
			auto s = segs.at(i);//, s in enumerate(segs):
			p_s.at(i) = p_top - 0.5 * (nodes.at(s.x).z + nodes.at(s.y).z) ;
		}else{
			p_s.at(i) = photo->psi_air;
		}
	}	
	
	
	photo->solve_photosynthesis(ea, es,
				1., p_s , false);
	std::cout <<path + name + ".xml"<<std::endl;
    return 0;
}