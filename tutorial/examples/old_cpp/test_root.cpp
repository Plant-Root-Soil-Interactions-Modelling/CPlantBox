#include "../../../src/structural/MappedOrganism.h"
#include <iostream>
#include <string> 

int main(int argc, char **argv) 
{
	auto plant = std::make_shared<CPlantBox::Organism>(2);
	auto p0 = std::make_shared<CPlantBox::RootRandomParameter>(plant);
	p0->name = "taproot";
	p0->subType=1;
	p0->la=10.; p0->lb=1.; p0->lmax=100.; p0->ln=1.; p0->r=1.5; p0->dx =0.5;
	std::vector<std::vector<int>> vec(1, std::vector<int> (1, 2));
	std::vector<std::vector<double>> vecd(1, std::vector<double> (1, 1.));
	p0->successorST = vec;
	p0->successorP = vecd; 
	auto p1 = std::make_shared<CPlantBox::RootRandomParameter>(plant);
	p1->name ="lateral";
	p1->subType=2; p1->lmax=25; p1->r=2.; p1->dx = 0.1;
	//self.p0, self.p1 = p0, p1  # needed at later point
	plant->setOrganRandomParameter(p0);//  # the organism manages the type parameters and takes ownership
	plant->setOrganRandomParameter(p1);
	auto srp = std::make_shared<CPlantBox::SeedRandomParameter>(plant);
	plant->setOrganRandomParameter(srp);

	auto param0 = std::static_pointer_cast<CPlantBox::RootSpecificParameter>(p0->realize()) ;// # set up root by hand (without a root system)
	param0->la=0; param0->lb =  0 ;// # its important parent has zero length, otherwise creation times are messed up
	auto parentroot = std::make_shared<CPlantBox::Root>(1, param0, true, true, 0., 0., CPlantBox::Vector3d(0, 0, -1), 0, false, 0) ;// # takes ownership of param0
	parentroot->setOrganism(plant);
	parentroot->addNode(CPlantBox::Vector3d(0, 0, -3), 0) ;// # there is no nullptr in Python

	//self.parentroot = parentroot ;// # store parent (not owned by child Organ)
	auto root =std::make_shared<CPlantBox::Root>(plant, p0->subType,  0, parentroot , 0);
	root->setOrganism(plant);
	
	bool verbose = true;
	root->simulate(.5, verbose);
	std::cout<<"r.hasMoved() "<<root->hasMoved()<<" is False?"<<std::endl;
	root->simulate(1e-1, verbose);
	std::cout<<"r.hasMoved() "<<root->hasMoved()<<" is True?"<<std::endl;
	//non = r.getNumberOfNodes();
	root->simulate(2.4, verbose);
	
    return 0;
}