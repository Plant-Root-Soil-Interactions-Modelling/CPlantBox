#include "Seed.h"

namespace CPlantBox {


Seed::Seed(Plant* plant) :Organ(plant, nullptr, 0, 0), seed_pos(Vector3d(0,0,0))
{
	param = (SeedParameter*)plant->getParameter(Organ::ot_seed,0)->realize();
}

SeedParameter* Seed::initializeparam()
{
	// Create root system

	SeedTypeParameter* stp = (SeedTypeParameter*)plant->getParameter(Organ::ot_seed, 0);
	param = stp->realize(); // throw the dice
	SeedParameter* sparam = (SeedParameter*)param;
	std::cout << "maxb = " << sparam->maxB << std::endl;
	std::cout << "maxTi in type parameter = " << sparam->maxTi << std::endl;
	return sparam;
}



/**
 *
 */
void Seed::initialize(SeedParameter* sparam)
{
	const double maxT = 365.; // maximal simulation time
	Vector3d iheading(0,0,-1);
	if (Plant::noParamFile[1] == 1) {
		std::cout<<"no root param"<<std::endl;
	} else {
		Root* taproot = new Root(plant, this, 1, 0, iheading ,0., 0.); // tap root has subtype 1
		taproot->addNode(sparam->seedPos,0);
		std::cout << "sparam->seedPos is" << sparam->seedPos.toString() << "\n" ;
		children.push_back(taproot);
    		if (sparam->maxB>0) {
			if (plant->getParameter(Organ::ot_root, basalType)->subType=4) { // if the type is not defined, copy tap root
				std::cout << "Basal root type #" << basalType << " was not defined, using tap root parameters instead\n";
//				RootTypeParameter* tapParam = (RootTypeParameter*)plant->getParameter(Organ::ot_root, 1);
//				RootTypeParameter* brtp = new RootTypeParameter(*tapParam);
//				brtp->subType = basalType;
//				plant->setParameter(brtp);
                for (int i=0; i<sparam->maxB; i++) {
				Root* basalroot = new Root(plant, this, 1, sparam->firstB, iheading ,0., 0.);
//				basalroot->r_nodes.push_back(mainstem->r_nodes.at(0)); // node
//				basalroot->nodeIDs.push_back(mainstem->nodeIDs.at(0)); // tap root ID
//				basalroot->nctimes.push_back(delay); // exact creation time
//				children.push_back(mainstem);
//				delay += sparam->delayB;
                basalroot->addNode(sparam->seedPos,0);
                children.push_back(basalroot);
								std::cout << "no basal root type 4= " << sparam->maxB << "\n";
                }
			//	std::cout << "maxb = " << sparam->maxB << "\n";

			}else{
			int maxB = sparam->maxB;
			if (sparam->delayB>0) {
				maxB = std::min(maxB,int(std::ceil((maxT-sparam->firstB)/sparam->delayB))); // maximal for simtime maxT
			}
			double delay = sparam->firstB;
			for (int i=0; i<maxB; i++) {
				Root* basalroot = new Root(plant, this, 1, delay, iheading ,0., 0.);
//				basalroot->r_nodes.push_back(mainstem->r_nodes.at(0)); // node
//				basalroot->nodeIDs.push_back(mainstem->nodeIDs.at(0)); // tap root ID
//				basalroot->nctimes.push_back(delay); // exact creation time
//				children.push_back(mainstem);
//				delay += sparam->delayB;
                basalroot->addNode(sparam->seedPos,0);
                children.push_back(basalroot);
			}
			}
		}
	}


	//main stem initialation

	if (Plant::noParamFile[2] == 1) {
		std::cout<<"no stem parameter, maybe the XML based parameter file is directly converted from the old rparam file"<<std::endl;
	} else {
		//	Vector3d isheading(0,0,1);//Initial Stem heading
		Vector3d isheading(0, 0, 1);//Initial Stem heading
		Stem* mainstem = new Stem(plant, this, 1, 0., isheading, 0., 0.); // tap root has subtype 1
		mainstem->addNode(sparam->seedPos, 0);
		children.push_back(mainstem);


		if (sparam->maxTi>0) {
			if (plant->getParameter(Organ::ot_stem, tillerType)->subType<1) { // if the type is not defined, copy tap root
				std::cout << "Basal root type #" << basalType << " was not defined, using tap root parameters instead\n";
				StemTypeParameter* tillParam = (StemTypeParameter*)plant->getParameter(Organ::ot_stem, 1);
//				StemTypeParameter* titp = new StemTypeParameter(*tillParam);
//				titp->subType = tillerType;
//				plant->setParameter(titp);
				std::cout << "default maxT type is main stem = " << sparam->maxTi << "\n";

			} else{
			int maxTi = sparam->maxTi;
			if (sparam->delayB>0) {
				maxTi = std::min(maxTi,int(std::ceil((maxT-sparam->firstB)/sparam->delayB))); // maximal for simtime maxT
			}
			std::cout << "maxT = " << sparam->maxTi << "\n";
			double delay = sparam->firstB;
			StemTypeParameter* tillParam = (StemTypeParameter*)plant->getParameter(Organ::ot_stem, 4);
			for (int i=0; i<maxTi; i++) {
				Stem* tiller = new Stem(plant, this, 4, delay, isheading ,0., 0.);
                tiller->addNode(sparam->seedPos,0);
                children.push_back(tiller);

				std::cout << "new maxT type is main stem = " << sparam->maxTi << "\n";
			}

		}

		//	Stem* tiller1 = new Stem(plant, this, 1, 2, isheading ,0., 0.); // tap root has subtype 1
		//	tiller1->addNode(sparam->seedPos,0);
		//	children.push_back(tiller1);
		//
		//	Stem* tiller2 = new Stem(plant, this, 1, 4 , isheading ,0., 0.); // tap root has subtype 1
		//	tiller2->addNode(sparam->seedPos,0);
		//	children.push_back(tiller2);
		//	Stem* tiller3 = new Stem(plant, this, 1, 5, isheading ,0., 0.); // tap root has subtype 1
		//	tiller3->addNode(sparam->seedPos,0);
		//	children.push_back(tiller3);
	}

}
}

/**
 * Quick info about the object for debugging
 */
std::string Seed::toString() const
{
	std::stringstream str;
	str << "Seed #"<< id <<": type "<< param->subType << ", length: "<< length << ", age: " << age;
	return str.str();
}

} // namespace CPlantBox
