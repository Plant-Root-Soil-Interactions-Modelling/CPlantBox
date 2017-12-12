#include "Seed.h"

Seed::Seed(Plant* plant) :Organ(plant, nullptr, 0, 0), seed_pos(Vector3d(0,0,-3))
{
	param = (SeedParameter*)plant->getParameter(Organ::ot_seed,0)->realize();
}

/**
 *
 */
void Seed::initialize()
{
	// Create root system
	const double maxT = 365.; // maximal simulation time
	SeedParameter* sparam = (SeedParameter*)param;  // rename

	// Taproot
	Vector3d iheading(0,0,-1);
	Root* taproot = new Root(plant, this, 1, 0., iheading ,0., 0.); // tap root has subtype 1
	taproot->addNode(sparam->seedPos,0);
	children.push_back(taproot);

	Vector3d isheading(0,0,1);//Initial Stem heading

	Stem* mainstem = new Stem(plant, this, 1, 0., isheading ,0., 0.); // tap root has subtype 1
	mainstem->addNode(sparam->seedPos,0);
	children.push_back(mainstem);

	//	 Basal roots
	if (sparam->maxB>0) {
		if (plant->getParameter(Organ::ot_root,basalType)->subType<1) { // if the type is not defined, copy tap root
			std::cout << "Basal root type #" << basalType << " was not defined, using tap root parameters instead\n";
			RootTypeParameter* tapParam = (RootTypeParameter*)plant->getParameter(Organ::ot_root, 0);
			RootTypeParameter* brtp = new RootTypeParameter(*tapParam);
			brtp->subType = basalType;
			plant->setParameter(brtp);
		}
		int maxB = sparam->maxB;
		if (sparam->delayB>0) {
			maxB = std::min(maxB,int(std::ceil((maxT-sparam->firstB)/sparam->delayB))); // maximal for simtime maxT
		}
		double delay = sparam->firstB;
		for (int i=0; i<maxB; i++) {
			Root* basalroot = new Root(plant, this, basalType, delay, iheading ,0., 0.);
			basalroot->r_nodes.push_back(sparam->seedPos); // node
			basalroot->nodeIDs.push_back(taproot->nodeIDs.at(0)); // tap root ID
			basalroot->nctimes.push_back(delay); // exact creation time
			children.push_back(basalroot);
			delay += sparam->delayB;
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
