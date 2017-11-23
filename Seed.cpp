#include "Seed.h"

Seed::Seed(Plant* plant) :Organ(plant, nullptr, 0, 0)
{
	root_param = (SeedParameter*)plant->getOrganTypeParameter(Plant::ot_seed,0)->realize();
}

/**
 * returns the type of the organ
 */
int Seed::organType() const
{
	return Plant::ot_seed;
}

/**
 *
 */
void Seed::initialize()
{
	// Create root system
	const double maxT = 365.; // maximal simulation time
	SeedParameter* sparam = (SeedParameter*)root_param;  // rename

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
		if (plant->getOrganTypeParameter(Plant::ot_root,basalType)->subType<1) { // if the type is not defined, copy tap root
			std::cout << "Basal root type #" << basalType << " was not defined, using tap root parameters instead\n";
			RootTypeParameter* tapParam = (RootTypeParameter*)plant->getOrganTypeParameter(Plant::ot_root, 0);
			RootTypeParameter* brtp = new RootTypeParameter(*tapParam);
			brtp->subType = basalType;
			plant->setOrganTypeParameter(brtp);
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
 *
 */
void Seed::simulate(double dt, bool silence)
{
	for (auto& c : children)  {
		c->simulate(dt);
	}
}


/**
 * Quick info about the object for debugging
 */
std::string Seed::toString() const
{
	std::stringstream str;
	str << "Seed #"<< id <<": type "<< root_param->subType << ", length: "<< length << ", age: " << age;
	return str.str();
}
