// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "Seed.h"

#include "Root.h"
#include "Stem.h"

namespace CPlantBox {

/**
 * todo docme
 */
Seed::Seed(std::shared_ptr<Organism> plant)
:Organ(plant, nullptr, Organism::ot_seed, 0, 0., Vector3d(), 0., 0)
{
	// realize() is called in Organ constructor
	addNode(param()->seedPos, 0.);
}

/**
 * Deep copies the organ into the new plant @param rs.
 * All laterals are deep copied, plant and parent pointers are updated.
 *
 * @param plant     the plant the copied organ will be part of
 */
std::shared_ptr<Organ> Seed::copy(std::shared_ptr<Organism> rs)
{
	auto s = std::make_shared<Seed>(*this); // shallow copy
	s->parent = std::weak_ptr<Organ>();
	s->plant = rs;
	s->param_ = std::make_shared<SeedSpecificParameter>(*param()); // copy parameters
	for (size_t i=0; i< children.size(); i++) {
		s->children[i] = children[i]->copy(rs); // copy laterals
		s->children[i]->setParent(s);
	}
	return s;
}

/**
 * Creates the initial organs,
 * i.e. taproot, basal root, (if needed) shoot borne root (if plant) stem
 */
void Seed::initialize(bool verbose)
{
	auto p = plant.lock();
	auto stemP = p->getOrganRandomParameter(Organism::ot_stem);
	bool plantBox = stemP.size()>0;
	if (verbose) {
		if (plantBox) {
			std::cout << "Seed::initialize: Plant \n";
		} else {
			std::cout << "Seed::initialize: RootSystem \n";
		}
	}

	/*
	 * Create roots
	 */
	const double maxT = 365.; // maximal simulation time
	auto sp = this->param(); // rename
	Vector3d iheading(0,0,-1);

	// Taproot
	std::shared_ptr<Organ> taproot = createRoot(plant.lock(), tapType, iheading ,0); // tap root has root type 1
	taproot->addNode(getNode(0), getNodeId(0), 0);
	this->addChild(taproot);

	// Basal roots
	int bt = getParamSubType(Organism::ot_root, "basal");
	if (bt>0) {
		basalType = bt;
	} // otherwise stick with default
	if (sp->maxB>0) {
		try {
			p->getOrganRandomParameter(Organism::ot_root, basalType); // if the type is not defined an exception is thrown
		} catch (...) {
			if (verbose) {
				std::cout << "Seed::initialize: Basal root type #" << basalType << " was not defined, using tap root parameters instead\n" << std::flush;
			}
			auto brtp = p->getOrganRandomParameter(Organism::ot_root, 1)->copy(plant.lock());
			brtp->subType = basalType;
			p->setOrganRandomParameter(brtp);
		}
		int maxB = (sp->maxB);
		if (sp->delayB > 0) { // limit if possible
			maxB = std::min(maxB,int(ceil((maxT-sp->firstB)/sp->delayB))); // maximal for simtime maxT
		}
		double delay = sp->firstB;
		for (int i=0; i<maxB; i++) {
			std::shared_ptr<Organ> basalroot = createRoot(plant.lock(), basalType, iheading, delay);
			basalroot->addNode(getNode(0), getNodeId(0), delay);
			this->addChild(basalroot);
			delay += sp->delayB;
		}
	}

	// Shoot borne roots
	if (!plantBox) { // use CRootBox initialization
		int st = getParamSubType(Organism::ot_root, "shootborne");
		if (st>0) {
			shootborneType = st;
		} // otherwise stick with default
		if ((sp->nC>0) && (sp->firstSB+sp->delaySB<maxT)) { // only if there are any shootborne roots
			try {
				p->getOrganRandomParameter(Organism::ot_root, shootborneType); // if the type is not defined an exception is thrown
			} catch (...) {
				if (verbose) {
					std::cout << "Seed::initialize:Shootborne root type #" << shootborneType << " was not defined, using tap root parameters instead\n";
				}
				auto srtp =  p->getOrganRandomParameter(Organism::ot_root, 1)->copy(plant.lock());
				srtp->subType = shootborneType;
				p->setOrganRandomParameter(srtp);
			}
			Vector3d sbpos = sp->seedPos;
			sbpos.z=sbpos.z/2.; // half way up the mesocotyl
			numberOfRootCrowns = ceil((maxT-sp->firstSB)/sp->delayRC); // maximal number of root crowns
			double delay = sp->firstSB;
			for (int i=0; i<numberOfRootCrowns; i++) {
				std::shared_ptr<Organ>  shootborne0 = createRoot(plant.lock(), shootborneType, iheading ,delay);
				// TODO fix the initial radial heading
				shootborne0->addNode(sbpos,delay);
				this->addChild(shootborne0);
				delay += sp->delaySB;
				for (int j=1; j<sp->nC; j++) {
					std::shared_ptr<Organ>  shootborne = createRoot(plant.lock(), shootborneType, iheading ,delay);
					// TODO fix the initial radial heading
					shootborne->addNode(shootborne0->getNode(0), shootborne0->getNodeId(0),delay);
					this->addChild(shootborne);
					delay += sp->delaySB;
				}
				sbpos.z+=sp->nz;  // move up, for next root crown
				delay = sp->firstSB + i*sp->delayRC; // reset age
			}
		} else {
			numberOfRootCrowns = 0;
		}
	}

	/*
	 * Create Stem
	 */
	if (plantBox) { // i.e. if a stem is defined
		// Stem
		Vector3d isheading(0, 0, 1);//Initial Stem heading
		std::shared_ptr<Organ> mainstem = createStem(plant.lock(), mainStemType, isheading, 0.); // main stem has subtype 1
		mainstem->addNode(getNode(0), getNodeId(0), 0);
		children.push_back(mainstem);
		// Optional tillers
		if (sp->maxTil>0) {
			try {
				p->getOrganRandomParameter(Organism::ot_stem, tillerType);
			} catch (...) {
				if (verbose) {
					std::cout << "Tiller stem type #" << tillerType << " was not defined, using main stem parameters instead, ";
				}
				auto tillParam = p->getOrganRandomParameter(Organism::ot_stem, 1)->copy(plant.lock());
				tillParam->subType = basalType;
				p->setOrganRandomParameter(tillParam);
			}
			int maxTi = sp->maxTil;
			if (sp->delayB>0) {
				maxTi = std::min(maxTi,int(std::ceil((maxT-sp->firstB)/sp->delayB))); // maximal for simtime maxT
			}
			double delay = sp->firstB;
			for (int i=0; i<maxTi; i++) {
				std::shared_ptr<Organ> tiller = createStem(plant.lock(), tillerType, isheading, delay);
				tiller->addNode(getNode(0), getNodeId(0), 0);
				children.push_back(tiller);
				delay += sp->delayB;
			}
		}
	}
}

/**
 * Creates a shallow copy of the seeds child organs.
 * Use: the rootsystem manages base root itself, and just uses Seed class for initialization
 */
std::vector<std::shared_ptr<Organ>> Seed::copyBaseOrgans()
{
	std::vector<std::shared_ptr<Organ>> organs;
	for (auto& o : children) {
		organs.push_back(o->copy(plant.lock()));
	}
	return organs;
}

/*
 * Searches for a parameter with name @param str, and returns its subtype
 */
int Seed::getParamSubType(int organtype, std::string str)
{
	auto orp = plant.lock()->getOrganRandomParameter(organtype);
	for (auto& o :orp) {
		if (o->name == str) {
			return o->subType;
		}
	}
	return -1;
}

/**
 * Quick info about the object for debugging
 * additionally, use getParam()->toString() and getOrganRandomParameter()->toString() to obtain all information.
 */
std::string Seed::toString() const
{
	std::string str = Organ::toString();
	str.replace(0, 5, "Seed");
	std::stringstream newstring;
	newstring << "; maximal number of basals: " << param()->maxB << ", of shootborne " << param()->nC  << ", of tillers " << param()->maxTil << ".";
	str.replace(str.size()-1, 1, newstring.str());
	return str;
}

/**
 * todo doc
 */
std::shared_ptr<Organ> Seed::createRoot(std::shared_ptr<Organism> plant, int type, Vector3d heading, double delay)
{
	return std::make_shared<Root>(plant, type, heading, delay, shared_from_this(), 0, 0);
}

/**
 * todo doc// overwrite if you want to change the types
 */
std::shared_ptr<Organ> Seed::createStem(std::shared_ptr<Organism> plant, int type, Vector3d heading, double delay)
{
	return std::make_shared<Stem>(plant, type, heading, delay, shared_from_this(), 0, 0);
}

} // namespace CPlantBox
