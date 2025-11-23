// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "Seed.h"

#include "Root.h"
#include "rootparameter.h"
#include "Stem.h"

namespace CPlantBox {

/**
 * todo docme
 */
Seed::Seed(int id, std::shared_ptr<const OrganSpecificParameter> param, bool alive, bool active, double age, double length, bool moved, int oldNON)
		:Organ(id, param, alive, active, age, length, Vector3d(0,0,1), 0, moved, oldNON)
{ }


/**
 * todo docme
 */
Seed::Seed(std::shared_ptr<Organism> plant)
		:Organ(plant, nullptr, Organism::ot_seed, 0, 0., 0)
{
	addNode(param()->seedPos, 0.); // realize() is called in Organ constructor
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
	bool plantBox = stemP.size()>1; // prototype + a real parameter definition
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
	const double maxT = 300.; // maximal simulation time
	auto sp = this->param(); // rename (SeedSpecificParameter)

	// Taproot
    try {
        auto rrp = p->getOrganRandomParameter(Organism::ot_root, tapType); // fix tap root angle to 0
        auto rrp1 = std::static_pointer_cast<RootRandomParameter>(rrp);
        double theta_bu = rrp1->theta;
        rrp1->theta = 0;
        std::shared_ptr<Organ> taproot = createRoot(plant.lock(), tapType, 0); // tap root has root type 1
        taproot->addNode(getNode(0), getNodeId(0), 0);
        this->addChild(taproot);
        rrp1->theta = theta_bu; // in case we copy parameter set for basal or shootborne
    } catch (...) {
        std::cout << "Seed::initialize: Parameter set for tap root of subType " << tapType << " not found \n" << std::flush;
        throw;
    }

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
			basalType = tapType;
		}
		int maxB = (sp->maxB);
		if (sp->delayB > 0) { // limit if possible
			maxB = std::min(maxB,int(ceil((maxT-sp->firstB)/sp->delayB))); // maximal for simtime maxT
		}
		double delay = sp->firstB;
		for (int i=0; i<maxB; i++) {
			std::shared_ptr<Organ> basalroot = createRoot(plant.lock(), basalType,  delay);
			basalroot->addNode(getNode(0), getNodeId(0), delay);
			this->addChild(basalroot);
			delay += sp->delayB;
		}
	}

	// Shoot borne roots
    int st = getParamSubType(Organism::ot_root, "shootborne");
    if (st>0) {
        shootborneType = st;
    } // otherwise stick with default
    if ((sp->nC>0) && (sp->firstSB+sp->delaySB<maxT)) { // only if there are any shootborne roots
        std::cout << "Seed::initialize: Shoot borne definition is DEPRICATED, shoot borne roots will be handeled like basal roots \n";
        try {
            p->getOrganRandomParameter(Organism::ot_root, shootborneType); // if the type is not defined an exception is thrown
        } catch (...) {
            if (verbose) {
                std::cout << "Seed::initialize: Shootborne root type #" << shootborneType << " was not defined, using tap root parameters instead\n";
            }
            shootborneType = tapType;
        }
		numberOfRootCrowns = ceil((maxT-sp->firstSB)/sp->delayRC); // maximal number of root crowns
//        Vector3d sbpos = sp->seedPos;
//        sbpos.z=sbpos.z/2.; // half way up the mesocotyl
		double delay = sp->firstSB;
		for (int i=0; i<numberOfRootCrowns; i++) {
			std::shared_ptr<Organ>  shootborne0 = createRoot(plant.lock(), shootborneType, delay);
//			 //TODO fix the initial radial heading
//			 shootborne0->addNode(sbpos,delay);
			shootborne0->addNode(getNode(0), getNodeId(0), delay);
            this->addChild(shootborne0);
            delay += sp->delaySB;
            for (int j=1; j<sp->nC; j++) {
                std::shared_ptr<Organ>  shootborne = createRoot(plant.lock(), shootborneType, delay);
//                // TODO fix the initial radial heading
//                shootborne->addNode(shootborne0->getNode(0), shootborne0->getNodeId(0),delay);
                shootborne->addNode(getNode(0), getNodeId(0), delay);
                this->addChild(shootborne);
                delay += sp->delaySB;
            }
//				sbpos.z+=sp->nz;  // move up, for next root crown
				delay = sp->firstSB + i*sp->delayRC; // reset age
			}
		} else {
			numberOfRootCrowns = 0;
		}


	/*
	 * Create Stem
	 */
	if (plantBox) { // i.e. if a stem is defined
	    std::shared_ptr<Organ> mainstem;
		mainstem = createStem(plant.lock(), mainStemType,0.); // main stem has subtype 1
		mainstem->addNode(Vector3d(0.,0.,0.), getNodeId(0), 0);
		children.push_back(mainstem);
		if (sp->maxTil>0) { // Optional tillers
			try {
				p->getOrganRandomParameter(Organism::ot_stem, tillerType);
			} catch (...) {
				if (verbose) {
					std::cout << "Tiller stem type #" << tillerType << " was not defined, using main stem parameters instead, ";
				}
				auto tillParam = p->getOrganRandomParameter(Organism::ot_stem, 1)->copy(plant.lock());
				tillParam->subType = tillerType;
				p->setOrganRandomParameter(tillParam);
			}
			int maxTi = sp->maxTil;
			double delay = sp->firstTil;
			for (int i=0; i<maxTi; i++) {
				std::shared_ptr<Organ> tiller = createStem(plant.lock(), tillerType,  delay);
				tiller->addNode(Vector3d(0.,0.,0.), getNodeId(0), 0);
				children.push_back(tiller);
				delay += sp->delayTil;
			}
		}
	}
}

/**
 * Creates a deep copy of the seeds child organs.
 * Use: the rootsystem manages base root itself, and just uses Seed class for initialization
 */
std::vector<std::shared_ptr<Organ>> Seed::copyBaseOrgans(std::shared_ptr<Organism> plant)
{
	std::vector<std::shared_ptr<Organ>> organs;
	for (auto& o : children) {
		organs.push_back(o->copy(plant));
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
 * additionally, use param()->toString() and getOrganRandomParameter()->toString() to obtain all information.
 */
std::string Seed::toString() const
{
	std::stringstream newstring;
	newstring << "; maximal number of basals: " << param()->maxB << ", of shootborne " << param()->nC  << ", of tillers " << param()->maxTil << ".";
	return Organ::toString() + newstring.str();
}

/**
 * todo doc
 */
std::shared_ptr<Organ> Seed::createRoot(std::shared_ptr<Organism> plant, int type, double delay)
{
	return std::make_shared<Root>(plant, type,delay, shared_from_this(), 0);
}

/**
 * todo doc// overwrite if you want to change the types
 */
std::shared_ptr<Organ> Seed::createStem(std::shared_ptr<Organism> plant, int type,  double delay)
{
	return std::make_shared<Stem>(plant, type,delay, shared_from_this(), 0);
}


} // namespace CPlantBox
