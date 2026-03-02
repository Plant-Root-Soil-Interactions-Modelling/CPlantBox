// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "Seed.h"

#include "Root.h"
#include "Stem.h"
#include "rootparameter.h"

namespace CPlantBox {

/**
 * Creates a seed from scratch by given parameters
 */
Seed::Seed(int id, std::shared_ptr<const OrganSpecificParameter> param, bool alive, bool active, double age, double length, bool moved, int oldNON)
    : Organ(id, param, alive, active, age, length, Vector3d(0, 0, 1), 0, moved, oldNON) {}

/**
 * Creates a seed from a Organism using its parameter set
 */
Seed::Seed(std::shared_ptr<Organism> plant) : Organ(plant, nullptr, Organism::ot_seed, 0, 0., 0) {
    addNode(param()->seedPos, 0.); // realize() is called within the Organ constructor
}

/**
 * Deep copies the organ into the new plant @param rs.
 * All laterals are deep copied, plant and parent pointers are updated.
 *
 * @param plant     the plant the copied organ will be part of
 */
std::shared_ptr<Organ> Seed::copy(std::shared_ptr<Organism> rs) {
    auto s = std::make_shared<Seed>(*this); // shallow copy
    s->parent = std::weak_ptr<Organ>();
    s->plant = rs;
    s->param_ = std::make_shared<SeedSpecificParameter>(*param()); // copy parameters
    for (size_t i = 0; i < children.size(); i++) {
        s->children[i] = children[i]->copy(rs); // copy laterals
        s->children[i]->setParent(s);
    }
    return s;
}

/**
 * Initialization creates the initial organs,
 * i.e. taproot, basal root, (if needed) shoot borne root (if plant) stem
 */
void Seed::initialize(bool verbose) {

    auto p = plant.lock();
    bool plantBox;
    try {
        auto stemP = p->getOrganRandomParameter(Organism::ot_stem, mainStemType);
        plantBox = true;
    } catch (...) {
        plantBox = false;
    }
    if (verbose) {
        if (plantBox) {
            std::cout << "Seed::initialize: Plant \n";
        } else {
            std::cout << "Seed::initialize: RootSystem \n";
        }
    }

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
    if (bt > 0) {
        basalType = bt;
    } // otherwise stick with default
    if (basalRoots()) {
        try {
            p->getOrganRandomParameter(Organism::ot_root, basalType); // if the type is not defined an exception is thrown
        } catch (...) {
            if (verbose) {
                std::cout << "Seed::initialize: Basal root type #" << basalType << " was not defined, using tap root parameters instead\n" << std::flush;
            }
            auto brtp = p->getOrganRandomParameter(Organism::ot_root, tapType)->copy(plant.lock());
            brtp->subType = basalType;
            p->setOrganRandomParameter(brtp);
        }
        int maxB = (sp->maxB);
        if (sp->delayB > 0) {                                                              // limit if possible
            maxB = std::min(maxB, int(ceil((getMaxSimTime() - sp->firstB) / sp->delayB))); // maximal for simtime getMaxT()
        }
        double delay = sp->firstB;
        for (int i = 0; i < maxB; i++) {
            std::shared_ptr<Organ> basalroot = createRoot(plant.lock(), basalType, delay);
            basalroot->addNode(getNode(0), getNodeId(0), delay);
            this->addChild(basalroot);
            delay += sp->delayB;
        }
    }

    // Shoot borne roots
    int st = getParamSubType(Organism::ot_root, "shootborne");
    if (st > 0) {
        shootborneType = st;
    } // otherwise stick with default
    if (shootBorneRoots()) { // only if there are any shootborne roots
        // std::cout << "Seed::initialize: Shoot borne definition is DEPRICATED\n"; // but we are stuck with it for now
        try {
            p->getOrganRandomParameter(Organism::ot_root, shootborneType); // if the type is not defined an exception is thrown
        } catch (...) {
            if (verbose) {
                std::cout << "Seed::initialize: Shootborne root type #" << shootborneType << " was not defined, using tap root parameters instead\n";
            }
            auto srtp = p->getOrganRandomParameter(Organism::ot_root, 1)->copy(plant.lock());
            srtp->subType = shootborneType;
            p->setOrganRandomParameter(srtp);
        }
        rootCrowns_ = 1; // ceil((getMaxSimTime() - sp->firstSB) / sp->delayRC); (old model) // maximal number of root crowns
        double delay = sp->firstSB;
        // Vector3d sbpos = sp->seedPos; // (old model)
        // sbpos.z=sbpos.z/2.; // (old model) half way up the mesocotyl
        for (int i = 0; i < rootCrowns_; i++) {
            std::shared_ptr<Organ> shootborne0 = createRoot(plant.lock(), shootborneType, delay);
            // TODO fix the initial radial heading
            // shootborne0->addNode(sbpos,delay); // (old model)
            shootborne0->addNode(getNode(0), getNodeId(0), delay);
            this->addChild(shootborne0);
            delay += sp->delaySB;
            for (int j = 1; j < sp->nC; j++) {
                std::shared_ptr<Organ> shootborne = createRoot(plant.lock(), shootborneType, delay);
                // TODO fix the initial radial heading
                shootborne->addNode(getNode(0), getNodeId(0), delay);
                this->addChild(shootborne);
                delay += sp->delaySB;
            }
            //	sbpos.z+=sp->nz;  // (old model) move up, for next root crown
            delay = sp->firstSB + i * sp->delayRC; // reset age
        }
    } else {
        rootCrowns_ = 0;
    }

    // Stem
    if (plantBox) { // i.e. if a stem is defined
        std::shared_ptr<Organ> mainstem;
        mainstem = createStem(plant.lock(), mainStemType, 0.); // main stem has subtype 1
        mainstem->addNode(Vector3d(0., 0., 0.), getNodeId(0), 0);
        children.push_back(mainstem);
        if (tillers()) { // Optional tillers
            try {
                p->getOrganRandomParameter(Organism::ot_stem, tillerType);
            } catch (...) {
                if (verbose) {
                    std::cout << "Seed::initialize: Tiller stem type #" << tillerType << " was not defined, using main stem parameters instead";
                }
                auto tillParam = p->getOrganRandomParameter(Organism::ot_stem, 1)->copy(plant.lock());
                tillParam->subType = tillerType;
                p->setOrganRandomParameter(tillParam);
            }
            int maxTi = sp->maxTil;
            if (sp->delayTil > 0) {                                                                  // limit if possible
                maxTi = std::min(maxTi, int(ceil((getMaxSimTime() - sp->firstTil) / sp->delayTil))); // maximal for simtime getMaxT()
            }
            double delay = sp->firstTil;
            for (int i = 0; i < maxTi; i++) {
                std::shared_ptr<Organ> tiller = createStem(plant.lock(), tillerType, delay);
                tiller->addNode(Vector3d(0., 0., 0.), getNodeId(0), 0);
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
std::vector<std::shared_ptr<Organ>> Seed::copyBaseOrgans(std::shared_ptr<Organism> plant) {
    std::vector<std::shared_ptr<Organ>> organs;
    for (auto &o : children) {
        organs.push_back(o->copy(plant));
    }
    return organs;
}

/**
 * Checks if there is a basal root definition in the parameter set. 
 * This works only before initialization, because on initialization the basal root type is copied from to the tap root type
 */
bool Seed::basalRootDefined() {
    int bt = getParamSubType(Organism::ot_root, "basal");
    if (bt <= 0) {
        bt = basalType; // stick with default if not defined
    }
    try {
        plant.lock()->getOrganRandomParameter(Organism::ot_root, basalType); // if the type is not defined an exception is thrown
        return true;
    } catch (...) {
        return false;
    }
}

/**
 * Checks if there is a shoot borne root definition in the parameter set. 
 * This works only before initialization, because on initialization the shoot borne root type is copied from to the basal root type
 */
bool Seed::shootBorneRootDefined() { 
	int st = getParamSubType(Organism::ot_root, "shootborne");
	if (st <= 0) {
		st = shootborneType; // stick with default if not defined
	}
	try {
		plant.lock()->getOrganRandomParameter(Organism::ot_root, shootborneType); // if the type is not defined an exception is thrown
		return true;
	} catch (...) {
		return false;
	}
}

/**
 * Checks if there is a tiller definition in the parameter set. 
 * This works only before initialization, because on initialization the tiller type is copied from to the main stem type
 */
bool Seed::tillersDefined() { 
	int tt = getParamSubType(Organism::ot_stem, "tiller");
	if (tt <= 0) {
		tt = tillerType; // stick with default if not defined
	}
	try {
		plant.lock()->getOrganRandomParameter(Organism::ot_stem, tt); // if the type is not defined an exception is thrown
		return true;
	} catch (...) {
		return false;
	}
}

/*
 * Searches for a parameter with name @param str, and returns its subtype, or -1 if not found.
 */
int Seed::getParamSubType(int organtype, std::string str) {
    auto orp = plant.lock()->getOrganRandomParameter(organtype);
    for (auto &o : orp) {
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
std::string Seed::toString() const {
    std::stringstream newstring;
    newstring << "; maximal number of basals: " << param()->maxB << ", of shootborne " << param()->nC << ", of tillers " << param()->maxTil << ".";
    return Organ::toString() + newstring.str();
}

/**
 * todo doc
 */
std::shared_ptr<Organ> Seed::createRoot(std::shared_ptr<Organism> plant, int type, double delay) {
    return std::make_shared<Root>(plant, type, delay, shared_from_this(), 0);
}

/**
 * todo doc// overwrite if you want to change the types
 */
std::shared_ptr<Organ> Seed::createStem(std::shared_ptr<Organism> plant, int type, double delay) {
    return std::make_shared<Stem>(plant, type, delay, shared_from_this(), 0);
}

} // namespace CPlantBox
