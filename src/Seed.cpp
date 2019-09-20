// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "Seed.h"
#include "Root.h"
#include "Stem.h"

namespace CPlantBox {

/**
 *
 */
void Seed::initialize()
{
    auto stemP = plant->getOrganRandomParameter(Organism::ot_stem);
    bool plantBox = stemP.size()>1;
    if (plantBox) {
        std::cout << "Seed::initialize: initializing plant \n";
    } else {
        std::cout << "Seed::initialize: initializing root system \n";
    }

    /*
     * Create roots
     */
    const double maxT = 365.; // maximal simulation time
    const SeedSpecificParameter* sp = this->param(); // rename
    Vector3d iheading(0,0,-1);

    // Taproot
    Organ* taproot = createRoot(plant, tapRootType, iheading ,0, nullptr, 0, 0); // tap root has root type 1
    taproot->addNode(sp->seedPos,0);
    this->addChild(taproot);

    // Basal roots
    int bt = getParamSubType(Organism::ot_root, "basal");
    if (bt>0) {
        basalType = bt;
    } // otherwise stick with default
    if (sp->maxB>0) {
        try {
            plant->getOrganRandomParameter(Organism::ot_root, basalType); // if the type is not defined an exception is thrown
        } catch (...) {
            std::cout << "Seed::initialize: Basal root type #" << basalType << " was not defined, using tap root parameters instead\n" << std::flush;
            RootRandomParameter* brtp = (RootRandomParameter*)plant->getOrganRandomParameter(Organism::ot_root, 1)->copy(plant);
            brtp->subType = basalType;
            plant->setOrganRandomParameter(brtp);
        }
        int maxB = (sp->maxB);
        if (sp->delayB > 0) { // limit if possible
            maxB = std::min(maxB,int(ceil((maxT-sp->firstB)/sp->delayB))); // maximal for simtime maxT
        }
        double delay = sp->firstB;
        for (int i=0; i<maxB; i++) {
            Organ* basalroot = createRoot(plant, basalType, iheading, delay, nullptr, 0, 0);
            basalroot->addNode(taproot->getNode(0), taproot->getNodeId(0), delay);
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
                plant->getOrganRandomParameter(Organism::ot_root, shootborneType); // if the type is not defined an exception is thrown
            } catch (...) {
                std::cout << "Seed::initialize:Shootborne root type #" << shootborneType << " was not defined, using tap root parameters instead\n";
                RootRandomParameter* srtp =  (RootRandomParameter*)plant->getOrganRandomParameter(Organism::ot_root, 1)->copy(plant);
                srtp->subType = shootborneType;
                plant->setOrganRandomParameter(srtp);
            }
            Vector3d sbpos = sp->seedPos;
            sbpos.z=sbpos.z/2.; // half way up the mesocotyl
            numberOfRootCrowns = ceil((maxT-sp->firstSB)/sp->delayRC); // maximal number of root crowns
            double delay = sp->firstSB;
            for (int i=0; i<numberOfRootCrowns; i++) {
                Organ* shootborne0 = createRoot(plant, shootborneType, iheading ,delay, nullptr, 0, 0);
                // TODO fix the initial radial heading
                shootborne0->addNode(sbpos,delay);
                this->addChild(shootborne0);
                delay += sp->delaySB;
                for (int j=1; j<sp->nC; j++) {
                    Organ* shootborne = createRoot(plant, shootborneType, iheading ,delay, nullptr, 0, 0);
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
        Organ* mainstem = createStem(plant, mainStemType, isheading, 0, this ,0., 0.); // main stem has subtype 1
        mainstem->addNode(sp->seedPos, 0);
        children.push_back(mainstem);
        // Optional tillers
        if (sp->maxTil>0) {
            if (plant->getOrganRandomParameter(Organism::ot_stem, tillerType)->subType<1) { // if the type is not defined, copy tap root
                std::cout << "Tiller stem type #" << tillerType << " was not defined, using main stem parameters instead, ";
                std::cout << "default maxT = " << sp->maxTil << "\n";
                OrganRandomParameter* tillParam = plant->getOrganRandomParameter(Organism::ot_stem, 1)->copy(plant);
                tillParam->subType = basalType;
                plant->setOrganRandomParameter(tillParam);

            } else{
                int maxTi = sp->maxTil;
                if (sp->delayB>0) {
                    maxTi = std::min(maxTi,int(std::ceil((maxT-sp->firstB)/sp->delayB))); // maximal for simtime maxT
                }
                std::cout << "maxT = " << sp->maxTil << "\n";
                double delay = sp->firstB;
                for (int i=0; i<maxTi; i++) {
                    Organ* tiller = createStem(plant, tillerType, isheading, delay, this ,0., 0.);
                    tiller->addNode(sp->seedPos,0);
                    children.push_back(tiller);
                    delay += sp->delayB;
                }
            }
        }
    }
}

/**
 * Creates a shallow copy of the seeds child organs.
 * Use: the rootsystem manages base root itself, and just uses Seed class for initialization
 */
std::vector<Organ*> Seed::copyBaseOrgans()
{
    std::vector<Organ*> organs;
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
    auto orp = plant->getOrganRandomParameter(organtype);
    for (auto& o:orp) {
        if (o->name == str) {
            return o->subType;
        }
    }
    return -1;
}

/**
 * Quick info about the object for debugging (TODO)
 */
std::string Seed::toString() const
{
    std::stringstream str;
    str << "Seed #"<< id <<": type "<< param_->subType << ", length: "<< length << ", age: " << age;
    return str.str();
}

/**
 *
 */
Organ* Seed::createRoot(Organism* plant, int type, Vector3d iheading, double delay, Organ* parent, double pbl, int pni)
{
    return new Root(plant, type, iheading, delay, parent, pbl, pni);
}

/**
 *
 */
Organ* Seed::createStem(Organism* plant, int type, Vector3d iheading, double delay, Organ* parent, double pbl, int pni) // overwrite if you want to change the types
{
    return new Stem(plant, type, iheading, delay, parent, pbl, pni);
}

} // namespace CPlantBox
