// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "Seed.h"

#include "Root.h"

namespace CRootBox {

/**
 *
 */
void Seed::initialize()
{
    /*
     * Create base roots
     */
    const double maxT = 365.; // maximal simulation time
    const SeedSpecificParameter* rs = this->param(); // rename
    Vector3d iheading(0,0,-1);

    // Taproot
    Root* taproot = new Root(plant, 1, iheading ,0, nullptr, 0, 0); // tap root has root type 1
    taproot->addNode(rs->seedPos,0);
    this->addChild(taproot);

    auto pmap = plant->getOrganRandomParameter(Organism::ot_root);
    // Basal roots
    int bt = getParamSubType(Organism::ot_root, "basal");
    if (bt>0) {
        basalType = bt;
    } // otherwise stick with default
    if (rs->maxB>0) {
        try {
            plant->getOrganRandomParameter(Organism::ot_root, basalType); // if the type is not defined an exception is thrown
        } catch (...) {
            std::cout << "Seed::initialize: Basal root type #" << basalType << " was not defined, using tap root parameters instead\n" << std::flush;
            RootRandomParameter* brtp = (RootRandomParameter*)plant->getOrganRandomParameter(Organism::ot_root, 1)->copy(plant);
            brtp->subType = basalType;
            plant->setOrganRandomParameter(brtp);
        }
        int maxB = rs->maxB;
        if (rs->delayB>0) { // limit if possible
            maxB = std::min(maxB,int(ceil((maxT-rs->firstB)/rs->delayB))); // maximal for simtime maxT
        }
        double delay = rs->firstB;
        for (int i=0; i<maxB; i++) {
            Root* basalroot = new Root(plant, basalType, iheading, delay, nullptr, 0, 0);
            basalroot->addNode(taproot->getNode(0), taproot->getNodeId(0), delay);
            this->addChild(basalroot);
            delay += rs->delayB;
        }
    }
    // Shoot borne roots
    int st = getParamSubType(Organism::ot_root, "shootborne");
    if (st>0) {
        shootborneType = st;
    } // otherwise stick with default
    if ((rs->nC>0) && (rs->firstSB+rs->delaySB<maxT)) { // only if there are any shootborne roots
        try {
            plant->getOrganRandomParameter(Organism::ot_root, shootborneType); // if the type is not defined an exception is thrown
        } catch (...) {
            std::cout << "Seed::initialize:Shootborne root type #" << shootborneType << " was not defined, using tap root parameters instead\n";
            RootRandomParameter* srtp =  (RootRandomParameter*)plant->getOrganRandomParameter(Organism::ot_root, 1)->copy(plant);
            srtp->subType = shootborneType;
            plant->setOrganRandomParameter(srtp);
        }
        Vector3d sbpos = rs->seedPos;
        sbpos.z=sbpos.z/2.; // half way up the mesocotyl
        numberOfRootCrowns = ceil((maxT-rs->firstSB)/rs->delayRC); // maximal number of root crowns
        double delay = rs->firstSB;
        for (int i=0; i<numberOfRootCrowns; i++) {
            Root* shootborne0 = new Root(plant, shootborneType, iheading ,delay, nullptr, 0, 0);
            // TODO fix the initial radial heading
            shootborne0->addNode(sbpos,delay);
            this->addChild(shootborne0);
            delay += rs->delaySB;
            for (int j=1; j<rs->nC; j++) {
                Root* shootborne = new Root(plant, shootborneType, iheading ,delay, nullptr, 0, 0);
                // TODO fix the initial radial heading
                shootborne->addNode(shootborne0->getNode(0), shootborne0->getNodeId(0),delay);
                this->addChild(shootborne);
                delay += rs->delaySB;
            }
            sbpos.z+=rs->nz;  // move up, for next root crown
            delay = rs->firstSB + i*rs->delayRC; // reset age
        }
    } else {
        numberOfRootCrowns = 0;
    }

    /*
     * Create Stem
     */
//    if (Plant::noParamFile[2] == 1) {
//        std::cout<<"no stem parameter, maybe the XML based parameter file is directly converted from the old rparam file"<<std::endl;
//    } else {
//        //	Vector3d isheading(0,0,1);//Initial Stem heading
//        Vector3d isheading(0, 0, 1);//Initial Stem heading
//        Stem* mainstem = new Stem(plant, this, 1, 0., isheading, 0., 0.); // tap root has subtype 1
//        mainstem->addNode(sparam->seedPos, 0);
//        children.push_back(mainstem);
//
//
//        if (sparam->maxTi>0) {
//            if (plant->getParameter(Organism::ot_stem, tillerType)->subType<1) { // if the type is not defined, copy tap root
//                std::cout << "Basal root type #" << basalType << " was not defined, using tap root parameters instead\n";
//                StemTypeParameter* tillParam = (StemTypeParameter*)plant->getParameter(Organism::ot_stem, 1);
//                //				StemTypeParameter* titp = new StemTypeParameter(*tillParam);
//                //				titp->subType = tillerType;
//                //				plant->setParameter(titp);
//                std::cout << "default maxT type is main stem = " << sparam->maxTi << "\n";
//
//            } else{
//                int maxTi = sparam->maxTi;
//                if (sparam->delayB>0) {
//                    maxTi = std::min(maxTi,int(std::ceil((maxT-sparam->firstB)/sparam->delayB))); // maximal for simtime maxT
//                }
//                std::cout << "maxT = " << sparam->maxTi << "\n";
//                double delay = sparam->firstB;
//                StemTypeParameter* tillParam = (StemTypeParameter*)plant->getParameter(Organism::ot_stem, 4);
//                for (int i=0; i<maxTi; i++) {
//                    Stem* tiller = new Stem(plant, this, 4, delay, isheading ,0., 0.);
//                    tiller->addNode(sparam->seedPos,0);
//                    children.push_back(tiller);
//
//                    std::cout << "new maxT type is main stem = " << sparam->maxTi << "\n";
//                }
//
//            }
//
//            //	Stem* tiller1 = new Stem(plant, this, 1, 2, isheading ,0., 0.); // tap root has subtype 1
//            //	tiller1->addNode(sparam->seedPos,0);
//            //	children.push_back(tiller1);
//            //
//            //	Stem* tiller2 = new Stem(plant, this, 1, 4 , isheading ,0., 0.); // tap root has subtype 1
//            //	tiller2->addNode(sparam->seedPos,0);
//            //	children.push_back(tiller2);
//            //	Stem* tiller3 = new Stem(plant, this, 1, 5, isheading ,0., 0.); // tap root has subtype 1
//            //	tiller3->addNode(sparam->seedPos,0);
//            //	children.push_back(tiller3);
//        }
//
//    }
}

/**
 * Creates a shallow copy of the seeds child organs
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

} // namespace CPlantBox
