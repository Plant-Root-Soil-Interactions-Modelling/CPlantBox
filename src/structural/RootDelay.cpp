// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "RootDelay.h"

namespace CPlantBox {

/**
 * Deep copies the organ into the new plant @param rs.
 * All laterals are deep copied, plant and parent pointers are updated.
 *
 * @param plant     the plant the copied organ will be part of
 */
std::shared_ptr<Organ> RootDelay::copy(std::shared_ptr<Organism> rs)
{
    auto r = std::make_shared<RootDelay>(*this); // shallow copy
    r->parent = std::weak_ptr<Organ>();
    r->plant = rs;
    r->param_ = std::make_shared<RootSpecificParameter>(*param()); // copy parameters
    for (size_t i=0; i< children.size(); i++) {
        r->children[i] = children[i]->copy(rs); // copy laterals
        r->children[i]->setParent(r);
    }
    return r;
}

/**
 * @return Quick info about the object for debugging
 * additionally, use param()->toString() and getOrganRandomParameter()->toString() to obtain all information.
 */
std::string RootDelay::toString() const
{
    std::string str = Organ::toString();
    str.replace(0, 6, "Delay");
    std::stringstream newstring;
    newstring << "; initial heading: " << getiHeading0().toString() << ", parent node index" << parentNI << ".";
    return str+newstring.str();
}

/**
 * Creates a new lateral root and passes time overhead
 *
 * Overwrite this method to implement more spezialized root classes.
 *
 * @param verbose   turns console output on or off
 */
void RootDelay::createLateral(double dt, bool verbose)
{
	// std::cout<< "create delayed root\n";
	auto rrp = getRootRandomParameter(); // rename
	for(int i = 0; i < rrp->successorST.size(); i++)
	{//go through each successor rule
		//found id
		bool applyHere = getApplyHere(i);
		if(applyHere)
		{
			int numlats = 1;//how many laterals? default = 1

			if (rrp->successorNo.size()>i) {
				numlats =  rrp->successorNo.at(i);
			}

			for(int nn = 0; nn < numlats; nn++) {
				const Vector3d& pos = Vector3d();
				int p_id = rrp->getLateralType(pos, i);
				if (p_id>=0) {
					int lt = rrp->successorST.at(i).at(p_id);
					double delay = std::max(rrp->ldelay + plant.lock()->randn()*rrp->ldelays, 0.);
					auto lateral = std::make_shared<RootDelay>(plant.lock(), lt,   delay,  shared_from_this(), nodes.size()-1);
					children.push_back(lateral);
					double ageLN = this->calcAge(length); // age of root when lateral node is created
					ageLN = std::max(ageLN, age-dt); // dt_*(1-dl/dl0) are ready
					lateral->simulate(age-ageLN,verbose); // pass time overhead (age we want to achieve minus current age)
				}
			}
		}
	}
	created_linking_node ++;
}

} // end namespace CPlantBox
