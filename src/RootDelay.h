// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef ROOTDELAY_H_
#define ROOTDELAY_H_

#include "Root.h"

namespace CPlantBox {

class RootState;

/**
 * RootDelay
 *
 * Laterals emerge after a specific prefdefined delay
 * (in contrast to emerging after the apical zone reaches a predefined length)
 */
class RootDelay :public Root
{
public:

    using Root::Root;
    std::shared_ptr<Organ> copy(std::shared_ptr<Organism> rs) override;  ///< deep copies the root tree
    std::string toString() const override;

	protected:
    void createLateral(double dt, bool silence) override; ///< creates a new lateral based on a delay

};

} // end namespace CPlantBox

#endif /* ROOTDELAY_H_ */
