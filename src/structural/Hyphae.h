// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef HYPHAE_H_
#define HYPHAE_H_

#include "mymath.h"
#include "Organ.h"
#include "Organism.h"
#include "hyphaeparameter.h"

namespace CPlantBox {

/**
 * Hyphae
 */
class Hyphae :public Organ
{
public:

    Hyphae(int id, std::shared_ptr<const OrganSpecificParameter> param, bool alive, bool active, double age, double length,
        Vector3d partialIHeading_, int pni, bool moved= false, int oldNON = 0); // ///< creates everything from scratch

    Hyphae(std::shared_ptr<Organism> rs, int type, double delay, std::shared_ptr<Organ> parent, int pni); ///< used within simulation

    virtual ~Hyphae() { }; ///< no need to do anything, children are deleted in ~Organ()

    std::shared_ptr<Organ> copy(std::shared_ptr<Organism> rs) override;  ///< deep copies the root tree

    int organType() const override { return Organism::ot_hyphae; }; ///< returns the organs type

    void simulate(double dt, bool silence = false) override; ///< root growth for a time span of @param dt

    double getParameter(std::string name) const override;

    void createLateral(double pni); ///< creates a lateral hyphae

    std::string toString() const override;

    // void makeanastomosis(std::shared_ptr<Hyphae> a, std::shared_ptr<Hyphae> b); ///< creates an anastomosis with another hyphae

//    /* From analytical equations */
    double calcLength(double age); ///< analytical length of the root
    double calcAge(double length) const; ///< analytical age of the root

    /* Abbreviations */
    std::shared_ptr<HyphaeRandomParameter> getHyphaeRandomParameter() const;  ///< root type parameter of this root
    std::shared_ptr<const HyphaeSpecificParameter> param() const; ///< root parameter


};

} // end namespace CPlantBox

#endif /* HYPHAE_H_ */
