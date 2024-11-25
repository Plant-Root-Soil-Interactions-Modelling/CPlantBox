#include "MycorrhizalRoot.h"

namespace CPlantBox {
    // Konstruktor 1
    // Konstruktor 2

    std::shared_ptr<Organ> MycorrhizalRoot::copy(std::shared_ptr<Organism> rs) {
        auto r = std::make_shared<MycorrhizalRoot>(*this);
        r->parent = std::weak_ptr<Organ>();
        r->plant = rs;
        r->param_ = std::make_shared<MycorrhizalRootSpecificParameter>(*param());
        for (size_t i=0; i< children.size(); i++) {
        r->children[i] = children[i]->copy(rs); // copy laterals
        r->children[i]->setParent(r);
        }
        return r;
    }
    //simulate
    // organType ?
    // simualte
    double MycorrhizalRoot::getParameter(std::string name) const {
        if (name=="vi") {return param()->vi;}
        return Root::getParameter(name);
    }
    //toString
    //getMycRootRandomParameter
    //Mycparam
}