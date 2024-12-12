#include "MycorrhizalRoot.h"
#include "Root.h"
#include "Organ.h"

namespace CPlantBox {

void MycorrhizalRoot::simulate(double dt, bool verbose)
{
    // std::cout << "\nstart" << getId() <<  std::flush;
    firstCall = true;
    moved = false;
    oldNumberOfNodes = nodes.size();

    const MycorrhizalRootSpecificParameter& p = *param(); // rename

    if (alive) { // dead roots wont grow

        // increase age
        if (age+dt>p.rlt) { // root life time
            dt=p.rlt-age; // remaining life span
            alive = false; // this root is dead
        }
        age+=dt;

        // probabilistic branching model
        if ((age>0) && (age-dt<=0)) { // the root emerges in this time step
            double P = getRootRandomParameter()->f_sbp->getValue(nodes.back(),shared_from_this());
            if (P<1.) { // P==1 means the lateral emerges with probability 1 (default case)
                double p = 1.-std::pow((1.-P), dt); //probability of emergence in this time step
                if (plant.lock()->rand()>p) { // not rand()<p
                    age -= dt; // the root does not emerge in this time step
                }
            }
        }

        if (age>0) { // unborn  roots have no children

            // children first (lateral roots grow even if base root is inactive)
            for (auto l:children) {
                l->simulate(dt,verbose);
            }


            if (active) {

                // length increment
                double age_ = calcAge(length); // root age as if grown unimpeded (lower than real age)
                double dt_; // time step
                if (age<dt) { // the root emerged in this time step, adjust time step
                    dt_= age;
                } else {
                    dt_=dt;
                }

                double targetlength = calcLength(age_+dt_)+ this->epsilonDx;

                double e = targetlength-length; // unimpeded elongation in time step dt
                double scale = getRootRandomParameter()->f_se->getValue(nodes.back(), shared_from_this());
                double dl = std::max(scale*e, 0.);//  length increment = calculated length + increment from last time step too small to be added
                length = getLength();
                this->epsilonDx = 0.; // now it is "spent" on targetlength (no need for -this->epsilonDx in the following)

                // create geometry
                if (p.laterals ) { // root has children
                    /* basal zone */
                    if ((dl>0)&&(length<p.lb)) { // length is the current length of the root
                        if (length+dl<=p.lb) {
                            createSegments(dl,dt_,verbose);
                            length+=dl; // - this->epsilonDx;
                            dl=0;
                        } else {
                            double ddx = p.lb-length;
                            createSegments(ddx,dt_,verbose);
                            dl-=ddx; // ddx already has been created
                            length=p.lb;
                            //							if(this->epsilonDx != 0){//this sould not happen as p.lb was redefined in rootparameter::realize to avoid this
                            //								throw std::runtime_error("Root::simulate: p.lb - length < dxMin");
                            //							} // this could happen, if the tip ends in this section
                        }
                    }

                    /* branching zone */
                    if ((dl>0)&&(length>=p.lb)) {
                        double s = p.lb; // summed length
                        for (size_t i=0; ((i<p.ln.size()) && (dl > 0)); i++) {
                            s+=p.ln.at(i);
                            if (length<=s) {//need "<=" instead of "<" => in some cases ln.at(i) == 0 when adapting ln to dxMin (@see rootrandomparameter::realize())

                                if (i==created_linking_node) { // new lateral
                                    createLateral(dt_, verbose);
                                }

                                if(length < s)//because with former check we have (length<=s)
                                {
                                    if (length+dl<=s) { // finish within inter-lateral distance i
                                        createSegments(dl,dt_,verbose);
                                        length+=dl; //- this->epsilonDx;
                                        dl=0;
                                    } else { // grow over inter-lateral distance i
                                        double ddx = s-length;
                                        createSegments(ddx,dt_,verbose);
                                        dl-=ddx;
                                        length=s;
                                        //									if(this->epsilonDx != 0){//this sould not happen as p.lb was redefined in rootparameter::realize to avoid this
                                        //										throw std::runtime_error( "Root::simulate: p.ln.at(i) - length < dxMin");
                                        //									} // this could happen, if the tip ends in this section
                                    }

                                }
                            }
                        }

                        if ((p.ln.size()==created_linking_node)&& (getLength(true)-s>-1e-9)){
                            createLateral(dt_, verbose);
                        }
                    }

                    /* apical zone */
                    if (dl>0) {
                        createSegments(dl,dt_,verbose);
                        length+=dl; // - this->epsilonDx;
                    }
                } else { // no laterals

                    if (dl>0) {
                        createSegments(dl,dt_,verbose);
                        length+=dl; //- this->epsilonDx;
                    }
                } // if lateralgetLengths
            } // if active
            active = getLength(false)<=(p.getK()*(1 - 1e-11)); // become inactive, if final length is nearly reached
        }
    } // if alive
    // std::cout << "end" << getId() << "\n" << std::flush;
}
    
std::shared_ptr<const MycorrhizalRootSpecificParameter> MycorrhizalRoot::param() const
{
    return std::static_pointer_cast<const MycorrhizalRootSpecificParameter>(param_);
}

std::shared_ptr<MycorrhizalRootRandomParameter> MycorrhizalRoot::getRootRandomParameter() const
{
    return std::static_pointer_cast<MycorrhizalRootRandomParameter>(plant.lock()->getOrganRandomParameter(Organism::ot_root, param_->subType));
}

double MycorrhizalRoot::getParameter(std::string name) const {
    // if (name == "infected") {return this->infected;}
    return Root::getParameter(name);
}

}