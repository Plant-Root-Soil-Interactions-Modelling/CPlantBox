#include "Leaf.h"

namespace CPlantBox {

std::vector<int> Leaf::leafphytomerID = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

/**
 * Constructs a root from given data.
 * The organ tree must be created, @see Organ::setPlant, Organ::setParent, Organ::addChild
 * Organ geometry must be created, @see Organ::addNode, ensure that this->getNodeId(0) == parent->getNodeId(pni)
 *
 * @param id        the organ's unique id (@see Organ::getId)
 * @param param     the organs parameters set, ownership transfers to the organ
 * @param alive     indicates if the organ is alive (@see Organ::isAlive)
 * @param active    indicates if the organ is active (@see Organ::isActive)
 * @param age       the current age of the organ (@see Organ::getAge)
 * @param length    the current length of the organ (@see Organ::getLength)
 * @param iheading  the initial heading of this root
 * @param pbl       base length of the parent root, where this root emerges
 * @param pni       local node index, where this root emerges
 * @param moved     indicates if nodes were moved in the previous time step (default = false)
 * @param oldNON    the number of nodes of the previous time step (default = 0)
 */
Leaf::Leaf(int id, const std::shared_ptr<const OrganSpecificParameter> param, bool alive, bool active, double age, double length,
    Vector3d iheading, double pbl, int pni, bool moved, int oldNON)
:Organ(id, param, alive, active, age, length, moved,  oldNON ), parentBaseLength(pbl), parentNI(pni)
{ }

/**
 * Constructor
 * This is a Copy Paste of the Root.cpp but it works independently, it has its own parameter file (in .lParam file) tropism, f_gf function, txt and vtp writing syleaf.
 * All of those can be modified to fit the real f_gf of the Plant.
 *
 * Typically called by the Plant::Plant(), or Leaf::createNewLeaf().
 * For leaf the initial node and node emergence time (netime) must be set from outside
 *
 * @param plant 		points to the plant
 * @param parent 		points to the parent organ
 * @param subtype		sub type of the leaf
 * @param delay 		delay after which the organ starts to develop (days)
 * @param rheading		relative heading (within parent organ)
 * @param pni			parent node index
 * @param pbl			parent base length
 */
Leaf::Leaf(std::shared_ptr<Organism> plant, int type, Vector3d iheading, double delay,  std::shared_ptr<Organ> parent, double pbl, int pni)
        :Organ(plant, parent,Organism::ot_leaf, type, delay), parentBaseLength(pbl), parentNI(pni)
{
    assert(parent!=nullptr && "Stem::Stem parent must be set");
    //  std::cout << "Leaf pni = "<< pni<< std::endl;
    //  std::cout << "Organism* plant ="<< plant <<" "<< parent<<std::endl;
    //  std::cout <<", "<<(LeafParameter*) param<< "\n";
    // std::cout <<"subtype ="<<param()->subType <<"getleafphytomerID =" <<getleafphytomerID(param()->subType)<< "\n";

    addleafphytomerID(param()->subType);
    double beta = getleafphytomerID(param()->subType)*M_PI*getLeafRandomParameter()->RotBeta
        + M_PI*plant->rand()*getLeafRandomParameter()->BetaDev ;  //+ ; //2 * M_PI*plant->rand(); // initial rotation
    Matrix3d ons = Matrix3d::ons(iheading);
    //	if (getLeafRandomParameter()->InitBeta >0 && getleafphytomerID(param()->subType)==0 ){
    beta = beta + getLeafRandomParameter()->InitBeta;
    //	}

    if (getLeafRandomParameter()->InitBeta >0 && getLeafRandomParameter()->subType==2 && getLeafRandomParameter()->lnf==5 && getleafphytomerID(2)%4==2 )
    {beta = beta + getLeafRandomParameter()->InitBeta*M_PI;}
    else if (getLeafRandomParameter()->InitBeta >0 && getLeafRandomParameter()->subType==2 && getLeafRandomParameter()->lnf==5 && getleafphytomerID(2)%4==3 )
    {beta = beta + getLeafRandomParameter()->InitBeta*M_PI + M_PI;}
    //ons.times(Matrix3d::rotX(beta));

    double theta = M_PI*param()->theta;
    if (parent->organType()!=Organism::ot_seed) { // scale if not a base root
        double scale = getLeafRandomParameter()->f_sa->getValue(parent->getNode(pni), parent);
        theta *= scale;
    }
    //ons.times(Matrix3d::rotZ(theta));
    this->iHeading = ons.times(Vector3d::rotAB(theta,beta)); // new initial heading

    // initial node
    if (parent->organType()!=Organism::ot_seed) { // the first node of the base roots must be created in RootSystem::initialize()
        // assert(pni+1 == parent->getNumberOfNodes() && "at object creation always at last node"); // ?????
        addNode(parent->getNode(pni), parent->getNodeId(pni), parent->getNodeCT(pni)+delay);
    }
}

/**
 * Deep copies the organ into the new plant @param rs.
 * All laterals are deep copied, plant and parent pointers are updated.
 *
 * @param plant     the plant the copied organ will be part of
 */
std::shared_ptr<Organ> Leaf::copy(std::shared_ptr<Organism> p)
{
    auto l = std::make_shared<Leaf>(*this); // shallow copy
    l->parent = std::weak_ptr<Organ>();
    l->plant = p;
    l->param_ = std::make_shared<LeafSpecificParameter>(*param()); // copy parameters
    for (size_t i=0; i< children.size(); i++) {
        l->children[i] = children[i]->copy(p); // copy laterals
        l->children[i]->setParent(shared_from_this());
    }
    return l;
}

/**
 * Simulates f_gf of this leaf for a time span dt
 *
 * @param dt       time step [day]
 * @param verbose  indicates if status messages are written to the console (cout) (default = false)
 */
void Leaf::simulate(double dt, bool verbose)
{
    firstCall = true;
    moved = false;
    oldNumberOfNodes = nodes.size();

    const LeafSpecificParameter& p = *param(); // rename

    if (alive) { // dead roots wont grow

        // increase age
        if (age+dt>p.rlt) { // root life time
            dt=p.rlt-age; // remaining life span
            alive = false; // this root is dead
        }
        age+=dt;

        // probabilistic branching model (todo test)
        if ((age>0) && (age-dt<=0)) { // the root emerges in this time step
            double P = getLeafRandomParameter()->f_sbp->getValue(nodes.back(),shared_from_this());
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

                double targetlength = calcLength(age_+dt_);
                double e = targetlength-length; // unimpeded elongation in time step dt
                double scale = getLeafRandomParameter()->f_sa->getValue(nodes.back(),shared_from_this());
                double dl = std::max(scale*e, 0.); // length increment

                // create geometry
                if (p.ln.size()>0) { // leaf has laterals
                    /* basal zone */
                    if ((dl>0)&&(length<p.lb)) { // length is the current length of the root
                        if (length+dl<=p.lb) {
                            createSegments(dl,verbose);
                            length+=dl;
                            dl=0;
                        } else {
                            double ddx = p.lb-length;
                            createSegments(ddx,verbose);
                            dl-=ddx; // ddx already has been created
                            length=p.lb;
                        }
                    }
                    /* branching zone */
                    if ((dl>0)&&(length>=p.lb)) {
                        double s = p.lb; // summed length
                        for (size_t i=0; ((i<p.ln.size()) && (dl>0)); i++) {
                            s+=p.ln.at(i);
                            if (length<s) {
                                if (i==children.size()) { // new lateral
                                    createLateral(verbose);
                                }
                                if (length+dl<=s) { // finish within inter-lateral distance i
                                    createSegments(dl,verbose);
                                    length+=dl;
                                    dl=0;
                                } else { // grow over inter-lateral distance i
                                    double ddx = s-length;
                                    createSegments(ddx,verbose);
                                    dl-=ddx;
                                    length=s;
                                }
                            }
                        }
                        if (p.ln.size()==children.size()) { // new lateral (the last one)
                            createLateral(verbose);
                        }
                    }
                    /* apical zone */
                    if (dl>0) {
                        createSegments(dl,verbose);
                        length+=dl;
                    }
                } else { // no laterals
                    if (dl>0) {
                        createSegments(dl,verbose);
                        length+=dl;
                    }
                } // if lateralgetLengths
            } // if active
            active = length<(p.getK()-dx()/10); // become inactive, if final length is nearly reached
        }
    } // if alive

}

/**
 *
 */
double Leaf::getParameter(std::string name) const {
    double r = Organ::getParameter(name);
    if (name=="basal zone") { r = param()->lb; }
    if (name=="apical zone") { r = param()->la; }
    if (name=="initial f_gf rate") { r = param()->r; }
    if (name=="radius") { r = param()->a; }
    if (name=="insertion angle") { r = param()->theta; }
    if (name=="root life time") { r = param()->rlt; }
    //	if (name=="mean internodal distance") {
    //		r = std::accumulate(param()->ln.begin(), param()->ln.end(), 0.0) / param()->ln.size();
    //	}
    //	if (name=="sd internodal distance") {
    //		const std::vector<double>& v_ = param()->ln;
    //		double mean = std::accumulate(v_.begin(), v_.end(), 0.0) / v_.size();
    //		double sq_sum = std::inner_product(v_.begin(), v_.end(), v_.begin(), 0.0);
    //		r = std::sqrt(sq_sum / v_.size() - mean * mean);
    //	}
    //	if (name=="surface") { r = param()->a*param()->a*M_PI*length; }
    //	if (name=="number of branches") { r = param()->ln.size() +1; }
    return r;
}

/**
 * Analytical creation (=emergence) time of a node at a length along the leaf
 *
 * @param length   length of the leaf [cm]
 */
double Leaf::calcCreationTime(double length)
{
    assert(length >= 0 && "Leaf::getCreationTime() negative length");
    double leafage = calcAge(length);
    leafage = std::min(leafage, age);
    assert(leafage >= 0 && "Leaf::getCreationTime() negative root age");
    return leafage+nodeCTs[0];
}

/**
 * Analytical length of the leaf at a given age
 *
 * @param age          age of the leaf [day]
 */
double Leaf::calcLength(double age)
{
    //    std::string organ_name = std::string(getLeafRandomParameter()->organName);
    //    //std::cout<<"organName is "<<name()<<"\n";
    //    if (name()  == "maize1eaf"){
    assert(age>=0);
    std::cout<<"organName is ? \n";
    //return getLeafRandomParameter()->f_gf->LeafgetLength(age,getLeafRandomParameter()->r,getLeafRandomParameter()->getK(),this);
    return getLeafRandomParameter()->f_gf->getLength(age,getLeafRandomParameter()->r,getleafphytomerID(getLeafRandomParameter()->subType)*3,shared_from_this());
    //	}else {
    //assert(age>=0);
    //	    return getLeafRandomParameter()->f_gf->LeafgetLength(age,getLeafRandomParameter()->r,getLeafRandomParameter()->getK(),this);
    //	    }
}

/**
 * Analytical age of the leaf at a given length
 *
 * @param length   length of the leaf [cm]
 */
double Leaf::calcAge(double length)
{
    //    std::string organ_name = std::string(getLeafRandomParameter()->organName);
    //    //std::cout<<getLeafRandomParameter()->name<< "\n";
    //     if ( name() == "maize1eaf"){
    assert(length>=0);
    //std::cout<<"length subtype is"<<getLeafRandomParameter()->subType<<"\n";
    return getLeafRandomParameter()->f_gf->getAge(length,getLeafRandomParameter()->r,getLeafRandomParameter()->getK(),shared_from_this());
    //        return getLeafRandomParameter()->f_gf->LeafgetAge(length,getLeafRandomParameter()->r,getleafphytomerID(getLeafRandomParameter()->subType)*3,this);
    //        }else {
    //assert(age>=0);
    //	    return getLeafRandomParameter()->f_gf->LeafgetAge(length,getLeafRandomParameter()->r,getLeafRandomParameter()->getK(),this);
    //	    }
}

/**
 *
 */
double Leaf::dx() const
{
    return getLeafRandomParameter()->dx;
}

/**
 * Creates a new lateral by calling Leaf::createNewleaf().
 *
 * Overwrite this method to implement more specialized leaf classes.
 */
void Leaf::createLateral(bool silence)
{

    int lt = getLeafRandomParameter()->getLateralType(getNode(nodes.size()-1));

    if (lt>0) {

        if (getLeafRandomParameter()->lnf==2&& lt>0) {
            double ageLN = this->calcAge(length); // age of Leaf when lateral node is created
            double ageLG = this->calcAge(length+param()->la); // age of the Leaf, when the lateral starts growing (i.e when the apical zone is developed)
            double delay = ageLG-ageLN; // time the lateral has to wait
            Vector3d h = heading(); // current heading
            auto lateral = std::make_shared<Leaf>(plant.lock(), lt, h, delay, shared_from_this(), nodes.size() - 1, length);
            //		lateral->setRelativeOrigin(nodes.back());
            children.push_back(lateral);
            lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
            auto lateral2 = std::make_shared<Leaf>(plant.lock(), lt, h, delay, shared_from_this(), nodes.size() - 1, length);
            //		lateral2->setRelativeOrigin(nodes.back());
            children.push_back(lateral2);
            lateral2->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
        }
        else if (getLeafRandomParameter()->lnf==3&& lt>0) { //ln equal and both side leaf
            double ageLN = this->calcAge(length); // age of Leaf when lateral node is created
            double ageLG = this->calcAge(length+param()->la); // age of the Leaf, when the lateral starts growing (i.e when the apical zone is developed)
            double delay = ageLG-ageLN; // time the lateral has to wait
            Vector3d h = heading(); // current heading
            auto lateral = std::make_shared<Leaf>(plant.lock(),  lt, h, delay, shared_from_this(), nodes.size() - 1, length);
            //		lateral->setRelativeOrigin(nodes.back());
            children.push_back(lateral);
            lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
            auto lateral2 = std::make_shared<Leaf>(plant.lock(),  lt, h, delay, shared_from_this(), nodes.size() - 1, length);
            //		lateral2->setRelativeOrigin(nodes.back());
            children.push_back(lateral2);
            lateral2->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
        }
        else if (getLeafRandomParameter()->lnf==4 && lt>0) {//ln exponential decreasing and one side leaf
            double ageLN = this->calcAge(length); // age of Leaf when lateral node is created
            double ageLG = this->calcAge(length+param()->la); // age of the Leaf, when the lateral starts growing (i.e when the apical zone is developed)
            double delay = ageLG-ageLN; // time the lateral has to wait
            Vector3d h = heading(); // current heading
            auto lateral = std::make_shared<Leaf>(plant.lock(),  lt, h, delay, shared_from_this(), nodes.size() - 1, length);
            //		lateral->setRelativeOrigin(nodes.back());
            children.push_back(lateral);
            lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)

        } else if (getLeafRandomParameter()->lnf==5&& lt>0) { //ln exponential decreasing and both side leaf
            double ageLN = this->calcAge(length); // age of Leaf when lateral node is created
            double ageLG = this->calcAge(length+param()->la); // age of the Leaf, when the lateral starts growing (i.e when the apical zone is developed)
            double delay = ageLG-ageLN; // time the lateral has to wait
            Vector3d h = heading(); // current heading
            auto lateral = std::make_shared<Leaf>(plant.lock(), lt,  h, delay,  shared_from_this(), nodes.size() - 1, length);
            //		lateral->setRelativeOrigin(nodes.back());
            children.push_back(lateral);
            lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)

            addleafphytomerID(getLeafRandomParameter()->subType);

            auto lateral2 = std::make_shared<Leaf>(plant.lock(), lt, h, delay,  shared_from_this(), nodes.size() - 1, length);
            //		lateral2->setRelativeOrigin(nodes.back());
            children.push_back(lateral2);
            lateral2->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
        } else if (lt>0) {
            double ageLN = this->calcAge(length); // age of Leaf when lateral node is created
            double ageLG = this->calcAge(length+param()->la); // age of the Leaf, when the lateral starts growing (i.e when the apical zone is developed)
            double delay = ageLG-ageLN; // time the lateral has to wait
            Vector3d h = heading(); // current heading

            auto lateral = std::make_shared<Leaf>(plant.lock(), lt, h, delay, shared_from_this(), nodes.size() - 1, length);
            //		lateral->setRelativeOrigin(nodes.back());
            children.push_back(lateral);
            lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
        }else {


            double ageLN = this->calcAge(length); // age of leaf when lateral node is created
            double ageLG = this->calcAge(length+param()->la); // age of the leaf, when the lateral starts growing (i.e when the apical zone is developed)
            double delay = ageLG-ageLN; // time the lateral has to wait
            Vector3d h = heading(); // current heading
            auto lateral = std::make_shared<Leaf>(plant.lock(), lt, h, delay, shared_from_this(), nodes.size()-1, length);
            //		lateral->setRelativeOrigin(nodes.back());
            children.push_back(lateral);
            lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
        }
    }
}

/**
 * @return Current root heading
 */
Vector3d Leaf::heading() const
{
    if (nodes.size()>1) {
        auto h = nodes.back().minus(nodes.at(nodes.size()-2)); // a->b = b-a
        h.normalize();
        return h;
    } else {
        return iHeading;
    }
}

/**
 * Returns the increment of the next segments
 *
 *  @param p       coordinates of previous node
 *  @param sdx     length of next segment [cm]
 *  @return        the vector representing the increment
 */
Vector3d Leaf::getIncrement(const Vector3d& p, double sdx)
{
    Vector3d h = heading();
    Matrix3d ons = Matrix3d::ons(h);
    Vector2d ab = getLeafRandomParameter()->f_tf->getHeading(p, ons, sdx,shared_from_this());
    Vector3d sv = ons.times(Vector3d::rotAB(ab.x,ab.y));
    return sv.times(sdx);
}

/**
 *  Creates nodes and node emergence times for a length l
 *
 *  Checks that each new segments length is <= dx but >= smallDx
 *
 *  @param l        total length of the segments that are created [cm]
 *  @param verbose  turns console output on or off
 */
void Leaf::createSegments(double l, bool verbose)
{
    if (l==0) {
        std::cout << "Root::createSegments: zero length encountered \n";
        return;
    }
    if (l<0) {
        std::cout << "Root::createSegments: negative length encountered \n";
    }

    // shift first node to axial resolution
    double shiftl = 0; // length produced by shift
    int nn = nodes.size();
    if (firstCall) { // first call of createSegments (in Root::simulate)
        firstCall = false;
        if ((nn>1) && (children.empty() || (nn-1 != std::static_pointer_cast<Leaf>(children.back())->parentNI)) ) { // don't move a child base node
            Vector3d n2 = nodes[nn-2];
            Vector3d n1 = nodes[nn-1];
            double olddx = n1.minus(n2).length(); // length of last segment
            if (olddx<dx()*0.99) { // shift node instead of creating a new node
                shiftl = std::min(dx()-olddx, l);
                double sdx = olddx + shiftl; // length of new segment
                Vector3d newdxv = getIncrement(n2, sdx);
                nodes[nn-1] = Vector3d(n2.plus(newdxv));
                double et = this->calcCreationTime(length+shiftl);
                nodeCTs[nn-1] = et; // in case of impeded growth the node emergence time is not exact anymore, but might break down to temporal resolution
                moved = true;
                l -= shiftl;
                if (l<=0) { // ==0 should be enough
                    return;
                }
            } else {
                moved = false;
            }
        } else {
            moved = false;
        }
    }
    // create n+1 new nodes
    double sl = 0; // summed length of created segment
    int n = floor(l/dx());
    for (int i = 0; i < n + 1; i++) {

        double sdx; // segment length (<=dx)
        if (i<n) {  // normal case
            sdx = dx();
        } else { // last segment
            sdx = l-n*dx();
            if (sdx<smallDx) { // quit if l is too small
                if (verbose) {
                    std::cout << "skipped small segment ("<< sdx <<" < "<< smallDx << ") \n";
                }
                return;
            }
        }
        sl += sdx;
        Vector3d newdx = getIncrement(nodes.back(), sdx);
        Vector3d newnode = Vector3d(nodes.back().plus(newdx));
        double et = this->calcCreationTime(length+shiftl+sl);
        // in case of impeded growth the node emergence time is not exact anymore,
        // but might break down to temporal resolution
        addNode(newnode, et);
    }
}

/**
 * @return The RootTypeParameter from the plant
 */
std::shared_ptr<LeafRandomParameter> Leaf::getLeafRandomParameter() const
{
    return std::static_pointer_cast<LeafRandomParameter>(plant.lock()->getOrganRandomParameter(Organism::ot_root, param_->subType));
}

/**
 * @return Parameters of the specific root
 */
std::shared_ptr<const LeafSpecificParameter>  Leaf::param() const
{
    return std::static_pointer_cast<const LeafSpecificParameter>(param_);
}

/**
 * Quick info about the object for debugging
 * additionally, use getParam()->toString() and getOrganRandomParameter()->toString() to obtain all information.
 */
std::string Leaf::toString() const
{
    std::string str = Organ::toString();
    str.replace(0, 5, "Leaf");
    std::stringstream newstring;
    newstring << "; initial heading: " << iHeading.toString() << ", parent base length " << parentBaseLength << ", parent node index" << parentNI << ".";
    return str+newstring.str();
}

} // namespace CPlantBox
