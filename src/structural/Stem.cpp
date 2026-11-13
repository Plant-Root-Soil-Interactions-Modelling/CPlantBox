#include "Stem.h"

#include "Leaf.h"
#include "Plant.h"
#include "Root.h"
#include "Seed.h"

#include <algorithm>

namespace CPlantBox {

/**
 * Constructs a root from given data.
 * The organ tree must be created, @see Organ::setPlant, Organ::setParent, Organ::addChild
 * Organ geometry must be created, @see Organ::addNode, ensure that this->getNodeId(0) == parent->getNodeId(pni)
 *
 * @param id        		the organ's unique id (@see Organ::getId)
 * @param param     		the organs parameters set, ownership transfers to the organ
 * @param alive     		indicates if the organ is alive (@see Organ::isAlive)
 * @param active    		indicates if the organ is active (@see Organ::isActive)
 * @param age       		the current age of the organ (@see Organ::getAge)
 * @param length    		the current length of the organ (@see Organ::getLength)
 * @param partialIHeading 	the initial partial heading of this root
 * @param pbl       		base length of the parent root, where this root emerges
 * @param pni       		local node index, where this root emerges
 * @deprecated moved		indicates if nodes were moved in the previous time step (default = false)
 * @param oldNON    		the number of nodes of the previous time step (default = 0)
 */
Stem::Stem(int id, std::shared_ptr<const OrganSpecificParameter> param, bool alive, bool active, double age, double length, Vector3d partialIHeading_, int pni, bool moved, int oldNON)
    : Organ(id, param, alive, active, age, length, partialIHeading_, pni, moved, oldNON) {}

/**
 * Constructor
 * This is a Copy Paste of the Root.cpp but it works independently, it has its own parameter file (in .stparam file) tropism, growth function, txt and vtp
 * writing system. All of those can be modified to fit the real growth of the Plant.
 *
 * Typically called by the Plant::Plant(), or Stem::createNewStem().
 * For stem the initial node and node emergence time (netime) must be set from outside
 *
 * @param plant 		points to the plant
 * @param parent 		points to the parent organ
 * @param subtype		sub type of the stem
 * @param delay 		delay after which the organ starts to develop (days)
 * @param rheading		relative heading (within parent organ)
 * @param pni			parent node index
 * @param pbl			parent base length
 */
Stem::Stem(std::shared_ptr<Organism> plant, int type, double delay, std::shared_ptr<Organ> parent, int pni) : Organ(plant, parent, Organism::ot_stem, type, delay, pni) {
    assert(parent != nullptr && "Stem::Stem parent must be set");
    auto p = this->param();
    int phytomerId = 0;
    if (parent->organType() == Organism::ot_stem) {
        std::static_pointer_cast<Stem>(parent)->addPhytomerId(p->subType);
        phytomerId = std::static_pointer_cast<Stem>(parent)->getphytomerId(p->subType);
    }
    double beta = phytomerId * M_PI * getStemRandomParameter()->rotBeta + M_PI * plant->rand() * getStemRandomParameter()->betaDev;
    beta = beta + getStemRandomParameter()->initBeta * M_PI;
    if (getStemRandomParameter()->initBeta > 0 && phytomerId) {
        beta = beta + getStemRandomParameter()->initBeta * M_PI;
    }
    double theta = p->theta;                        // M_PI*p->theta;
    if (parent->organType() != Organism::ot_seed) { // scale if not a base organ, to delete?
        double scale = getStemRandomParameter()->f_sa->getValue(parent->getNode(pni), parent);
        theta *= scale;
    }
    // used when computing actual heading, @see Stem::getIHeading
    this->partialIHeading = Vector3d::rotAB(theta, beta);
    if (parent->organType() != Organism::ot_seed) { // initial node
        // if lateral of stem, initial creation time:
        // time when stem reached end of basal zone (==CT of parent node of first lateral) + delay
        //  @see stem::createLateral
        double creationTime;
        if (parent->getNumberOfChildren() == 0) {
            creationTime = parent->getNodeCT(pni) + delay;
        } else {
            creationTime = parent->getChild(0)->getParameter("creationTime") + delay;
        }
        addNode(Vector3d(0., 0., 0.), parent->getNodeId(pni), creationTime);
    }
}

/**
 * Deep copies the organ into the new plant @param rs.
 * All laterals are deep copied, plant and parent pointers are updated.
 *
 * @param plant     the plant the copied organ will be part of
 */
std::shared_ptr<Organ> Stem::copy(std::shared_ptr<Organism> p) {
    auto s = std::make_shared<Stem>(*this); // shallow copy
    s->parent = std::weak_ptr<Organ>();
    s->plant = p;
    s->param_ = std::make_shared<StemSpecificParameter>(*param()); // copy parameters
    for (size_t i = 0; i < children.size(); i++) {
        s->children[i] = children[i]->copy(p); // copy laterals
        s->children[i]->setParent(s);
    }
    return s;
}

/**
 * Simulates growth of this stem for a time span dt
 *
 * @param dt       time step [day]
 * @param verbose  indicates if status messages are written to the console (cout) (default = false)
 */
void Stem::simulate(double dt, bool verbose) {
    const StemSpecificParameter &p = *param(); // rename
    firstCall = true;
    oldNumberOfNodes = nodes.size();

    if (alive) { // dead roots wont grow

        // increase age
        if (age + dt > p.rlt) { // root life time
            dt = p.rlt - age;   // remaining life span
            alive = false;      // this root is dead
        }
        age += dt;

        // we currently do not use a probabilistic branching model

        if (age > 0) { // unborn  roots have no children

            // children first (lateral roots grow even if base root is inactive)
            for (auto l : children) {
                l->simulate(dt, verbose);
            }

            if (active) {

                double dt_;     // time step
                if (age < dt) { // the root emerged in this time step, adjust time step
                    dt_ = age;
                } else {
                    dt_ = dt;
                }
                double age__ = age;         // length increment
                if (age > p.delayNGStart) { // simulation ends after start of growth pause
                    if (age < p.delayNGEnd) {
                        age__ = p.delayNGStart; // during growth pause
                    } else {
                        age__ = age - (p.delayNGEnd - p.delayNGStart); // simulation ends after end of growth pause
                    }
                } // delay to apply

                // we currently do not implement impeded growth for stem and leaves
                double targetlength = calcLength(age__) + this->epsilonDx;
                double e = targetlength - length; // store value of elongation to add
                double dl = e;                    // length increment = calculated length + increment from last time step too small to be added
                length = getLength(true);
                this->epsilonDx = 0.; // now it is "spent" on targetlength (no need for -this->epsilonDx in the following)

                // create geometry
                if (p.laterals) { // stem has laterals
                    /* basal zone */
                    if ((dl > 0) && (length < p.lb)) { // length is the current length of the root
                        if (length + dl <= p.lb) {
                            // std::cout << "Basal zone growth: dl = " << dl << std::endl << std::flush;
                            createSegments(dl, dt_);
                            length += dl;
                            dl = 0;
                        } else {
                            double ddx = p.lb - length;
                            // std::cout << "Basal zone growth: ddx = " << ddx << std::endl << std::flush;
                            createSegments(ddx, dt_);
                            dl -= ddx; // ddx already has been created
                            length = p.lb;
                        }
                    }
                    /* branching zone */
                    if (((created_linking_node) < (p.ln.size() + 1)) && (length >= p.lb)) {
                        for (size_t i = 0; (i < p.ln.size()); i++) {
                            createLateral(dt_, verbose);
                            if (p.ln.at(created_linking_node - 1) > 0) {
                                // std::cout << "Phytomer zone growth: " << this->dxMin() << std::endl << std::flush;
                                createSegments(this->dxMin(), dt_);
                                dl -= this->dxMin();
                                length += this->dxMin();
                            }
                        }
                        createLateral(dt_, verbose);
                    }
                    // we can have (p.ln.size()+1)>(created_linking_node) if one ln == 0cm
                    if ((length >= p.lb) && ((p.ln.size() + 1) < (created_linking_node))) {
                        std::stringstream errMsg;
                        errMsg << "Stem::simulate(): higher number of realized linking nodes (" << created_linking_node << ") than of max laterals (" << p.ln.size() + 1 << ")";
                        throw std::runtime_error(errMsg.str().c_str());
                    }
                    // internode elongation
                    double maxInternodeDistance = p.getK() - p.la - p.lb; // maximum length of branching zone
                    if ((dl > 0) && (length >= p.lb) && (maxInternodeDistance > 0)) {
                        int nn = localId_linking_nodes.back();                                      // node carrying the last lateral == end of branching zone
                        double currentInternodeDistance = getLength(nn) - p.lb;                     // actual length of branching zone
                        double ddx = std::min(maxInternodeDistance - currentInternodeDistance, dl); // length to add to branching zone
                        if (ddx > 0) {
                            // std::cout << "Phytomer elongation: " << ddx << std::endl << std::flush;
                            internodalGrowth(ddx, dt_);
                            dl -= ddx;
                            length += ddx;
                        }
                    }
                    /* apical zone */
                    // only grows once the basal and branching nodes are developped (TODO)
                    if ((dl > 0) && (length - (maxInternodeDistance + p.lb) > -1e-9)) {
                        // std::cout << "Apical zone growth: dl = " << dl << std::endl << std::flush;
                        createSegments(dl, dt_);
                        length += dl;
                    }
                } else { // no laterals
                    if (dl > 0) {
                        // std::cout << "No laterals: " << dl << std::endl << std::flush;
                        createSegments(dl, dt_);
                        length += dl;
                    }
                }
                if (dl < 0) {                  // to keep in memory that realised length is too long, as created nodes to carry children
                    this->epsilonDx = dl;      // targetlength + e - length;
                    length += this->epsilonDx; // go back to having length = theoratical length
                }
            } // if active
            // set limit below 1e-10, as the test files see if correct length
            // once rounded at the 10th decimal
            active = getLength(false) <= (p.getK() * (1 - 1e-11)); // become inactive, if final length is nearly reached
        }
    } // if alive
}

/**
 * @brief Returns the growth increment vector for the next segment.
 *
 * Queries the tropism function (f_tf) to obtain heading angles alpha/beta,
 * then rotates the heading frame accordingly.
 *
 * @param p    Coordinates of the previous node.
 * @param sdx  Length of the next segment [cm].
 * @param n    Node index used to determine the current heading (-1 = last node).
 * @return     Increment vector of length @p sdx in the growth direction, relative to the heading frame (Matrix3d::ons of the heading).
 */
Vector3d Stem::getIncrement(const Vector3d &p, Vector3d &h, double sdx) {
    Matrix3d ons = Matrix3d::ons(h);
    Vector2d ab = getOrganRandomParameter()->f_tf->getHeading(p, ons, dx(), shared_from_this());
    return Vector3d::rotAB(ab.x, ab.y).times(sdx); // relative to the heading frame (ons)
}

/**
 * @brief Creates polyline segments totalling length @p l.
 *
 * Tip growth (phytoIdx < 0) appends nodes at the end of the stem. Internodal elongation (phytoIdx >= 0)
 * inserts nodes directly above node phytoIdx, i.e. at the base of the phytomere.
 * If the segment that is to be extended is shorter than dx() (first call of Stem::simulate for tip growth,
 * always for elongation), it is lengthened first, before new nodes are created at intervals of dx().
 * Residual growth smaller than dxMin() is stored in epsilonDx.
 *
 * @param l         Total length of new segments to create [cm].
 * @param dt        Current time step [days] (used for creation-time calculation).
 * @param phytoIdx  Node index below which the phytomere grows (internodal elongation); -1 for tip growth.
 */
void Stem::createSegments(double l, double dt, int phytoIdx) {
    if (l == 0) {
        std::cout << "Stem::createSegments: zero length encountered \n";
        return;
    }
    if (l < 0) {
        std::cout << "Stem::createSegments: negative length encountered \n";
        return;
    }
    const bool elongation = (phytoIdx >= 0);
    const int nn = elongation ? phytoIdx + 1 : nodes.size(); // index of the first new node

    // lengthen an existing segment up to dx() instead of creating a new node
    double shiftl = 0; // length produced by shift
    if (elongation || firstCall) {
        if (!elongation) {
            firstCall = false;
        }
        int shiftIdx; // node (relative coordinate) whose segment is lengthened
        bool canShift;
        if (elongation) { // first segment inside the phytomere (the segment below belongs to the previous phytomere)
            shiftIdx = nn;
            canShift = (shiftIdx < static_cast<int>(nodes.size()));
        } else { // last segment, unless it ends in the base node of a lateral
            shiftIdx = nn - 1;
            canShift = (nn > 1) && (children.empty() || (shiftIdx != children.back()->getParentNI()));
        }
        moved = false;
        if (canShift) {
            Vector3d h = nodes.at(shiftIdx);
            double olddx = h.length();
            if (olddx < dx() * (1 - 1e-10)) {
                shiftl = std::min(dx() - olddx, l);
                h.normalize();
                nodes.at(shiftIdx) = h.times(olddx + shiftl);
                nodeCTs.at(shiftIdx) =
                    this->calcCreationTime(getLength(true) + shiftl, dt); // in case of impeded growth the node emergence time is not exact anymore, but might break down to temporal resolution
                moved = true;
                l -= shiftl;
                if (l <= 0) {
                    return;
                }
            }
        }
    }

    // create n+1 new nodes
    double sl = 0; // summed length of created segment
    int n = floor(l / dx());
    Vector3d h = heading(nn - 1); // heading and position of the node preceding the new one, updated incrementally
    Vector3d p = getNode(nn - 1);
    for (int i = 0; i < n + 1; i++) {

        double sdx;  // segment length (<=dx)
        if (i < n) { // normal case
            sdx = dx();
        } else { // last segment
            sdx = l - n * dx();
            if (sdx < dxMin() * (1 - 1e-10)) { // quit if l is too small
                if (elongation) {
                    this->epsilonDx += sdx;
                } else {
                    this->epsilonDx = sdx;
                }
                return;
            }
            this->epsilonDx = 0; // no residual
        }
        sl += sdx;
        // the new node is inserted at index nn+i, its predecessor is the node nn+i-1
        Vector3d newnode = getIncrement(p, h, sdx);                            // relative coordinate
        double et = this->calcCreationTime(getLength(true) + shiftl + sl, dt); // here length or getLength? it's the same because epsilonDx was set back to 0 at beginning of simulate?
        if (elongation) {
            // std::cout << "Inserted node at index " << phytoIdx << ", " << nn << ", " << et << ", " << sl << std::flush << std::endl;
            insertNode(newnode, et, nn);
            // insertion point p and heading h remain unchanged
        } else {
            addNode(newnode, et);
            Vector3d s = Matrix3d::ons(h).times(newnode); // absolute segment vector
            p = p.plus(s);                                // update position and heading
            h = s.normalized();
        }
    }
}

/**
 * Creates the laterals of a new linking node (@see Organ::createLateral) and stores the local id of the node
 */
void Stem::createLateral(double dt, bool verbose) {
    Organ::createLateral(dt, verbose);
    storeLinkingNodeLocalId(created_linking_node);
}

/**
 *  @see Organ::createLateral
 *  @param dt       time step recieved by parent organ [day]
 *  @return growth period to send to lateral after creation
 */
double Stem::getLatInitialGrowth(double dt) {
    double ageLN = this->calcAge(param()->lb); // MINIMUM age of root when lateral node is created
    ageLN = std::max(ageLN, age - dt);
    return age - ageLN;
}

/**
 * Overwritten to support different delay times for laterals along the organ,
 * even though all laterals are created at once when the stem emerges.
 *
 *  @see Organ::createLateral
 *  @param ot_lat       organType of lateral to create
 *  @param st_lat       subType of lateral to create
 *  @param dt       time step recieved by parent organ [day]
 *  @return emergence delay to send to lateral after creation
 */
double Stem::getLatGrowthDelay(int ot_lat, int st_lat, double dt) const { // override for stems
    bool verbose = false;
    auto rp = getOrganRandomParameter(); // rename
    double forDelay;                     // store necessary variables to define lateral growth delay
    int delayDefinition = getOrganism()->getDelayDefinition(ot_lat);

    assert(delayDefinition >= 0);

    if (verbose) {
        std::cout << "create lat, delay def " << delayDefinition << " " << getId() << " " << (nodes.size() - 1) << " " << age << " " << getNodeId(nodes.size() - 1) << " " << getNodeId(0) << std::endl;
    }

    // count the number of laterals of subtype st already created on this organ std::function<double(int, int, std::shared_ptr<Organ>)>
    auto correctST = [ot_lat, st_lat](std::shared_ptr<Organ> org) -> double {
        return double((org->getParameter("organType") == ot_lat) && (org->getParameter("subType") == st_lat));
    }; // return 1. if organ of correct type and subtype, 0. otherwise
    double multiplyDelay = double(std::count_if(children.begin(), children.end(),
                                                correctST)); // coun how may children of the correct type and subtype exists (they do not have to be emerged)

    switch (delayDefinition) {
    case Organism::dd_distance: {
        double meanLn = getParameter("lnMean");                             // mean inter-lateral distance
        double effectiveLa = std::max(getParameter("la") - meanLn / 2, 0.); // effective apical distance, observed apical distance is in [la-ln/2, la+ln/2]
        if (param()->lb + effectiveLa == param()->getK()) {
            effectiveLa /= 2;
        } // otherwise the growth delay == parent root life time.
        if (verbose) {
            std::cout << "case Organism::dd_distance " << organType() << " " << getParameter("subType") << " " << getLength(true) << " " << effectiveLa << " " << getParameter("la") << " " << meanLn
                      << std::endl;
        }
        double ageLN = this->calcAge(param()->lb); // age of root when lateral node is created
        ageLN = std::max(ageLN, age - dt);
        double ageLG = this->calcAge(param()->lb + effectiveLa); // age of the root, when the lateral starts growing (i.e when the apical zone is developed)
        forDelay = ageLG - ageLN;                                // time the lateral has to wait
        multiplyDelay = 1;                                       // in this case, even for stems, it does not matter how many laterals there were before.
        break;
    }
    case Organism::dd_time_lat: { // lateral has to wait for a fixed time (defined by parent organ)
        // time the lateral has to wait
        forDelay = std::max(rp->ldelay + plant.lock()->randn() * rp->ldelays, 0.);
        multiplyDelay *= rp->multDelay;
        if (verbose) {
            std::cout << "Organism::dd_time_lat " << rp->ldelay << " " << rp->ldelays << " " << forDelay << std::endl;
        }
        break;
    }
    case Organism::dd_time_self: {                                          // lateral has to wait for a fixed time (defined by lateral organ itself)
        auto latRp = plant.lock()->getOrganRandomParameter(ot_lat, st_lat); // random parameter of lateral to create
        forDelay = std::max(latRp->ldelay + plant.lock()->randn() * latRp->ldelays, 0.);
        multiplyDelay *= latRp->multDelay;
        if (verbose) {
            std::cout << "create lat, delay output " << forDelay << std::endl;
            std::cout << "						 " << ot_lat << ", " << st_lat << " " << latRp->ldelay << " " << latRp->ldelays << " " << forDelay << " " << nodes.size()
                      << std::endl;
        }
        break;
    }
    default: {
        std::cout << "delayDefinition " << delayDefinition << " " << Organism::dd_distance << " ";
        std::cout << Organism::dd_time_lat << " " << Organism::dd_time_self << std::endl << std::flush;
        std::cout << "				" << (delayDefinition == Organism::dd_distance) << " ";
        std::cout << (delayDefinition == Organism::dd_time_lat) << " " << (delayDefinition == Organism::dd_time_self) << std::endl << std::flush;
        throw std::runtime_error("Delay definition type (delayDefinition) not recognised");
    }
    }
    if (verbose) {
        std::cout << "create lat, delay defEND " << forDelay << " " << multiplyDelay << std::endl;
    }
    return forDelay * multiplyDelay;
}
/**
 * Simulates internodal growth of dl for this stem
 * divid total stem growth between the phytomeres
 * currently two option:
 * (p.nodalGrowth == 0) the phytomere grow sequentially
 * (p.nodalGrowth == 1) growth devided equally between the phytomeres or
 *
 * @param 	dl			total length of the segments that are created [cm]
 */
void Stem::internodalGrowth(double dl, double dt) {
    const StemSpecificParameter &p = *param(); // rename
    std::vector<double> toGrow(p.ln.size());
    double dl_;
    const int ln_0 = std::count(p.ln.cbegin(), p.ln.cend(), 0); // number of laterals which grow from the same branching point as the one before
    if (p.nodalGrowth == 0) {                                   // sequentiall growth
        // std::cout << "Sequential growth: dl = " << dl << std::endl << std::flush;
        toGrow[0] = dl;
        std::fill(toGrow.begin() + 1, toGrow.end(), 0);
    }
    if (p.nodalGrowth == 1) { // equal growth
        // std::cout << "Equal growth: dl = " << dl << std::endl << std::flush;
        std::fill(toGrow.begin(), toGrow.end(), dl / (p.ln.size() - ln_0));
    }
    int loopId = 0;
    size_t phytomerId = 0;
    while ((dl > 0) && (loopId < 2)) { // do the loop at most twice over the children
        // if the phytomere can do a growth superior to the mean phytomere growth, we add the value of "missing"
        //(i.e., length left to grow to get the predefined total growth of the branching zone)
        int nn1 = localId_linking_nodes.at(phytomerId);     // node at the beginning of phytomere
        int nn2 = localId_linking_nodes.at(phytomerId + 1); // node at end of phytomere (if nn1 != nn2)

        double length1 = getLength(nn1);
        double availableForGrowth = p.ln.at(phytomerId) - (getLength(nn2) - length1); // difference between maximum and current length of the phytomer
        if (availableForGrowth < -1e-3) {
            std::stringstream errMsg;
            errMsg << "Stem::internodalGrowth phytomere " << phytomerId << " is too long: " << availableForGrowth << " " << p.ln.at(phytomerId) << " " << getLength(nn2) << " " << length1 << std::endl;
            throw std::runtime_error(errMsg.str().c_str());
        }
        dl_ = std::max(0., std::min(std::min(toGrow[phytomerId], availableForGrowth), dl));
        if (dl_ > 0) {
            createSegments(dl_, dt, nn1);
            dl -= dl_;
        }
        if ((phytomerId + 1) < p.ln.size()) {
            toGrow.at(phytomerId + 1) += toGrow.at(phytomerId) - dl_;
            phytomerId++;
        } else {
            toGrow.at(0) += toGrow.at(phytomerId) - dl_;
            loopId++;
            phytomerId = 0;
        } // loop twice
    }
    if (std::abs(dl) > 1e-6) { // this sould not happen as computed dl to be <= sum(availableForGrowth)
        std::stringstream errMsg;
        errMsg << "Stem::internodalGrowth length left to grow: " << dl;
        throw std::runtime_error(errMsg.str().c_str());
    }
}

/**
 * Returns a parameter per organ
 *
 * @param name 		parameter name (returns nan if not available)
 *
 */
double Stem::getParameter(std::string name, std::map<std::string, double> addParams) const {
    if (name == "lb") {
        return param()->lb;
    } // basal zone [cm]
    if (name == "delayNGStart") {
        return param()->delayNGStart;
    } // delay for nodal growth [day]
    if (name == "delayNGEnd") {
        return param()->delayNGEnd;
    } // delay for nodal growth [day]
    if (name == "la") {
        return param()->la;
    } // apical zone [cm]
    if (name == "nob") {
        return param()->nob();
    } // number of branching points
    if (name == "r") {
        return param()->r;
    } // initial growth rate [cm day-1]
    if (name == "radius") {
        return param()->a;
    } // root radius [cm]
    if (name == "a") {
        return param()->a;
    } // root radius [cm]
    if (name == "theta") {
        return param()->theta;
    } // angle between root and parent root [rad]
    if (name == "rlt") {
        return param()->rlt;
    } // root life time [day]
    if (name == "k") {
        return param()->getK();
    }; // maximal root length [cm]
    if (name == "lnMean") { // mean lateral distance [cm]
        auto &v = param()->ln;
        if (v.size() > 0) {
            return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
        } else {
            return 0;
        }
    }
    if (name == "lnDev") { // standard deviation of lateral distance [cm]
        auto &v = param()->ln;
        double mean = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
        double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
        return std::sqrt(sq_sum / v.size() - mean * mean);
    }
    if (name == "volume") {
        return param()->a * param()->a * M_PI * getLength(true);
    } // // root volume [cm^3]
    if (name == "surface") {
        return 2 * param()->a * M_PI * getLength(true);
    }
    if (name == "type") {
        return this->param_->subType;
    } // delete to avoid confusion?
    if (name == "subType") {
        return this->param_->subType;
    } // organ sub-type [-]
    if (name == "parentNI") {
        return parentNI;
    } // local parent node index where the lateral emerges
    return Organ::getParameter(name, addParams);
}

/**
 * Analytical length of the stem at a given age
 *
 * @param age          age of the stem [day]
 */
double Stem::calcLength(double age) {
    assert(age >= 0 && "Stem::calcLength() negative root age");
    return getStemRandomParameter()->f_gf->getLength(age, param()->r, param()->getK(), shared_from_this());
}

/**
 * Analytical age of the stem at a given length
 * no scaling of organ growth , so can return age directly
 * otherwise cannot compute exact age between delayNGStart and delayNGEnd
 * @param length   length of the stem [cm]
 */
double Stem::calcAge(double length) const {
    assert(length >= 0 && "Stem::calcAge() negative root age");
    double age__ = getStemRandomParameter()->f_gf->getAge(length, param()->r, param()->getK(), shared_from_this());
    if (age__ > param()->delayNGStart) {
        age__ += (param()->delayNGEnd - param()->delayNGStart);
    }
    return age__;
}

/**
 * stores the local id of the linking node. used by @see Stem::internodalGrowth()
 */
void Stem::storeLinkingNodeLocalId(int numCreatedLN) {
    localId_linking_nodes.push_back(nodes.size() - 1);
    if (numCreatedLN != localId_linking_nodes.size()) {
        throw std::runtime_error("wrong number of linking nodes in stem: " + std::to_string(numCreatedLN) + " against " + std::to_string(localId_linking_nodes.size()));
    }
}

/**
 * Inserts a node between existing nodes because of internodal growth, and updates the linking nodes.
 *
 * @see Organ::insertNode
 *
 * @param n        new node
 * @param id       global node index
 * @param t        exact creation time of the node
 * @param index	   position were new node is to be added
 */
void Stem::insertNode(Vector3d n, double t, size_t index) {
    Organ::insertNode(n, t, index);
    for (int i = 0; i < localId_linking_nodes.size(); i++) { // update the local ids of the linking nodes
        if ((localId_linking_nodes.at(i) >= index) && (localId_linking_nodes.at(i) > 0)) {
            localId_linking_nodes.at(i) += 1;
        }
    }
}

/**
 * @return The StemTypeParameter from the plant
 */
std::shared_ptr<StemRandomParameter> Stem::getStemRandomParameter() const {
    return std::static_pointer_cast<StemRandomParameter>(plant.lock()->getOrganRandomParameter(Organism::ot_stem, param_->subType));
}

/**
 * @return Parameters of the specific root
 */
std::shared_ptr<const StemSpecificParameter> Stem::param() const { return std::static_pointer_cast<const StemSpecificParameter>(param_); }

/*
 * Quick info about the object for debugging
 * additionally, use param()->toString() and getOrganRandomParameter()->toString() to obtain all information.
 */
std::string Stem::toString() const {
    std::stringstream newstring;
    newstring << "; initial heading: " << getInitialHeading().toString() << ", parent node index" << parentNI << ".";
    return Organ::toString() + newstring.str();
}

OrganState Stem::getState() const {
    OrganState s = Organ::getState();
    stateSet(s, "stemphytomerId", stemphytomerId);
    stateSet(s, "localId_linking_nodes", localId_linking_nodes);
    if (auto p = param()) {
        stateSet(s, "param.lb", p->lb);
        stateSet(s, "param.la", p->la);
        stateSet(s, "param.ln", p->ln);
        stateSet(s, "param.r", p->r);
        stateSet(s, "param.theta", p->theta);
        stateSet(s, "param.rlt", p->rlt);
        stateSet(s, "param.laterals", p->laterals);
        stateSet(s, "param.nodalGrowth", p->nodalGrowth);
        stateSet(s, "param.delayNGStart", p->delayNGStart);
        stateSet(s, "param.delayNGEnd", p->delayNGEnd);
    }
    return s;
}

void Stem::setState(const OrganState &state, const std::map<int, std::shared_ptr<Organ>> &organs) {
    Organ::setState(state, organs);
    if (param_) {
        auto p = std::make_shared<StemSpecificParameter>(*param());
        stateGet(state, "subType", p->subType);
        stateGet(state, "param.a", p->a);
        stateGet(state, "param.lb", p->lb);
        stateGet(state, "param.la", p->la);
        stateGet(state, "param.ln", p->ln);
        stateGet(state, "param.r", p->r);
        stateGet(state, "param.theta", p->theta);
        stateGet(state, "param.rlt", p->rlt);
        stateGet(state, "param.laterals", p->laterals);
        stateGet(state, "param.nodalGrowth", p->nodalGrowth);
        stateGet(state, "param.delayNGStart", p->delayNGStart);
        stateGet(state, "param.delayNGEnd", p->delayNGEnd);
        param_ = p;
    }
    stateGet(state, "stemphytomerId", stemphytomerId);
    stateGet(state, "localId_linking_nodes", localId_linking_nodes);
}

} // namespace CPlantBox
