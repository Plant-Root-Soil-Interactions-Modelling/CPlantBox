#include "Stem.h"
#include "Leaf.h"
#include "Root.h"
/**
* Constructor
* This is a Copy Paste of the Root.cpp but it works independently, it has its own parameter file (in .stparam file) tropism, growth function, txt and vtp writing system.
* All of those can be modified to fit the real growth of the Plant.
*
* Typically called by the Plant::Plant(), or Stem::createNewStem().
* For stem the initial node and node emergence time (netime) must be set from outside
*
* @param rs 			points to Stem
* @param type 		    type of stem that is created
* @param pheading		heading of parent stem at emergence
* @param delay 		to give apical zone of parent time to develop
* @param parent		parent stem
* @param pni			parent node index
* @param pbl			parent base length
*/
Stem::Stem(Plant* plant, Organ* parent, int type, double delay, Vector3d isheading ,int pni, double pbl) :Organ(plant,parent,type,delay), pni(pni), pbl(pbl)
{

  initialstemHeading=isheading;
  std::cout << "stem pni = "<< pni<< std::endl;
//  std::cout << "Stem constructor \n";
  StemTypeParameter* sttp = (StemTypeParameter*) plant->getOrganTypeParameter(Plant::ot_stem, type);


  stem_param = sttp->realize(); // throw the dice
  StemParameter* stem_p = (StemParameter*) stem_param;
//  std::cout <<", "<<(StemParameter*) stem_param<< "\n";
  double beta = 2*M_PI*plant->rand(); // initial rotation
  Matrix3d ons = Matrix3d::ons(initialstemHeading);
  ons.times(Matrix3d::rotX(beta));
  double theta = stem_p->theta;
  if (parent->organType()!=Plant::ot_seed) { // scale if not a base stem
    double scale = sttp->sa->getValue(parent->getNode(pni),this);
    theta*=scale;
  }
  ons.times(Matrix3d::rotZ(theta));
  // initial node
  if (parent->organType()!=Plant::ot_seed) { // the first node of the base stems must be created in Seed::initialize()
    // otherwise, don't use addNode for the first node of the stem,
    // since this node exists already and does not need a new identifier
    r_nodes.push_back(parent->getNode(pni));
    nodeIDs.push_back(parent->getNodeID(pni));
    nctimes.push_back(parent->getNodeCT(pni)+delay);
  }
}

int Stem::organType() const
{
  return Plant::ot_stem;
}

/**
* Simulates growth of this stem for a time span dt
*
* @param dt       time step [day]
* @param silence  indicates if status messages are written to the console (cout) (default = false)
*/
void Stem::simulate(double dt, bool silence)
{

  old_non = 0; // is set in Stem:createSegments, the zero indicates the first call to createSegments
  Vector3d ilheading(0,0,1);


  Leaf* LeafGrow = new Leaf(plant, this , 2, 0., ilheading ,r_nodes.size()-1, 20.);

//    Vector3d ilheading(0,0,1);//Initial Stem heading
//	Leaf* mainleaf = new Leaf(plant, this, 1, 0., ilheading ,0., 0.); // tap root has subtype 1
//	mainleaf->addNode(sparam->seedPos,0);
//	children.push_back(mainleaf);




  const StemParameter* sp = sParam(); // rename
  const StemTypeParameter* sttp = stParam();

  // increase age
  if (age+dt>sp->rlt) { // stem life time
    dt=sp->rlt-age; // remaining life span
    alive = false; // this stem is dead
  }
  age+=dt;

  if (alive) { // dead stem wont grow
//createShootborneroot(silence);
    // probabilistic branching model (todo test)
    if ((age>0) && (age-dt<=0)) { // the stem emerges in this time step
      double stem_P = sttp->sbp->getValue(r_nodes.back(),this);
      if (stem_P<1.) { // P==1 means the lateral emerges with probability 1 (default case)
        double stem_p = 1.-std::pow((1.-stem_P), dt); //probability of emergence in this time step
//        std::cout <<stem_P<<", "<<stem_p<< "\n";
        if (plant->rand()>stem_p) { // not rand()<p
          age -= dt; // the stem does not emerge in this time step
        }
      }
    }

    if (age>0) {

      // children first (lateral stems grow even if base stem is inactive)
      for (auto c:children) {
        c->simulate(dt,silence);
      }

      if (active) {

        // length increment
        double length_ = StemgetLength(std::max(age-dt,0.)); // length of the stem for unimpeded growth (i.e. length_==length for unimpeded growth)
        double targetlength = StemgetLength(age);
        double e = targetlength-length_; //elongation in time step dt
        double scale = sttp->se->getValue(r_nodes.back(),this); // hope some of this is optimized out if not set
        double dl = std::max(scale*e, double(0)); // length increment, dt is not used anymore

        // create geometry
        if (sp->ln.size()>0) { // stem has laterals
          // basal zone
          if ((dl>0)&&(length<sp->lb)) { // length is the current length of the stem
            if (length+dl<=sp->lb) {
              createSegments(dl,silence);
              length+=dl;
              dl=0;
            } else {
              double ddx = sp->lb-length;
              createSegments(ddx,silence);
              dl-=ddx; // ddx already has been created
              length=sp->lb;
            }
          }
          // branching zone Condition will be changed later
          if ((dl>0)&&(length>=sp->lb)) {
            double s = sp->lb; // summed length
            for (size_t i=0; ((i<sp->ln.size()) && (dl>0)); i++) {
              s+=sp->ln.at(i);
              if (length<s) {
                if (i==children.size()) { // new lateral
LeafGrow->addNode(Stem::getNode(r_nodes.size()-1), 10.);
               children.push_back(LeafGrow);
//                 createLateral(silence);
                }
                if (length+dl<=s) { // finish within inter-lateral distance i
                  createSegments(dl,silence);
                  length+=dl;
                  dl=0;
                } else { // grow over inter-lateral distance i
                  double ddx = s-length;
                  createSegments(ddx,silence);
                  dl-=ddx;
                  length=s;
                }
              }
            }
            if (dl>0) {
              if (sp->ln.size()==children.size()) { // new lateral (the last one)
LeafGrow->addNode(Stem::getNode(children.size()), 10.);
                 children.push_back(LeafGrow);

//                LeafGrow->createLateral(silence);
//children.push_back(LeafGrow);
              }
            }
          }
          // apical zone
          if (dl>0) {
            createSegments(dl,silence);
            length+=dl;
          }
        } else { // no laterals
          if (dl>0) {
            createSegments(dl,silence);
            length+=dl;
          }
        } // if laterals
      } // if active
      active = StemgetLength(std::max(age,0.))<(sp->getK()-dx()/10); // become inactive, if final length is nearly reached
    }
  } // if alive



}

/**
*
*/
double Stem::getScalar(int stype) const {
  switch(stype) {
  // st_rlt, st_meanln, st_stdln , st_nob, st_surface, , // stem level
  case Plant::st_lb:
    return sParam()->lb;
  case Plant::st_la:
    return sParam()->la;
  case Plant::st_r:
    return sParam()->r;
  case Plant::st_radius:
    return sParam()->a;
  case Plant::st_theta:
    return sParam()->theta;
  case Plant::st_rlt:
    return sParam()->rlt;
  case Plant::st_meanln:
    return std::accumulate(sParam()->ln.begin(), sParam()->ln.end(), 0.0) / sParam()->ln.size();
  case Plant::st_stdln: {
    const std::vector<double>& v_ = sParam()->ln;
    double mean = std::accumulate(v_.begin(), v_.end(), 0.0) / v_.size();
    double sq_sum = std::inner_product(v_.begin(), v_.end(), v_.begin(), 0.0);
    return std::sqrt(sq_sum / v_.size() - mean * mean);
  }
  case Plant::st_surface:
    return sParam()->a*sParam()->a*M_PI*length;
  case Plant::st_nob:
    return sParam()->ln.size();
  default:
    return  Organ::getScalar(stype);
  }
}

/**
* Analytical creation (=emergence) time of a node at a length along the stem
*
* @param length   length of the stem [cm]
*/
double Stem::getCreationTime(double length)
{
  assert(length>=0);
  double stemage = StemgetAge(length);
  assert(stemage>=0);
  if (parent->organType()!=Plant::ot_seed) {
    if (parent->organType()==Plant::ot_shoot) {
      double pl = pbl+((Stem*)parent)->stParam()->la; // parent length, when this stem was created
      double pAge=((Stem*)parent)->getCreationTime(pl);
      return stemage+pAge;
    } else { // organ type is seed
      return stemage;
    }
  } else {
    return stemage;
  }
}

/**
* Analytical length of the stem at a given age
*
* @param age          age of the stem [day]
*/
double Stem::StemgetLength(double age)
{
  assert(age>=0);
  return stParam()->growth->StemgetLength(age,stParam()->r,stParam()->getK(),this);
}

/**
* Analytical age of the stem at a given length
*
* @param length   length of the stem [cm]
*/
double Stem::StemgetAge(double length)
{
  assert(length>=0);
  return stParam()->growth->StemgetAge(length,stParam()->r,stParam()->getK(),this);
}

/**
*
*/
StemTypeParameter* Stem::stParam() const {
  return (StemTypeParameter*)getOrganTypeParameter();
}

/**
*
*/
double Stem::dx() const
{
  return ((StemTypeParameter*)getOrganTypeParameter())->dx;
}

/**
* Creates a new lateral by calling Stem::createNewstem().
*
* Overwrite this method to implement more specialized stem classes.
*/
void Stem::createLateral(bool silence)
{
  const StemParameter* sp = sParam(); // rename
  int lt = stParam()->getLateralType(r_nodes.back());
  	std::cout << "Stem createLateral()\n";
  	std::cout << "Stem lateral type " << lt << "\n";

  if (lt>0) {
    double ageLN = this->StemgetAge(length); // age of stem when lateral node is created
    double ageLG = this->StemgetAge(length+sp->la); // age of the stem, when the lateral starts growing (i.e when the apical zone is developed)
    double delay = ageLG-ageLN; // time the lateral has to wait
    Vector3d h = heading(); // current heading
    Stem* lateral = new Stem(plant, this, lt, delay, h,  r_nodes.size()-1, length);
    children.push_back(lateral);
    lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
    }


}



/**
*  Creates nodes and node emergence times for length l,
*  and updates the stem heading
*
*  Cecks that each new segments length is <= dx but >= ddx
*
*  @param l       length the stem growth [cm]
*/
void Stem::createSegments(double l, bool silence)
{
   std::cout << "create Stem Segments("<< l << ")\n";
  assert(l>0);
  double sl=0; // summed length of created segment

  // shift first node to axial resolution
  int nn = r_nodes.size();
  if (old_non==0) { // first call of createSegments (in stem::simulate)
    if (nn>1) {
      auto n2 = r_nodes.at(nn-2);
      auto n1 = r_nodes.at(nn-1);
      double olddx = n1.minus(n2).length();
      if (olddx<dx()*0.99) { // shift node instead of creating a new node

        Vector3d h; // current heading
        if (nn>2) {
          h = n2.minus(r_nodes.at(nn-3));
          h.normalize();
        } else {
          h = initialstemHeading;
        }
        double sdx = std::min(dx()-olddx,l);

        Matrix3d ons = Matrix3d::ons(h);
        Vector2d ab = stParam()->tropism->getHeading(r_nodes.at(nn-2),ons,olddx+sdx,this);
        ons.times(Matrix3d::rotX(ab.y));
        ons.times(Matrix3d::rotZ(ab.x));
        Vector3d newdx = Vector3d(ons.column(0).times(olddx+sdx));

        Vector3d newnode = Vector3d(r_nodes.at(nn-2).plus(newdx));
        sl = sdx;
        double ct = this->getCreationTime(length+sl);
        r_nodes[nn-1] = newnode;
        nctimes[nn-1] = std::max(ct,plant->getSimTime()); // in case of impeded growth the node emergence time is not exact anymore, but might break down to temporal resolution
        old_non = nn-1;
        l -= sdx;
        if (l<=0) { // ==0 should be enough
          return;
        }
      }
    }
    old_non = nn; // CHECK
  }

  if (l<smallDx) {
    if (!silence) {
      std::cout << "skipped small segment l<Dx (<"<< smallDx << ") \n";
    }
    return;
  }

  int n = floor(l/dx());
  // create n+1 new nodes
  for (int i=0; i<n+1; i++) {

    double sdx; // segment length (<=dx)
    if (i<n) {  // normal case
      sdx = dx();
    } else { // last segment
      sdx = l-n*dx();
      if (sdx<smallDx) {
        if (!silence) {
          std::cout << "skipped small segment i<n (<"<< smallDx << ") \n";
        }
        return;
      }
    }
    sl+=sdx;

    Vector3d h= heading();
    Matrix3d ons = Matrix3d::ons(h);
    Vector2d ab = stParam()->tropism->getHeading(r_nodes.back(),ons,sdx,this);
    ons.times(Matrix3d::rotX(ab.y));
    ons.times(Matrix3d::rotZ(ab.x));
    Vector3d newdx = Vector3d(ons.column(0).times(sdx));
    Vector3d newnode = Vector3d(r_nodes.back().plus(newdx));
    double ct = this->getCreationTime(length+sl);
    ct = std::max(ct,plant->getSimTime()); // in case of impeded growth the node emergence time is not exact anymore, but might break down to temporal resolution
     std::cout<<"add node "<<newnode.toString()<<"\n";
    addNode(newnode,ct);

  } // for

}
//***********************stem heading*****************************
Vector3d Stem::heading() const {
  Vector3d h;
  if ((stem_param->subType == 1) || (this->r_nodes.size()<=1) ) {// Make heading upward if it is main stem and
    h = initialstemHeading;// getHeading(b-a)
  } else {
     h = r_nodes.back().minus(r_nodes.at(r_nodes.size()-2));
  }
  return h;
}

/**
* Adds the next node to the stem.
*
* Add nodes only with this function! For simplicity nodes can not be deleted, and stems can only become deactivated by dying
*
* @param n        the new node
* @param t        exact creation time of the node
*/
void Stem::addNode(Vector3d n, double t)
{
  assert(t>=0.);
  r_nodes.push_back(n); // node
  nodeIDs.push_back(plant->getNodeIndex()); // new unique id
  nctimes.push_back(t); // exact creation time
}

/**
* writes RSML stem tag
*
* @param cout      typically a file out stream
* @param indent    we care for looks
*/
void Stem::writeRSML(std::ostream & cout, std::string indent) const
{
  if (this->r_nodes.size()>1) {
    cout << indent << "<stem id=\"" <<  id << "\">\n";  // open stem

      /* geometry tag */
      cout << indent << "\t<geometry>\n"; // open geometry
      cout << indent << "\t\t<polyline>\n"; // open polyline
      // polyline nodes
      cout << indent << "\t\t\t" << "<point ";
      Vector3d v = r_nodes.at(0);
      cout << "x=\"" << v.x << "\" y=\"" << v.y << "\" z=\"" << v.z << "\"/>\n";
      int n = this->plant->rsmlReduction;
      for (size_t i = 1; i<r_nodes.size()-1; i+=n) {
        cout << indent << "\t\t\t" << "<point ";
        Vector3d v = r_nodes.at(i);
        cout << "x=\"" << v.x << "\" y=\"" << v.y << "\" z=\"" << v.z << "\"/>\n";
      }
      cout << indent << "\t\t\t" << "<point ";
      v = r_nodes.at(r_nodes.size()-1);
      cout << "x=\"" << v.x << "\" y=\"" << v.y << "\" z=\"" << v.z << "\"/>\n";
      cout << indent << "\t\t</polyline>\n"; // close polyline
      cout << indent << "\t</geometry>\n"; // close geometry

      /* properties */
      cout << indent <<"\t<properties>\n"; // open properties
      // TODO
      cout << indent << "\t</properties>\n"; // close properties

      cout << indent << "\t<functions>\n"; // open functions
      cout << indent << "\t\t<function name='emergence_time' domain='polyline'>\n"; // open functions
      cout << indent << "\t\t\t" << "<sample>" << nctimes.at(0) << "</sample>\n";
      for (size_t i = 1; i<nctimes.size()-1; i+=n) {
        cout << indent << "\t\t\t" << "<sample>" << nctimes.at(i) << "</sample>\n";

      }
      cout << indent << "\t\t\t" << "<sample>" << nctimes.at(nctimes.size()-1) << "</sample>\n";

      cout << indent << "\t\t</function>\n"; // close functions
      cout << indent << "\t</functions>\n"; // close functions

      /* laterals stems */
      for (size_t i = 0; i<children.size(); i++) {
        children[i]->writeRSML(cout,indent+"\t");
      }

      cout << indent << "</root>\n"; // close stem
  }
}

/**
 * Quick info about the object for debugging
 */
std::string Stem::toString() const
{
  std::stringstream str;
  str << "Root #"<< id <<": type "<<stem_param->subType << ", length: "<< length << ", age: " <<age<<" with "<< children.size() << " laterals\n";
  return str.str();
}


