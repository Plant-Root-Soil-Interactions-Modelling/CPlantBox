#include "Organ.h"

/**
 * Constructor
 */
Organ::Organ(Plant* plant, Organ* parent, int type, double delay) :plant(plant), parent(parent), age(delay)
{
  id = plant->getOrganIndex();

}

/**
 * Destructor, tell the kids
 */
Organ::~Organ()
{
  for(auto c : children) {
      delete c;
  }
  delete param;
}

/**
 * returns the type of the organ
 */
int Organ::organType() const {
	return Plant::ot_organ;
}


/**
 * todo
 */
OrganTypeParameter* Organ::getOrganTypeParameter() const
{
  return plant->getOrganTypeParameter(param->organType, param->type);
};

/**
 * todo test...
 */
Vector3d Organ::getAbsoluteOrigin() const {
  const Organ* p = this;
  Vector3d o0 = Vector3d();
  while (p->organType() != Plant::ot_seed) {
      Vector3d o;
      if (p->parent!=0) {
          o = parent->getRelativeInitialHeading().times(p->getRelativeOrigin());
      } else {
          o = p->getRelativeOrigin();
      }
      o0.plus(o);
      p = p->parent;
  }
  return o0;
}

/**
 * todo test, multiplication from left?
 */
Matrix3d Organ::getAbsoluteInitialHeading() const {
  const Organ* p = this;
  Matrix3d ah = Matrix3d();
  while (p->organType() != Plant::ot_seed) {
      Matrix3d iH = p->getRelativeInitialHeading();
      iH.times(ah);
      ah = iH;
      p = p->parent;
  }
  return ah;
}

/**
 *
 */
double Organ::getScalar(int stype) {
  switch(stype) {
    case Plant::st_id:
      return id;
    case Plant::st_otype:
      return param->organType;
    case Plant::st_type:
      return param->type;
    case Plant::st_alive:
      return alive;
    case Plant::st_active:
      return active;
    case Plant::st_age:
      return age;
    case Plant::st_length:
      return length;
    default:
      throw std::invalid_argument( "Organ::getScalar: unknown scalar type" );
  }
}

/**
 * Returns the organs as sequential list,
 * copies only organs with more than 1 node.
 *
 * \return sequential list of organs
 */
std::vector<Organ*> Organ::getOrgans(int otype)
{
  std::vector<Organ*> v = std::vector<Organ*>();
  this->getOrgans(otype, v);
  return v;
}

/**
 * Returns the organs as sequential list,
 * copies only organs with more than 1 node.
 *
 * @param v     adds the organ sub tree to this vector
 */
void Organ::getOrgans(int otype, std::vector<Organ*>& v)
{
  if (this->r_nodes.size()>1) {
      int ot = this->param->organType;
      switch(otype) {
        case Plant::ot_organ:
          v.push_back(this);
          break;
        case Plant::ot_seed:
          if (ot==Plant::ot_seed) {
              v.push_back(this);
          }
          break;
        case Plant::ot_root:
          if (ot==Plant::ot_root) {
              v.push_back(this);
          }
          break;
        case Plant::ot_stem:
          if (ot==Plant::ot_stem) {
              v.push_back(this);
          }
          break;
        case Plant::ot_leafe:
          if (ot==Plant::ot_leafe) {
              v.push_back(this);
          }
          break;
        case Plant::ot_shoot:
          if ((ot==Plant::ot_leafe)||(ot==Plant::ot_stem)) {
              v.push_back(this);
          }
          break;
        default:
          throw std::invalid_argument( "Organ::getOrgans: unknown organ type" );
      }
  }
  for (auto const& c:this->children) {
      c->getOrgans(otype,v);
  }
}

/**
 * writes RSML root tag (todo update for general organs)
 *
 * @param cout      typically a file out stream
 * @param indent    we care for looks
 */
void Organ::writeRSML(std::ostream & cout, std::string indent) const
{
  if (this->r_nodes.size()>1) {
      cout << indent << "<root id=\"" <<  id << "\">\n";  // open root

      /* geometry tag */
      cout << indent << "\t<geometry>\n"; // open geometry
      cout << indent << "\t\t<polyline>\n"; // open polyline
      // polyline nodes
      cout << indent << "\t\t\t" << "<point ";
      Vector3d v = this->getNode(0);
      cout << "x=\"" << v.x << "\" y=\"" << v.y << "\" z=\"" << v.z << "\"/>\n";
      int n = this->plant->rsmlReduction;
      for (size_t i = 1; i<r_nodes.size()-1; i+=n) {
          cout << indent << "\t\t\t" << "<point ";
          Vector3d v = this->getNode(i);
          cout << "x=\"" << v.x << "\" y=\"" << v.y << "\" z=\"" << v.z << "\"/>\n";
      }
      cout << indent << "\t\t\t" << "<point ";
      v = r_nodes.back();
      cout << "x=\"" << v.x << "\" y=\"" << v.y << "\" z=\"" << v.z << "\"/>\n";
      cout << indent << "\t\t</polyline>\n"; // close polyline
      cout << indent << "\t</geometry>\n"; // close geometry

      /* properties */
      cout << indent <<"\t<properties>\n"; // open properties
      // TODO
      cout << indent << "\t</properties>\n"; // close properties

      cout << indent << "\t<functions>\n"; // open functions
      cout << indent << "\t\t<function name='emergence_time' domain='polyline'>\n"; // open functions
      cout << indent << "\t\t\t" << "<sample>" << nctimes[0] << "</sample>\n";
      for (size_t i = 1; i<nctimes.size()-1; i+=n) {
          cout << indent << "\t\t\t" << "<sample>" << nctimes.at(i) << "</sample>\n";

      }
      cout << indent << "\t\t\t" << "<sample>" << nctimes.back() << "</sample>\n";

      cout << indent << "\t\t</function>\n"; // close functions
      cout << indent << "\t</functions>\n"; // close functions

      /* laterals roots */
      for (size_t i = 0; i<children.size(); i++) {
          children[i]->writeRSML(cout,indent+"\t");
      }

      cout << indent << "</root>\n"; // close root
  }
}

/**
 * Quick info about the object for debugging
 */
std::string Organ::toString() const
{
  std::stringstream str;
  str << "Organ #"<< id <<": type "<< param->type << ", length: "<< length << ", age: " << age
      <<" with "<< children.size() << " successors\n";
  return str.str();
}


