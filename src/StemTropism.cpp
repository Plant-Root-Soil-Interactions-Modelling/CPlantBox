#include "StemTropism.h"

namespace CPlantBox {




/**
* Applies angles a and b and goes dx [cm] into the new direction and returns the new position
*
* @param pos          root tip position
* @param old          rotation matrix, heading is old(:,1)
* @param a            angle alpha, describes the angular rotation of the axis old
* @param b            angle beta, describes the radial rotation of the axis old
* @param dx           distance to look ahead
*
* \return             position
*/
Vector3d StemTropismFunction::getPosition(const Vector3d& pos, Matrix3d old, double a, double b, double dx)
{
  old.times(Matrix3d::rotX(b)); // TODO not nice, there should be a better way, but no easy
  old.times(Matrix3d::rotZ(a));
  return pos.plus(old.column(0).times(dx));
}

/**
* Dices N times picking angles alpha and beta, takes the optimal direction according to the objective function
*
* @param pos          root tip position
* @param old          rotation matrix, heading is old(:,1)
* @param dx           distance to look ahead (e.g in case of hydrotropism)
* @param root         points to the root that called getHeading(...), just in case something else is needed (i.e. iheading for exotropism)
*
* \return             angle alpha and beta
*/


///---------------------------------------------------------------Tropism Modified and Used--------------------------------------------------------------///
Vector2d StemTropismFunction::getHeading(const Vector3d& pos, Matrix3d old, double dx,const Organ* stem)
{
  //std::cout<<"StemTropismFunction::getHeading()\n";
  double a = 0;//sigma*randn()*sqrt(dx);
  double b = 0;//rand()*2*M_PI;
  double v;

  double n_=n*sqrt(dx);
  if (n_>0) {
    double dn = n_-floor(n_);
    if (rand()<dn) {
      n_ = ceil(n_);
    } else {
      n_ = floor(n_);
    }
    double bestA = a;
    double bestB = b;
    double bestV = this->stemtropismObjective(pos,old,a,b,dx,stem);
    for (int i=0; i<n_; i++) {
      b = 2*rand()*M_PI;///rand()*2*M_PI;
      a = sigma*randn()*sqrt(dx);///sigma*randn()*sqrt(dx);
      v = this->stemtropismObjective(pos,old,a,b,dx,stem);
      if (v<bestV) {
        bestV=v;
        bestA=a;
        bestB=b;
      }
    }
    a = bestA;
    b = bestB;
  }

  return Vector2d(a,b);
}



/**
* Inside the geometric domaint baseTropism::getHeading() is returned.
* In case geometric boundaries are hit, the rotations are modified, so that growth stays inside the domain.
*
* @param pos        root tip postion
* @param old        rotation matrix, heading is old(:,1)
* @param dx         distance to look ahead (e.g in case of hydrotropism)
* @param root       points to the root that called getHeading(), just in case something else is needed (i.e. iheading for exotropism)
*
* \return           the rotations alpha and beta
*/
Vector2d ConfinedStemTropism::getHeading(const Vector3d& pos, Matrix3d old, double dx, const Organ* stem)
{
  std::cout<<"ConfinedStemTropism::getHeading()\n";

  Vector2d h = stemtropism->getHeading(pos, old, dx, stem);

  double a = h.x;
  double b = h.y;
  double d = geometry->getDist(this->getPosition(pos,old,a,b,dx));
  double dmin = d;

  double bestA = a;
  double bestB = b;
  int i=0; // counts change in alpha
  int j=0;    // counts change in beta

  while (d>0) { // not valid

    i++;
    j=0;
    while ((d>0) && j<betaN) { // change beta

      b = 2*M_PI*rand(); // dice
      d = geometry->getDist(this->getPosition(pos,old,a,b,dx));
      if (d<dmin) {
        dmin = d;
        bestA = a;
        bestB = b;
      }
      j++;
    }

    if (d>0) {
      a = a + M_PI/2./double(alphaN);
    }

    if (i>alphaN) {
      std::cout << "Could not respect geometry boundaries \n";
      a = bestA;
      b = bestB;
      break;
    }

  }
  return Vector2d(a,b);
}



/**
* getHeading() minimizes this function, @see TropismFunction::tropismObjective
*/
double StemExotropism::stemtropismObjective(const Vector3d& pos, Matrix3d old, double a, double b, double dx, const Organ* stem)
{
  old.times(Matrix3d::rotX(b));
  old.times(Matrix3d::rotZ(a));
  Vector3d isheading = ((Stem*)stem)->A.column(0);
  double s = isheading.times(old.column(0));
  s*=(1./isheading.length()); // iheading should be normed anyway?
  s*=(1./old.column(0).length());
  return acos(s)/M_PI; // 0..1
}



/**
* getHeading() minimizes this function, @see TropismFunction::tropismObjective
*/
double StemPhototropism::stemtropismObjective(const Vector3d& pos, Matrix3d old, double a, double b, double dx, const Organ* stem)
{
  assert(soil!=nullptr);
  Vector3d newpos = this->getPosition(pos,old,a,b,dx);
  double v = soil->getValue(newpos,stem);
//   std::cout << "\n" << newpos.getString() << ", = "<< v;
  return -v; ///< (-1) because we want to maximize the soil property
}


/**
* Constructs a combined tropism with two tropsims,
* the new objective funciton is (t1->tropismObjective()*w1)+(t2->tropismObjective()*w2)
*
* @param n             number of tries
* @param sigma         standard deviation of angular change [1/cm]
* @param t1            first tropism
* @param w1            first weight
* @param t2            first tropism
* @param w2            first weight
*/
CombinedStemTropism::CombinedStemTropism(double n, double sigma,StemTropismFunction* t1, double w1, StemTropismFunction* t2, double w2): StemTropismFunction(n,sigma)
{
  tropisms.push_back(t1);
  tropisms.push_back(t2);
  weights.push_back(w1);
  weights.push_back(w2);
}

/**
 * getHeading() minimizes this function, @see TropismFunction::tropismObjective
 */
double CombinedStemTropism::stemtropismObjective(const Vector3d& pos, Matrix3d old, double a, double b, double dx, const Organ* stem)
{
  //std::cout << "CombinedTropsim::tropismObjective\n";
  double v = tropisms[0]->stemtropismObjective(pos,old,a,b,dx,stem)*weights[0];
  for (size_t i = 1; i< tropisms.size(); i++) {
    v += tropisms[i]->stemtropismObjective(pos,old,a,b,dx,stem)*weights[i];
  }
  return v;
}





} // namespace CPlantBox
