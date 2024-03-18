// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "tropism.h"

#include "soil.h"
#include "sdf.h"
#include "Root.h"
#include <stdlib.h>
#include <cstdlib>

namespace CPlantBox {

/**
 * Copies this tropism
 */
std::shared_ptr<Tropism> Tropism::copy(std::shared_ptr<Organism> plant)
{
    auto nt = std::make_shared<Tropism>(*this); // default copy constructor
    nt->plant = plant;
    return nt;
}

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
Vector3d Tropism::getPosition(const Vector3d& pos, const Matrix3d& old, double a, double b, double dx)
{
    return pos.plus((old.times(Vector3d::rotAB(a,b))).times(dx));
}

/**
 * Dices N times picking angles alpha and beta, takes the optimal direction according to the objective function
 *
 * @param pos          root tip position
 * @param old          rotation matrix, heading is old(:,1)
 * @param dx           distance to look ahead (e.g in case of hydrotropism)
 * @param root         points to the root that called getHeading(...), just in case something else is needed (i.e. iheading for exotropism)
 * @param nodeIdx    index of the node on which to apply the tropism
 *
 * \return             angle alpha and beta
 */
Vector2d Tropism::getUCHeading(const Vector3d& pos, const Matrix3d& old, double dx, const std::shared_ptr<Organ> o, int nodeIdx)
{
    double a = sigma*randn(nodeIdx)*sqrt(dx);
    double b = rand(nodeIdx)*2*M_PI;
    double v;

    double n_=n*sqrt(dx);
    if (n_>0) {
        double dn = n_-floor(n_);
        if (rand(nodeIdx)<dn) {
            n_ = ceil(n_);
        } else {
            n_ = floor(n_);
        }
        double bestA = a;
        double bestB = b;
        double bestV = this->tropismObjective(pos,old,a,b,dx,o);
        for (int i=0; i<n_; i++) {
            b = rand(nodeIdx)*2*M_PI;
            a = sigma*randn(nodeIdx)*sqrt(dx);
            v = this->tropismObjective(pos,old,a,b,dx,o);
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
 * Inside the geometric domain baseTropism::getHeading() is returned.
 * In case geometric boundaries are hit, the rotations are modified, so that growth stays inside the domain.
 *
 * @param pos        root tip position
 * @param old        rotation matrix, heading is old(:,1)
 * @param dx         distance to look ahead (e.g in case of hydrotropism)
 * @param root       points to the root that called getHeading(), just in case something else is needed (i.e. iheading for exotropism)
 * @param nodeIdx    index of the node on which to apply the tropism
 *
 * \return           the rotations alpha and beta
 */
Vector2d Tropism::getHeading(const Vector3d& pos, const Matrix3d& old, double dx, const std::shared_ptr<Organ> o, int nodeIdx)
{
    if(nodeIdx > 0 ){gen =  std::mt19937(plant.lock()->getSeedVal() + nodeIdx + o->getId());}
    Vector2d h = this->getUCHeading(pos, old, dx, o, nodeIdx);
    double a = h.x;
    double b = h.y;

    if (!geometry.expired()) {
        double d = geometry.lock()->getDist(this->getPosition(pos,old,a,b,dx));
        double dmin = d;

        double bestA = a;
        double bestB = b;
        int i=0; // counts change in alpha
        int j=0;    // counts change in beta

        while ((d>0)&&(o->organType()== Organism::ot_root))  { // not valid
            i++;
            j=0;
            while ((d>0) && j<betaN) { // change beta

                b = 2*M_PI*rand(nodeIdx); // dice
                d = geometry.lock()->getDist(this->getPosition(pos,old,a,b,dx));
                if (d<dmin) {
                    dmin = d;
                    bestA = a;
                    bestB = b;
                }
                j++;
            }

            if ((d>0)&&(o->organType()== Organism::ot_root))  {
                a = a + M_PI/2./double(alphaN);
            }

            if (i>alphaN) {
                //std::cout << "Could not respect geometry boundaries \n";
                a = bestA;
                b = bestB;
                break;
            }

        }
    }
    return Vector2d(a,b);
}



/**
 * getHeading() minimizes this function, @see TropismFunction::tropismObjective
 */
double Exotropism::tropismObjective(const Vector3d& pos, const Matrix3d& old, double a, double b, double dx, const std::shared_ptr<Organ> o)
{
    Vector3d iheading =o->getiHeading0();
    double s = iheading.times(old.times(Vector3d::rotAB(a,b)));
    s*=(1./iheading.length()); // iheading should be normed anyway?
    s*=(1./old.column(0).length());
    return acos(s)/M_PI; // 0..1
}



/**
 * getHeading() minimizes this function, @see TropismFunction::tropismObjectiveMatrix3d
 */
double Hydrotropism::tropismObjective(const Vector3d& pos, const Matrix3d& old, double a, double b, double dx, const std::shared_ptr<Organ> o)
{
    assert(!soil.expired());
    Vector3d newpos = this->getPosition(pos,old,a,b,dx);
    double v = soil.lock()->getValue(newpos,o);
    // std::cout << "\n" << newpos.getString() << ", = "<< v;
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
CombinedTropism::CombinedTropism(std::shared_ptr<Organism> plant, double n, double sigma, std::shared_ptr<Tropism> t1, double w1, std::shared_ptr<Tropism> t2, double w2)
    :Tropism(plant,n,sigma)
{
    tropisms.push_back(t1);
    tropisms.push_back(t2);
    weights.push_back(w1);
    weights.push_back(w2);
}

/**
 * getHeading() minimizes this function, @see TropismFunction::tropismObjective
 */
double CombinedTropism::tropismObjective(const Vector3d& pos, const Matrix3d& old, double a, double b, double dx, const std::shared_ptr<Organ> o)
{
    //std::cout << "CombinedTropsim::tropismObjective\n";
    double v = tropisms[0]->tropismObjective(pos,old,a,b,dx,o)*weights[0];
    for (size_t i = 1; i< tropisms.size(); i++) {
        v += tropisms[i]->tropismObjective(pos,old,a,b,dx,o)*weights[i];
    }
    return v;
}

} // end namespace CPlantBox
