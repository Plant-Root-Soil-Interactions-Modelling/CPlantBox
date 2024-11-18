#ifndef _CPLANTBOX_CATMULLROMSPLINE_H
#define _CPLANTBOX_CATMULLROMSPLINE_H

#include "mymath.h"
#include "Quaternion.h"

namespace CPlantBox {


/**
 * Catmull-Rom spline interpolation
 * It is used to store the Catmull-Rom splines of depth 4
*/
class CatmullRomSpline
{
  public:
  CatmullRomSpline() = default;
  CatmullRomSpline(Vector3d y0, Vector3d y1, Vector3d y2, Vector3d y3, double t0, double t1) : y0(y0), y1(y1), y2(y2), y3(y3), t0(t0), t1(t1) {}
  CatmullRomSpline(std::vector<Vector3d> y, double t0, double t1) : y0(y[0]), y1(y[1]), y2(y[2]), y3(y[3]), t0(t0), t1(t1) {}
  Vector3d operator() (double t) const {
    double t_ = (t-t0)/(t1-t0);
    return 0.5 * ((2.0*y1) + (-y0 + y2) * t_ + (2.0*y0 - 5.0*y1 + 4.0*y2 - y3) * t_ * t_ + (-y0 + 3.0*y1 - 3.0*y2 + y3) * t_ * t_ * t_);
  }

  double getT0() const { return t0; }
  double getT1() const { return t1; }

  Vector3d derivative(double t) const {
    double t_ = (t-t0)/(t1-t0);
    return 0.5 * ((-y0 + y2) + (2.0*y0 - 5.0*y1 + 4.0*y2 - y3) * 2.0 * t_ + (-y0 + 3.0*y1 - 3.0*y2 + y3) * 3.0 * t_ * t_);
  }

  Vector3d operator[](int i) {
    switch(i) {
      case 0: return y0;
      case 1: return y1;
      case 2: return y2;
      case 3: return y3;
      default: throw std::runtime_error("CatmullRomSpline: index out of bounds");
    }
  }

  Quaternion computeOrientation(double t) const {
    Vector3d v = derivative(t);
    return Quaternion::FromForward(v);
  }

  private:
  // the control points of the spline
  Vector3d y0, y1, y2, y3;
  // start and end time of the spline
  double t0, t1;
  // the stages of the apline
  Vector3d a0, a1, a2, a3;
  Vector3d b0, b1, b2;
  // spline parameter
  double alpha = 0.5;
};

/**
 * A manager class for the Catmull-Rom spline interpolation
 * It is used to store the Catmull-Rom splines of depth 4
 * We use a weighted average of the splines to get a smooth transition
*/
class CatmullRomSplineManager
{
  public:
  CatmullRomSplineManager() = default;
  CatmullRomSplineManager(std::vector<Vector3d> y) : y(y) {
    computeT();
  }

  CatmullRomSpline spline(int i) const {
    return splines[i];
  }
  
  // a method that selects the most suitable spline depending on most equal distance to start and finish
  CatmullRomSpline selectSpline(double t) const {
    double min = std::numeric_limits<double>::max();
    int index = 0;
    for(int i = 0; i < splines.size(); i++)
    {
      double d = std::abs((splines[i].getT0() - t) - (splines[i].getT1() - t));
      if(d < min && t > splines[i].getT0() && t < splines[i].getT1())
      {
        min = d;
        index = i;
      }
    }
    return splines[index];
  }

  const std::vector<CatmullRomSpline> &getSplines() const
  {
    return splines;
  }

  Vector3d operator() (double t) const {
    Vector3d p(0,0,0);
    int sum = 0;
    for(int i = 0; i < splines.size(); i++)
    {
      // whether the spline has t in its interval
      bool in = t >= splines[i].getT0() && t <= splines[i].getT1();
      // we add the spline if it is in the interval
      if(in)
      {
        p = p + splines[i](t);
        sum++;
      }
    }
    return p / static_cast<double>(sum);
  }
  void setY(std::vector<Vector3d> y) {
    this->y = y;
    computeT();
  }

  // a method that returns the offset of the value to the last point in local space
  double getOffset(double t) const {
    int i = 0;
    while(i < yt.size() && yt[i] < t) i++;
    if(i + 1 >= yt.size())
    {
      return 1.0;
    }
    else
    {
      return (t-yt[i])/(yt[i+1]-yt[i]);
    }
    
  }

  // a method that returns the last index for t 
  int getIndex(double t) const {
    int i = 0;
    while(i < yt.size() && yt[i] < t) i++;
    return i;
  }

  private:

  void computeT()
  {
    yt.clear();
    yt.push_back(0);
    for(int i = 1; i < y.size(); i++)
    {
      yt.push_back(yt[i-1] + (y[i]-y[i-1]).length());
    }
    for(int i = 0; i < yt.size(); i++)
    {
      yt[i] /= yt.back();
    }
    splines.clear();
    for(int i = 0; i < y.size()-3; i++)
    {
      splines.push_back(CatmullRomSpline({y[i], y[i+1], y[i+2], y[i+3]}, yt[i], yt[i+3]));
    }
    this->t0 = splines[0].getT0();
    this->t1 = splines.back().getT1();
  }

  std::vector<Vector3d> y; // control points
  std::vector<double> yt; // t values
  std::vector<CatmullRomSpline> splines;
  float t0 = -1.0;
  float t1 = -1.0;
};

}

#endif // CATMULLROMSPLINE_H
