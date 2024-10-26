#ifndef _CPLANTBOX_CATMULLROMSPLINE_H
#define _CPLANTBOX_CATMULLROMSPLINE_H

#include "mymath.h"
#include "Quaternion.h"
#include <optional>
#include <vector>

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
    Vector3d quadp = 0.5 * ((2.0*y1) + (-y0 + y2) * t_ + (2.0*y0 - 5.0*y1 + 4.0*y2 - y3) * t_ * t_ + (-y0 + 3.0*y1 - 3.0*y2 + y3) * t_ * t_ * t_);
    Vector3d linp = y1 + t_ * (y2 - y1);
    return quadp * alpha + linp * (1.0 - alpha);
  }

  double getT0() const { return t0; }
  double getT1() const { return t1; }

  Vector3d derivative(double t) const {
    double t_ = (t-t0)/(t1-t0);
    Vector3d quadd = 0.5 * ((-y0 + y2) + (2.0*y0 - 5.0*y1 + 4.0*y2 - y3) * 2.0 * t_ + (-y0 + 3.0*y1 - 3.0*y2 + y3) * 3.0 * t_ * t_);
    Vector3d lind = y2 - y1;
    return quadd * alpha + lind * (1.0 - alpha);
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

  int closestId(double t) const {
    double t_ = (t-t0)/(t1-t0);
    if(t_ < 0.5)
    {
      return id1;
    }
    else
    {
      return id2;
    }
  }

  void IdPair(int id1, int id2) {
    this->id1 = id1;
    this->id2 = id2;
  }

  void setAlpha(double alpha) {
    this->alpha = alpha;
  }

  private:
  // start and end time of the spline
  double t0, t1;
  // the control points of the spline
  Vector3d y0, y1, y2, y3;
  // the original IDs of the control points
  int id1, id2;
  // spline parameter
  double alpha = 0.9;
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
  CatmullRomSplineManager(std::vector<Vector3d> y, std::optional<std::vector<int>> ids = std::nullopt) : y(y), indices(ids) {
    computeT();
  }
  CatmullRomSplineManager(std::initializer_list<Vector3d> y) : y(y), indices(std::nullopt) {
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

  double getT0() const {
    return t0;
  }
  double getT1() const {
    return t1;
  }

  std::vector<double> getT() const {
    return yt;
  }

  void setAlpha(double alpha, int i) {
    splines[i].setAlpha(alpha);
  }

  const std::vector<CatmullRomSpline> &getSplines() const
  {
    return splines;
  }

  Vector3d operator() (double t) const {
    return std::find_if(splines.begin(), splines.end(), [t](const CatmullRomSpline &s) { return t >= s.getT0() && t <= s.getT1(); })->operator()(t);
  }

  Vector3d derivative(double t) const {
    return std::find_if(splines.begin(), splines.end(), [t](const CatmullRomSpline &s) { return t >= s.getT0() && t <= s.getT1(); })->derivative(t);
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

  Vector3d help_lower() const {
    return y[0] - (y[1] - y[0]);
  }
  Vector3d help_upper() const {
    return y[y.size()-1] + (y[y.size()-1] - y[y.size()-2]);
  }
  Vector3d getControlPoint(int i) const {
    return y[i];
  }
  int size() const {
    return y.size();
  }
  int splineSize() const {
    return splines.size();
  }
  std::vector<double> getTValues(int spline) const {
    return {splines[spline].getT0(), splines[spline].getT1()};
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
    for(int i = 0; i < y.size()-1; i++)
    {
      if(i == 0)
      {
        auto helper = help_lower();
        splines.push_back(CatmullRomSpline({helper, y[0], y[1], y[2]}, 0.0, yt[1]));
        if(indices)
        {
          splines.back().IdPair(indices.value()[0], indices.value()[1]);
        }
      }
      else if(i == y.size()-1)
      {
        auto helper = help_upper();
        splines.push_back(CatmullRomSpline({y[i-1], y[i], y[i+1], helper}, yt[i], 1.0));
        if(indices)
        {
          splines.back().IdPair(indices.value()[i], indices.value()[i+1]);
        }
      }
      else
      {
        splines.push_back(CatmullRomSpline({y[i-1], y[i], y[i+1], y[i+2]}, yt[i], yt[i+1]));
        if(indices)
        {
          splines.back().IdPair(indices.value()[i], indices.value()[i+1]);
        }
      }
    }
    this->t0 = splines[0].getT0();
    this->t1 = splines.back().getT1();
  }

  std::vector<Vector3d> y; // control points
  std::vector<double> yt; // t values
  std::optional<std::vector<int>> indices;
  std::vector<CatmullRomSpline> splines;
  float t0 = -1.0;
  float t1 = -1.0;
};

}

#endif // CATMULLROMSPLINE_H
