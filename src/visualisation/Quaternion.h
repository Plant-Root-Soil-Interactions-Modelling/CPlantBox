#ifndef _CPLANTBOX_QUATERNION_H_
#define _CPLANTBOX_QUATERNION_H_
#pragma once

#include "mymath.h"

#include <limits>

namespace CPlantBox {

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// convenience functions that are used in this code only

static inline Vector3d operator+ (const Vector3d& a, const Vector3d& b) {
  return Vector3d(a.x + b.x, a.y + b.y, a.z + b.z);
}
static inline Vector3d operator- (const Vector3d& a, const Vector3d& b) {
  return Vector3d(a.x - b.x, a.y - b.y, a.z - b.z);
}
static inline Vector3d operator* (const Vector3d& a, double s) {
  return Vector3d(a.x * s, a.y * s, a.z * s);
}
static inline Vector3d operator* (double s, const Vector3d& a) {
  return Vector3d(a.x * s, a.y * s, a.z * s);
}
static inline Vector3d operator/ (const Vector3d& a, double s) {
  return Vector3d(a.x / s, a.y / s, a.z / s);
}
static inline Vector3d operator- (const Vector3d& a) {
  return Vector3d(-a.x, -a.y, -a.z);
}
static inline Vector3d vectorNormalized(const Vector3d& a) {
  return a / a.length();
}

/**
 * Quaternion class with common functions for CG
 * this is for general Quaternions, but we use it for H0 and H1
*/
class Quaternion {
  public:

    double w{};
    Vector3d v{};

    Quaternion() = default;
    ~Quaternion() = default;
    Quaternion(double w, Vector3d v) : w(w), v(v) {}
    Quaternion(const Quaternion& q) : w(q.w), v(q.v) {}
    Quaternion(Quaternion&& q) : w(q.w), v(q.v) {}
    Quaternion(std::initializer_list<double> l) {
      assert(l.size() == 4);
      auto it = l.begin();
      w = *it;
      ++it;
      v = Vector3d(*it, *(it+1), *(it+2));
    }

    Quaternion& operator=(const Quaternion& q) { w = q.w; v = q.v; return *this; }

    Quaternion& operator=(Quaternion&& q) { w = q.w; v = q.v; return *this; }

    Quaternion operator+(const Quaternion& q) const { return Quaternion(w + q.w, v + q.v); }

    Quaternion operator-(const Quaternion& q) const { return Quaternion(w - q.w, v - q.v); }

    Quaternion operator*(const Quaternion& q) const {
      return Quaternion(w * q.w - v.times(q.v), q.v.times(w) + v.times(q.w) + v.cross(q.v));
      //return Quaternion(1.0, q.v.times(w) + v.times(q.w) + v.cross(q.v));
    }

    std::string toString() const {
      std::ostringstream strs;
      strs << w << " + " << v.toString();
      return strs.str();
    }

    friend Quaternion operator*(double s, const Quaternion& q) { return Quaternion(q.w * s, q.v * s); }

    Quaternion operator*(double s) const { return Quaternion(w * s, v * s); }

    Quaternion operator/(double s) const { return Quaternion(w / s, v / s); }

    Quaternion operator-() const { return Quaternion(-w, -v); }

    double& operator[](int i) {
      assert((i >= 0) && (i < 4));
      switch (i) {
        case 0: return w;
        case 1: return v.x;
        case 2: return v.y;
        case 3: return v.z;
      }
      throw 0; // just to not produce a warning
    }

    Quaternion& operator+=(const Quaternion& q) {
      w += q.w;
      v = v + q.v;
      return *this;
    }

    Quaternion& operator-=(const Quaternion& q) {
      w -= q.w;
      v = v - q.v;
      return *this;
    }

    Quaternion& operator*=(const Quaternion& q) {
      w = w * q.w - v.times(q.v);
      v = w * q.v + v * q.w + v.cross(q.v);
      return *this;
    }

    Quaternion& operator*=(double s) {
      w *= s;
      v = v * s;
      return *this;
    }

    Quaternion& operator/=(double s) {
      w /= s;
      v = v / s;
      return *this;
    }

    Quaternion conjugate() const {
      return Quaternion(w, -v);
    }

    Quaternion inverse() const {
      return conjugate() / norm2();
    }

    Quaternion inverse2() const {
      return conjugate() / (w * w + v.times(v));
    }

    void invert() {
      *this = inverse();
    }
    void conjugateInPlace() {
      v = -v;
    }

    inline double norm() const {
      return std::sqrt(w * w + v.times(v));
    }

    inline double norm2() const {
      return w * w + v.times(v);
    }

    inline Quaternion normalized() const {
      return *this / norm();
    }

    void normalize() {
      *this /= norm();
    }

    double dot(const Quaternion& q) const {
      return w * q.w + v.times(q.v);
    }

    // a method that calculates a quaternion from a forward vector.
    static Quaternion FromForward(const Vector3d& forward) {
      Vector3d up(0, 0, 1);
      Vector3d right = up.cross(forward);
      up = forward.cross(right);
      right.normalize();
      up.normalize();
      return Quaternion::FromMatrix3d(Matrix3d(forward, right, up));
    }

    // a method that calculates a quaternion from a forward vector while respecting the previous up vector.
    static Quaternion FromForwardAndUp(const Vector3d& forward, const Vector3d& prev_up)
    {
      // initial copy of the previous up vector
      Vector3d up = vectorNormalized(prev_up);
      // calculate the right vector
      Vector3d right = up.cross(vectorNormalized(forward));
      // calculate the new up vector
      up = forward.cross(right);
      // normalize the vectors
      right.normalize();
      up.normalize();
      // return the quaternion as calculated from the matrix
      return Quaternion::FromMatrix3d(Matrix3d(vectorNormalized(forward), right, up));
    }

    static Quaternion FromForwardAndRight(const Vector3d& forward, const Vector3d& right)
    {
      // calculate the new up vector
      Vector3d up = vectorNormalized(forward).cross(vectorNormalized(right));
      // normalize the vectors
      up.normalize();
      // return the quaternion as calculated from the matrix
      return Quaternion::FromMatrix3d({forward, right, up});
    }


    /**
     * @brief Rotates a vector by the quaternion
     * @param v vector to rotate
     * @return rotated vector
     * @note the formula is v' = q * v * q^-1
    */
    Vector3d Rotate(const Vector3d& v) const {
      Quaternion standin = *this;
      standin.w = 1.0;
      auto result = (standin * Quaternion(0, v) * standin.conjugate()).v;
      return result / result.length() * v.length();
    }

    /**
     * @brief A method that rotates a vector, ensuring uniform scale
     * @param v vector to rotate
     * @return rotated vector
     */
    Vector3d RotateUniform(const Vector3d& v) const {
      Quaternion standin = *this;
      standin.w = 1.0;
      return (standin * Quaternion(0, v) * standin.conjugate()).v * norm();
    }

    /**
     * @brief Converts the quaternion to a 3x3 rotation matrix
     * @return 3x3 rotation matrix
     * @note this is essentially the same as rotating the axis vectors by the quaternion
    */
    Matrix3d ToMatrix3d() const {
      return {{1 - 2 * v.y * v.y - 2 * v.z * v.z, 2 * v.x * v.y - 2 * w * v.z, 2 * v.x * v.z + 2 * w * v.y},
      {2 * v.x * v.y + 2 * w * v.z, 1 - 2 * v.x * v.x - 2 * v.z * v.z, 2 * v.y * v.z - 2 * w * v.x},
      {2 * v.x * v.z - 2 * w * v.y, 2 * v.y * v.z + 2 * w * v.x, 1 - 2 * v.x * v.x - 2 * v.y * v.y}};
    }

    static Quaternion FromMatrix3d(const Matrix3d& m)
    {
      double w = std::sqrt(1 + m.r0.x + m.r1.y + m.r2.z) / 2.0;
      double w4 = (4 * w);
      double x = (m.r0.y - m.r1.z) / w4;
      double y = (m.r0.z - m.r2.x) / w4;
      double z = (m.r1.x - m.r0.y) / w4;
      return Quaternion(w, Vector3d(x, y, z)).normalized();
    }

    // a method that returns a look-at-direction representing the local x axis
    // TODO: PLEASE check with daniel about coordinate axis conventions!!
    Vector3d Forward() const 
    {
      return Rotate(Vector3d(1,0,0));
    }

    // a method that returns a look-at-direction representing the local y axis
    // TODO: SOON check with daniel about coordinate axis conventions!!
    Vector3d Right() const 
    {
      return Rotate(Vector3d(0,1,0));
    }

    // a method that returns a look-at-direction representing the local z axis
    // TODO: NOW check with daniel about coordinate axis conventions!!
    Vector3d Up() const 
    {
      return Rotate(Vector3d(0,0,1));
    }

    // a method that computes the normed linear sum of axis vectors
    // as I always use it with a set length, I am not normalizing the result
    Vector3d Axis(int x = 0, int y = 0, int z = 0) const
    {
      return (Forward() * x + Right() * y + Up() * z);
    }

    // a method that selects the shortest arc between two quaternions
    static Quaternion ShortestArc(const Quaternion& a, const Quaternion& b)
    {
      if (a.dot(b) < 0)
        return a * -1.0;
      else
        return a;
    }

    /**
     * @brief Static method to compute look-at-direction between two points
     * @param from the point to look from
     * @param to the point to look at
     * @param up the up direction
    */
    static Quaternion LookAt(const Vector3d& from, const Vector3d& to, const Vector3d& up)
    {
      // compute the forward direction
      auto forward = vectorNormalized(to - from);
      // compute the right direction
      auto right = vectorNormalized(forward.cross(up));
      // compute the up direction
      auto up2 = vectorNormalized(right.cross(forward));
      // compute the rotation matrix
      Matrix3d m = {-forward, right, up2};
      // return the quaternion
      return FromMatrix3d(m);
    }

    // This Method calculates the Quaternion that rotates a to b
    // this is basically the same as the implementation in the Unreal Engine
    inline static Quaternion geodesicRotation(Vector3d a, Vector3d b)
    {
      double w = 1.0;
      // normalize the vectors
      a.normalize();
      b.normalize();
      // if the vectors are parallel, the rotation axis is undefined
      if(w < std::numeric_limits<double>::epsilon() && std::abs(a.x) <= std::abs(a.z))
        return Quaternion(w, Vector3d(0, -a.z, a.y));
      // case 2 of the above
      else if(w < std::numeric_limits<double>::epsilon() && std::abs(a.x) > std::abs(a.z))
        return Quaternion(w, Vector3d(-a.y, a.x, 0));
      // otherwise, the rotation axis is the cross product of a and b
      else
        return Quaternion(w, a.cross(b)).normalized();
    }


    // this method provides a spherical interpolation between two quaternions
    inline static Quaternion SphericalInterpolation(Quaternion a, Quaternion b, double t = 0.5)
    {
      double angle = a.dot(b);
      if (angle >= 1.0) ///< if angle is 1.0, the quaternions are the same
      {
        return a;
      }
      double half = acos(angle);
      double sinHalf = sin(half);
      if(std::abs(sinHalf) < std::numeric_limits<double>::epsilon()) //< if angle is 0.0, the quaternions are opposite
      {
        return {0.5*a.w + 0.5*b.w, 0.5*a.v + 0.5*b.v};
      }
      //< otherwise, interpolate
      double ratioA = sin((1.0 - t) * half) / sinHalf;
      double ratioB = sin(t * half) / sinHalf;
      return {ratioA * a.w + ratioB * b.w, ratioA * a.v + ratioB * b.v};
    }
};

} // namespace CPlantBox

#endif // _CPLANTBOX_QUATERNION_H_
