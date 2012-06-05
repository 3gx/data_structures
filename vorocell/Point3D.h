#ifndef __POINT3D_H__
#define __POINT3D_H__

#include <vector>
#include <deque>
#include <cassert>
#include <cmath>
#include "vector3.h"

struct Point3D : public vector3<double> 
{
  typedef std::vector<Point3D> Vector;
  operator vector3<double>() const {return (vector3<double>)*this;}
  Point3D() : vector3<double>() {}
  Point3D(const vector3<double> &p) : vector3<double>(p) {}
};


struct Vector3D : public Point3D
{
  typedef std::vector<Vector3D> Vector;
  operator Point3D() const {return (Point3D)*this;}
  Vector3D() : Point3D() {}
  Vector3D(const Point3D &p) : Point3D(p) {}
  Vector3D(const vector3<double> &p) : Point3D(p) {}
};


#endif /* __POINT3D_H__ */

