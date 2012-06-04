#ifndef __POINT3D_H__
#define __POINT3D_H__

#include <vector>
#include <cassert>
#include "vector3.h"

struct Point3D : public vector3<double> 
{
  typedef std::vector<Point3D> Vector;
};


struct Vector3D : public vector3<double> 
{
  typedef std::vector<Vector3D> Vector;
};


#endif /* __POINT3D_H__ */

