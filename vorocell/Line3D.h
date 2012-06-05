#ifndef __LINE3D_H__
#define __LINE3D_H__

#include "Point3D.h"

struct Plane3D;
struct Line3D
{
  typedef std::vector<Line3D> Vector;
  
  private:
  Vector3D _tang;
   Point3D _orig;

  public:
  Line3D() {}
  Line3D(const Point3D &p1, const Point3D &p2)
  {
    _orig = p1;
    const Point3D diff = p2 - p1;
    _tang = diff;

    const double len2 = _tang.norm2();
    assert(len2 > 0.0);
    _tang *= 1.0/std::sqrt(len2);

  }

  const Vector3D& tang() const {return _tang;}
  const  Point3D& orig() const {return _orig;}

  Point3D operator()(const double lambda) const
  {
    return Point3D(_orig + lambda*_tang);
  }

  friend Point3D perpendicularFoot(const Line3D &l, const Point3D &b)
  {
    const double lambda = (b - l._orig)*l._tang;
    return l(lambda);
  }
  friend Point3D perpendicularFoot(const Point3D &p, const Line3D &l)
  {
    return perpendicularFoot(l, p);
  }

};

#endif  /* __LINE3D_H__ */
