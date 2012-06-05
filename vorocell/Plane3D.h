#ifndef __PLANE3D_H__
#define __PLANE3D_H__

#include "Line3D.h"

struct Plane3D
{
  typedef std::vector<Plane3D> Vector;

  private:        /* plane equations, n.(r-r0)=0 -> n.r = h, where h = n.r0 */
  Vector3D _n;    /* n = unit normal */
  double   _h;    /* h = n.r0 */

  public:
  Plane3D() {}
  Plane3D(const Point3D &p1, const Point3D &p2, const double loc = 0.5)
  {
    assert(loc >= 0.0 && loc <= 1.0);
    _n = p2 - p1;

    const double l2 = _n.norm2();
    assert(l2 > 0.0);
    _n*= 1.0/std::sqrt(l2);

    _h  = (p1 + loc * (p2 - p1))*_n;
  };

  /* intersection of two planes, returns a line */
  friend Line3D intersect(const Plane3D &p1, const Plane3D &p2)
  {
    const double n12 = p1._n * p2._n;
    assert(n12*n12 < 1.0);
    const double f   = 1.0/(1.0 - n12*n12);
    const double c1  = (p1._h - p2._h*n12)*f;
    const double c2  = (p2._h - p1._h*n12)*f;

    const  Point3D orig = p1._n*c1 + p2._n*c2;
    const Vector3D tang = p1._n%p2._n;

    return Line3D(orig, orig + tang);
  }

  /* intersection of a plane and a line, returns a point */
  friend Point3D intersect(const Plane3D &p, const Line3D &l)
  {
    const double dot = p._n*l.tang();
    assert(dot > 0.0);

    const double d = (p._h - p._n*l.orig())*(1.0/dot);

    return l(d);
  }
  friend Point3D intersect(const Line3D &l, const Plane3D &p)
  {
    return intersect(p, l);
  }

  /* computes distance between a plane and a line */
  friend double distance(const Plane3D &plane, const Point3D &point)
  {
    return plane._h - plane._n*point;
  }
  friend double distance(const Point3D &point, const Plane3D &plane)
  {
    return distance(plane, point);
  }

};

#endif /* __PLANE3D_H__ */

