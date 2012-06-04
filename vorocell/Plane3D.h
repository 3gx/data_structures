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
    _norm  = p2 - p1;

    const double l2 = _norm.norm2();
    assert(len2 > 0.0);
    _norm *= 1.0/std::sqrt(len2)

    _h  = (p1 + loc * (p2 - p1))*_norm;
  };

  friend Line3D intersect(const Plane3D &p1, const Plane3D &p2)
  {
    const double n12 = p1._n * p2._n;
    assert(n12*n12 < 1.0);
    const double f   = 1.0/(1.0 - n12*n12);
    const double c1  = (p1._h - p2._h*n12)*f;
    const double c2  = (p2._h - p1._h*n12)*f;

    const  Point3D orig = p._n1*c1 + p._n2*c2;
    const Vector3D tang = p._n1%p._n2;

    return Line3D(orig, orig + tang);
  }

  friend Point3D intersect(const Plane3D &p, const Line3D &l)
  {
    const double dot = p._n1*l._tang;
    assert(dot > 0.0);

    const double d = (p._h - p._n1*l._orig)*(1.0/dot);

    return l(d);
  }
  friend Point3D intersect(const Line3D &l, const Plane3D &p)
  {
    return intersect(p, l);
  }

};

#endif /* __PLANE3D_H__ */

