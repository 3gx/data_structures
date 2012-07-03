#ifndef __DIRECTPOLY_H__
#define __DIRECTPOLY_H__

#include <cassert>
#include "vector3.h"

template<const int N>
struct DirectPolyhedron
{
  typedef double real;
  typedef vector3<real> vec3;

  private:
    struct HalfSpace
    {
      vec3 p;
      real h;
      HalfSpace() {}
      HalfSpace(const vec3 &pos) : p(pos), h(pos*pos)
      {
        assert(sizeof(HalfSpace) == 4*sizeof(real));
      }

      bool inside(const vec3 &pos) const
      {
        return p*pos < h;
      }
    };

    int n;
    HalfSpace list[N];
    int         nb[N];

  public:
  DirectPolyhedron() : n(0) {}
  void clear() { n = 0; }
  int  nface() const { return n; }
  int operator[](const int i) const { return nb[i]; }

  bool push(const vec3 &pos, const int idx)
  {

    /* check if the current point is inside direct polyhedron */

    for (int i = 0; i < n; i++)
      if (!list[i].inside(pos))
        return false;

    /* if so, find out if there are reduntant points */

    const HalfSpace h(pos);
    int i = -1;
    while (++i < n)
      if (!h.inside(list[i].p))
      {
        std::swap(list[i  ], list[--n]);
        std::swap(  nb[i--],   nb[  n]);
      }

    /* now store the new point */
    assert(n < N);
    list[n  ] = h;
    nb  [n++] = idx;

    return true;
  }
};

#endif /* __DIRECTPOLY_H__ */

