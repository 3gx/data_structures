#pragma once

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
      HalfSpace(const vec3 &pos) : p(pos), h(p*p)
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

  bool push(const vec3 &pos, const int idx, const real f = 1.0)
  {

    /* check if the current point is inside direct polyhedron */

    const vec3 hpos = pos * (1.0/f);
    for (int i = 0; i < n; i++)
      if (!list[i].inside(hpos))
        return false;

    /* if so, find out if there are reduntant points */

#if 0
    const HalfSpace h(pos*f);
    int i = -1;
    while (++i < n)
      if (!h.inside(list[i].p))
      {
#if 0
        if (nb[i] == 20850 || nb[i] == 42255 || nb[i] == 0)
        {
          fprintf(stderr, "i= %d: pi= %g %g %g  | %g\n", nb[i], list[i].p.x, list[i].p.y, list[i].p.z, list[i].p.norm2());
          fprintf(stderr, "j= %d: pj= %g %g %g  | %g\n", idx, pos.x, pos.y, pos.z, pos.norm2());
          fprintf(stderr, "%g\n", pos*list[i].p);
        }
        //    assert(nb[i] != 20850);
#endif
        std::swap(list[i  ], list[--n]);
        std::swap(  nb[i--],   nb[  n]);
      }
#endif

    /* now store the new point */
    assert(n < N);
    list[n  ] = HalfSpace(pos);
    nb  [n++] = idx;

    return true;
  }
};

