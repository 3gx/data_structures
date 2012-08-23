#pragma once

#include <cassert>
#include <vector>
#include "vector3.h"


#if 0
struct HalfSpace
{
  typedef std::vector<HalfSpace> Vector;
  vec3 n;
  real h;
  HalfSpace() {}
  HalfSpace(const vec3 &_n, const vec3 &p) : n(_n)
  {
#if 0
    const real fn = n.abs();
    assert(fn > 0.0);
    n *= 1.0/fn;
#endif
    h  = p*n;
  }
  HalfSpace(const vec3 &p) : n(p)
  {
#if 0
    const real fn = n.abs();
    assert(fn >  0.0);
    n *= 1.0/fn;
#endif
    h = p*n;
  }


  bool outside(const vec3 &p) const 
  {
#if 0
    return n*p < h;
#else
    return n*p - h < -1.0e-10*std::abs(h);
#endif
  }
  friend std::pair<vec3, real> intersect
    (const HalfSpace &p1, const HalfSpace &p2, const HalfSpace &p3, const vec3 &c)
    {
      const vec3 w1 = p2.n%p3.n;
      const vec3 w2 = p3.n%p1.n;
      const vec3 w3 = p1.n%p2.n;
      const real  w = p1.n*w1;
      if (w == 0.0) return std::make_pair(0.0, -HUGE);
      const real iw = 1.0/w;
      const real d1 = p1.h * iw;
      const real d2 = p2.h * iw;
      const real d3 = p3.h * iw;
      const vec3 v  = w1*d1 + w2*d2 + w3*d3;
      return std::make_pair(v, c*v);
    }

#if 0
  real dist(const vec3 &p) const
  {
    return n*p - h;
  }
  vec3 project(const vec3 &p) const
  {
    return p - n*(p*n - h);
  }
#endif
};
#endif


template<const int N>
struct DirectPolyhedron
{
  typedef double real;
  typedef vector3<real> vec3;

  private:
#if 0
  struct HalfSpace
  {
    vec3 n;
    real h;
    HalfSpace() {}
    HalfSpace(const vec3 &_n, const vec3 &p) : n(_n)
    {
      assert(sizeof(HalfSpace) == 4*sizeof(real));
      h = p*n;
    }

    bool outside(const vec3 &p) const
    {
      return n*p< h;
    }
  };
#endif

  int n;
  HalfSpace list[N];
  int         nb[N];

  public:
  DirectPolyhedron() : n(0) {}
  void clear() { n = 0; }
  int  nface() const { return n; }
  int operator[](const int i) const { return nb[i]; }
  const HalfSpace& getHalfSpace(const int i) const { return list[i];}

  bool intersect(const vec3 &pos)
  {
    for (int i = 0; i < n; i++)
      if (list[i].outside(pos))
        return false;
    return true;
  }

  bool push(const vec3 &pos, const int idx, const real f = 1.0)
  {

    /* check if the current point is inside direct polyhedron */

    const vec3 hpos = pos * (1.0/f);
    for (int i = 0; i < n; i++)
      if (list[i].outside(hpos))
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
    list[n  ] = HalfSpace(-pos,pos);
    nb  [n++] = idx;

    return true;
  }
};

