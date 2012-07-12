#ifndef __SEIDELLP_H__
#define __SEIDELLP_H__

#include "vector3.h"
#include <vector>
#include <cassert>

struct HalfSpace
{
  typedef std::vector<HalfSpace> Vector;
  vec3 n;
  real h;
  HalfSpace() {}
  HalfSpace(const vec3 &_n, const vec3 &p) : n(_n)
  {
    const real fn = n.abs();
    assert(fn > 0.0);
    n *= 1.0/fn;
    h  = p*n;
  }
  bool outside(const vec3 &p) const 
  {
    return n*p < h;
  }
  friend vec3 intersect(const HalfSpace &p1, const HalfSpace &p2, const HalfSpace &p3)
  {
    const vec3 w1 = p2.n%p3.n;
    const vec3 w2 = p3.n%p1.n;
    const vec3 w3 = p1.n%p2.n;
    const real  w = p1.n*w1;
    if (w == 0.0) return vec3(HUGE);
    assert(w != 0.0);
    const real iw = 1.0/w;
    const real d1 = p1.h * iw;
    const real d2 = p2.h * iw;
    const real d3 = p3.h * iw;
    return w1*d1 + w2*d2 + w3*d3;
  }
};


template<const int N>
struct SeidelLP
{
  private:
    int n;
    HalfSpace halfSpaceList[N];
    int n1, n2, n3;

  public:
    SeidelLP() : n(0) {}
    SeidelLP(const HalfSpace::Vector &plist) 
    {
      assert((int)plist.size() <= N);
      n = plist.size();
      for (int i = 0; i < n; i++)
        halfSpaceList[i] = plist[i];
    }

    void clear() { n = 0; }
    bool push(const HalfSpace &p)
    {
      assert(n < N);
      halfSpaceList[n++] = p;
      return true;
    }
    int nspace() const { return n; }

    void rotate(const vec3 &line)
    {
    }

    vec3 solve()
    {
      n1 = n2 = n3 = 0;
      std::random_shuffle(halfSpaceList, halfSpaceList+n);
      const vec3 v = solve_lp3D(n);

      fprintf(stderr, "n3= %d\n", n3);
      fprintf(stderr, "n2= %d\n", n2);
      fprintf(stderr, "n1= %d\n", n1);

      return v;
    }

    vec3 solve_lp3D(const int n)
    {
      vec3 v(-HUGE);
      for (int i = 0; i < n; i++)
      {
        n3++;
        const HalfSpace &h = halfSpaceList[i];
        if (h.outside(v))
          v = solve_lp2D(i, v, h);
      }
      return v;
    }

    vec3 solve_lp2D(const int n, vec3 v, const HalfSpace h1)
    {
      for (int i = 0; i < n; i++)
      {
        n2++;
        const HalfSpace &h = halfSpaceList[i];
        if (h.outside(v))
          v = solve_lp1D(i, v, h1, h);
      }
      return v;
    }

    vec3 solve_lp1D(const int n, vec3 v, const HalfSpace h1, const HalfSpace h2)
    {
      for (int i = 0; i < n; i++)
      {
        n1++;
        const HalfSpace &h = halfSpaceList[i];
        if (h.outside(v))
          v = intersect(h, h1, h2);
      }
      return v;
    }

};

#endif /* __SEIDELLP_H__ */
