#ifndef __SEIDELLP_H__
#define __SEIDELLP_H__

#include "vector3.h"
#include <vector>
#include <algorithm>
#include <cassert>

typedef double real;
typedef vector3<real> vec3;

template<class T>
inline T __min(const T a, const T b) {return a < b ? a : b;}

template<class T>
inline T __max(const T a, const T b) {return a > b ? a : b;}

template<class T>
inline T __abs(const T a) {return a < T(0.0) ? -a : a;}

template<class T>
inline T __sign(const T a) {return a < T(0.0) ? (T)-1.0 : (T)+1.0;}

unsigned long long nflops = 0;

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
  HalfSpace(const vec3 &p) 
  {
#if 0
    const real fn = p.abs();
    assert(fn >  0.0);
    n = p*(1.0/fn);
#else
    n = p;
#endif
    h = p*n;
  }


  bool outside(const vec3 &p) const __attribute__((always_inline))
  {
    nflops += 6;
    return n*p < h;
  }
  friend vec3 intersect(const HalfSpace &p1, const HalfSpace &p2, const HalfSpace &p3) __attribute__((always_inline))
  {
    nflops += 9*3 + 5+1+3*3+3*3+3*2;
    const vec3 w1 = p2.n%p3.n;
    const vec3 w2 = p3.n%p1.n;
    const vec3 w3 = p1.n%p2.n;
    const real  w = p1.n*w1;
//    if (w == 0.0) return vec3(HUGE);
    assert(w != 0.0);
    const real iw = 1.0/w;
    const real d1 = p1.h * iw;
    const real d2 = p2.h * iw;
    const real d3 = p3.h * iw;
    return w1*d1 + w2*d2 + w3*d3;
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


template<const int N>
struct SeidelLP
{
  private:
    int n;
    vec3 bmax, cvec;
    HalfSpace halfSpaceList[N];
    HalfSpace bnd[3];

  public:
    SeidelLP(const vec3 &_bmax = 1.0) : n(0), bmax(_bmax) {}
    SeidelLP(const HalfSpace::Vector &plist, const vec3 &_bmax = 1.0) : bmax(_bmax)
    {
      assert((int)plist.size() <= N);
      n = plist.size();
      for (int i = 0; i < n; i++)
        halfSpaceList[i] = plist[i];
    }

    void clear(const vec3 &_bmax = 1.0) { n = 0; bmax = _bmax; }
    bool push(const HalfSpace &p)
    {
      assert(n < N);
      halfSpaceList[n++] = p;
      return true;
    }
    int nspace() const { return n; }

    vec3 solve(const vec3 &cvec, const bool randomize = false)
    {
      this->cvec = cvec;
      bnd[0] = HalfSpace(vec3(cvec.x > 0.0 ? bmax.x : -bmax.x, 0.0, 0.0));
      bnd[1] = HalfSpace(vec3(0.0, cvec.y > 0.0 ? bmax.y : -bmax.y, 0.0));
      bnd[2] = HalfSpace(vec3(0.0, 0.0, cvec.z > 0.0 ? bmax.z : -bmax.z));

      if (randomize)
        std::random_shuffle(halfSpaceList, halfSpaceList+n);

      return solve_lp3D(n);
    }

    inline vec3 solve_lp3D(const int n) const 
    {
      asm("#test1");
      vec3 v = intersect(bnd[0], bnd[1], bnd[2]);
      for (int i = 0; i < n; i++)
      {
        const HalfSpace &h = halfSpaceList[i];
        if (h.outside(v)) 
          v = solve_lp2D(i, h);
      }
      asm("#test2");
      return v;
    }

    inline vec3 solve_lp2D(const int n, const HalfSpace &h1) const __attribute__((always_inline))
    {
      const vec3 v1 = intersect(h1, bnd[0], bnd[1]);
      const vec3 v2 = intersect(h1, bnd[0], bnd[2]);
      const vec3 v3 = intersect(h1, bnd[1], bnd[2]);
      const real f1 = cvec*v1;
      const real f2 = cvec*v2;
      const real f3 = cvec*v3;
      vec3 v = v1;
      real f = f1;
      if (f2 > f) {f = f2; v = v2;}
      if (f3 > f) {f = f3; v = v3;}

      for (int i = 0; i < n; i++)
      {
        const HalfSpace &h = halfSpaceList[i];
        if (h.outside(v)) 
          v = solve_lp1D(i,  h1, h);
      }
      return v;
    }

#if 0
    inline vec3 solve_lp1D(const int n, const HalfSpace &h1, const HalfSpace &h2) const
    {
      const real norm2 = h1.n.norm2() * h2.n.norm2();
      const real n12   = h1.n * h2.n  * (1.0/std::sqrt(norm2));
//      if (n12*n12 == 1.0) return v;
      assert(n12*n12 < 1.0);
      const real f   = 1.0/(1.0 - n12*n12);
      const real c1  = (h1.h - h2.h*n12)*f;
      const real c2  = (h2.h - h1.h*n12)*f;

      const vec3 orig = h1.n*c1 + h2.n*c2;
      const vec3 tang = h1.n%h2.n;

      real tmin = -HUGE;
      real tmax = +HUGE;

      for (int i = 0; i < n; i++)
      {
        const HalfSpace &h = halfSpaceList[i];
        const real     dot = h.n * tang;
        if (dot*dot == 0.0) continue;
        assert(dot*dot > 0.0);
        const real tau = (h.h - h.n*orig)*(1.0/dot);
        if (dot > 0.0) tmin = __max(tmin, tau);
        else           tmax = __min(tmax, tau);
      }

      assert(tang.z != 0.0);
    
      const vec3 v = orig + tang * (tang*cvec > 0.0 ? tmax : tmin); 
      return v;
    }
#else
    inline vec3 solve_lp1D(const int n, const HalfSpace &h1, const HalfSpace &h2) const __attribute__((always_inline))
    {
      const vec3 v1 = intersect(h1, h2, bnd[0]);
      const vec3 v2 = intersect(h1, h2, bnd[1]);
      const vec3 v3 = intersect(h1, h2, bnd[2]);
      const real f1 = cvec*v1;
      const real f2 = cvec*v2;
      const real f3 = cvec*v3;
      
      vec3 v = v1;
      real f = f1;
      if (f2 > f) {f = f2; v = v2;}
      if (f3 > f) {f = f3; v = v3;}
      
      for (int i = 0; i < n; i++)
      {
        const HalfSpace &h = halfSpaceList[i];
        if (h.outside(v)) 
          v = intersect(h1, h2, h);
      }
      return v;
    }
#endif

};

#endif /* __SEIDELLP_H__ */
