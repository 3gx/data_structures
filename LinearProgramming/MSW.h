#pragma once

#include "vector3.h"
#include <vector>
#include <cassert>
#include <algorithm>
#include <cfloat>

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
    nflops += 6;
#if 0
    return n*p < h;
#else
    return n*p - h < -1.0e-10*std::abs(h);
#endif
  }
  friend std::pair<vec3, real> intersect(const HalfSpace &p1, const HalfSpace &p2, const HalfSpace &p3, const vec3 &c)
  {
    nflops += 9*3 + 5+1+3*3+3*3+3*2 + 5;
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

struct Basis
{
  int x, y, z;
  Basis() {}
  Basis(const int _x, const int _y, const int _z) : x(_x), y(_y), z(_z) {}
};


struct MSW
{
  enum {N = 1000};
  private:
    int n;
    vec3 bmax, cvec;
    HalfSpace halfSpaceList[N];
    int       halfSpaceFlag[N];

  public:
    MSW(const vec3 &_bmax = 1.0) : n(3), bmax(_bmax) 
    {
      for (int i = 0; i < N; i++)
        halfSpaceFlag[i] = 1;
    }
    MSW(const HalfSpace::Vector &plist, const vec3 &_bmax = 1.0) : bmax(_bmax)
    {
      assert((int)plist.size() <= N);
      n = plist.size() + 3;
      for (int i = 3; i < n; i++)
        halfSpaceList[i] = plist[i];
      for (int i = 0; i < N; i++)
        halfSpaceFlag[i] = 1;
    }

    void clear(const vec3 &_bmax = 1.0) { n = 3; bmax = _bmax; }
    bool push(const HalfSpace &p)
    {
      assert(n > 2);
      assert(n < N);
      halfSpaceList[n++] = p;
      return true;
    }
    int nspace() const { return n-3; }

    vec3 solve(const vec3 &cvec, const bool randomize = false)
    {
      this->cvec = cvec;
      const vec3 bvec(
          cvec.x > 0.0 ? bmax.x : -bmax.x,
          cvec.y > 0.0 ? bmax.y : -bmax.y,
          cvec.z > 0.0 ? bmax.z : -bmax.z);
      halfSpaceList[0] = HalfSpace(vec3(-bvec.x,0.0,0.0), vec3(bvec.x,0.0,0.0));
      halfSpaceList[1] = HalfSpace(vec3(0.0,-bvec.y,0.0), vec3(0.0,bvec.y,0.0));
      halfSpaceList[2] = HalfSpace(vec3(0.0,0.0,-bvec.z), vec3(0.0,0.0,bvec.z));
      halfSpaceFlag[0] = 0;
      halfSpaceFlag[1] = 0;
      halfSpaceFlag[2] = 0;

      if (randomize)
        std::random_shuffle(halfSpaceList+3, halfSpaceList+n);
      
      const vec3 v  = solve_lp3D(n);

      int i = 2;
      int j = 0;
      while (++i < n)
        if (halfSpaceFlag[i] == 0)
        {
          std::swap(halfSpaceList[i  ], halfSpaceList[j  ]);
          std::swap(halfSpaceFlag[i--], halfSpaceFlag[j++]);
        }

      return v;
    }

    vec3 solve_lp3D(const int n) __attribute__((always_inline))
    {
      vec3 v = intersect(halfSpaceList[0], halfSpaceList[1], halfSpaceList[2], cvec).first;

#if 1
      int i = 2;
      while (++i < n)
        if (halfSpaceList[i].outside(v) && halfSpaceFlag[i])
        {
          v = newBasis(i);
          i = 2;
        }
#else
      bool flag = true;
      while(flag)
      {
        flag = false;
        for (int i = 3; i < n; i++)
          if (halfSpaceList[i].outside(v) && halfSpaceFlag[i])
          {
            v = newBasis(i);
            flag = true;
            break;
          }
      }
#endif

      return v;
    }


    vec3 newBasis(const int i) __attribute__((always_inline))
    {
      HalfSpace *hs = halfSpaceList;

      std::pair<vec3, real> v[3] = {
        intersect(hs[i], hs[1], hs[2], cvec),
        intersect(hs[0], hs[i], hs[2], cvec),
        intersect(hs[0], hs[1], hs[i], cvec) };
      if (hs[0].outside(v[0].first)) v[0].second = -HUGE;
      if (hs[1].outside(v[1].first)) v[1].second = -HUGE;
      if (hs[2].outside(v[2].first)) v[2].second = -HUGE;

      int j = 0;
      if (v[1].second > v[j].second) j = 1;
      if (v[2].second > v[j].second) j = 2;

      assert(v[j].second > -HUGE);

#if 0
      {
        for (int j = 0; j < 3; j++)
          assert(!hs[i].outside(v[j].first));
      }

      {
        assert(!hs[i].outside(v[j].first));
        for (int k = 0; k < 3; k++)
          assert(!hs[k].outside(v[j].first));
      }
#endif

      std::swap(halfSpaceFlag[i], halfSpaceFlag[j]);
      std::swap(halfSpaceList[i], halfSpaceList[j]);
      return v[j].first;
    }

};


