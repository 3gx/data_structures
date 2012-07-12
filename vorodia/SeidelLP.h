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
  real dist(const vec3 &p) const
  {
    return n*p - h;
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

  vec3 project(const vec3 &p) const
  {
    return p - n*(p*n - h);
  }
};


template<const int N>
struct SeidelLP
{
  private:
    int n;
    vec3 bmax;
    vec3 cvec;
    HalfSpace halfSpaceList[N];
    HalfSpace boundaryBox  [6];
    int n1, n2, n3;

  public:
    SeidelLP() : n(0) { set_bnd(vec3(1.0)); }
    SeidelLP(const HalfSpace::Vector &plist, const vec3 bmax = 1.0)
    {
      assert((int)plist.size() <= N);
      set_bnd(bmax);
      n = plist.size();
      for (int i = 0; i < n; i++)
        halfSpaceList[i] = plist[i];
    }

    void clear(const vec3 bmax = 1.0) { set_bnd(); }
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

    void set_bnd(const vec3 &bmax)
    {
      this->bmax = bmax;
      const real x = bmax.x;
      const real y = bmax.y;
      const real z = bmax.z;
      boundaryBox[0] = HalfSpace(vec3(+1.0, 0.0, 0.0), vec3(0.0,0.0,0.0));
      boundaryBox[1] = HalfSpace(vec3(-1.0, 0.0, 0.0), vec3( x ,0.0,0.0));
      boundaryBox[2] = HalfSpace(vec3( 0.0,+1.0, 0.0), vec3(0.0,0.0,0.0));
      boundaryBox[3] = HalfSpace(vec3( 0.0,-1.0, 0.0), vec3(0.0, y ,0.0));
      boundaryBox[4] = HalfSpace(vec3( 0.0, 0.0,+1.0), vec3(0.0,0.0,0.0));
      boundaryBox[5] = HalfSpace(vec3( 0.0, 0.0,-1.0), vec3(0.0,0.0, z ));
    }

    vec3 solve(const vec3 &cvec)
    {
      this->cvec = cvec;

      n1 = n2 = n3 = 0;
      //std::random_shuffle(halfSpaceList, halfSpaceList+n);
      const vec3 v = solve_lp3D(n);

      fprintf(stderr, "n3= %d\n", n3);
      fprintf(stderr, "n2= %d\n", n2);
      fprintf(stderr, "n1= %d\n", n1);

#if 0  /* sanity check */
      for (int i = 0; i < n; i++)
        assert(!halfSpaceList[i].outside(v));
#endif

      const HalfSpace &h = halfSpaceList[n-1];
      fprintf(stderr, "%g \n",  h.dist(v));
      fprintf(stderr, "%g \n",  h.dist(vec3(0.3, 0.3, 0.0)));

      return v;
    }

    template<class T1, class T2>
      struct cmp_data
      {
        bool operator()(const std::pair<T1, T2> &lhs, const std::pair<T1, T2> &rhs) const
        {
          return lhs.first < rhs.first;
        }
      };

    vec3 solve_lp3D(const int n)
    {
      const real x = bmax.x;
      const real y = bmax.y;
      const real z = bmax.z;

      vec3 vtx[8];
      std::pair<real, int> fval[8];

      vtx[0] = vec3(0.0, 0.0, 0.0);
      vtx[1] = vec3( x , 0.0, 0.0);
      vtx[2] = vec3(0.0,  y , 0.0);
      vtx[3] = vec3( x,   y , 0.0);
      vtx[4] = vec3(0.0, 0.0,  z );
      vtx[5] = vec3( x , 0.0,  z );
      vtx[6] = vec3(0.0,  y ,  z );
      vtx[7] = vec3( x,   y ,  z );

      for (int i = 0; i < 8; i++)
        fval[i] = std::make_pair(vtx[i]*cvec, i);

      std::sort(fval, fval+8, cmp_data<real, int>());
      vec3 v = vtx[fval[7].second];

      for (int i = 0; i < 8; i++)
      {
        fprintf(stderr, "i= %d : c= %g  idx= %d\n",
            i, fval[i].first, fval[i].second);
      }

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
//      fprintf(stderr, " n2= %d  v= %g %g %g", n, v.x, v.y, v.z);
      v = h1.project(v);
//      fprintf(stderr, " vp= %g %g %g \n", v.x, v.y, v.z);
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
      v = h2.project(v);
//      fprintf(stderr, "  n1= %d\n", n);
      for (int i = 0; i < 6; i++)
      {
        const HalfSpace &h = boundaryBox[i];
        if (h.outside(v))
          v = intersect(h, h1, h2);
      }

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
