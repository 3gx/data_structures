#pragma once

#include <array>
#include <vector>
#include <list>
#include <stack>
#include <cassert>
#include <cmath>

template<typename real_t, int N>
class Vector_t
{
  private:
    std::array<real_t,N> x;
  public:
    real_t& operator[](const int i)       { return x[i]; }
    real_t  operator[](const int i) const { return x[i]; }

    friend real_t dot(const Vector_t &a, const Vector_t &b)
    {
      real_t sum = 0;
      for (int l = 0; l < N; l++)
        sum += a[l]*b[l];
      return sum;
    }
};


struct QHull
{
  enum {NDIM = 3};
  typedef float real_t;
  typedef int   id_t;

  typedef Vector_t <real_t,NDIM> vec_t;

  struct pos_t
  {
    typedef std::vector<pos_t> vector;
    vec_t pos;
    id_t  idx;
    real_t& operator[](const int i)       { return pos[i]; }
    real_t  operator[](const int i) const { return pos[i]; }
  };
  typedef std::array<pos_t,NDIM> vtx_t;

  struct Facet
  {
    typedef std::list<Facet> list;
    typedef list::iterator iterator;

    std::pair<vec_t,real_t> plane;
    vtx_t vtx;     /* vertecies are stored in right-handed orientation */
    real_t distance(const vec_t &pos) const
    {
      return dot(plane.first,pos) + plane.second;
    }

    static std::pair<vec_t,real_t> makePlane(const vtx_t &vtx)
    {
      vec_t n;
      real_t p;

      n[0] = vtx[0][0];
      n[1] = vtx[1][0];
      n[2] = vtx[2][0];
      p = 0;
      //std::array<real_t,NDIM> n = cross(vtx[0], vtx[1]);

      return std::make_pair(n,p);
    }

    Facet makeFace(const pos_t &p, const int facetIdx) const
    {
      assert(facetIdx >= 0 && facetIdx < NDIM);

      Facet f;
      f.vtx[facetIdx] = p;
      for (int i = 0; i < NDIM; i++)
      {
        if (facetIdx != i)
          f.vtx[i] = vtx[i];
      }

      f.plane = makePlane(f.vtx);

      return f;
    }
  };
  Facet::list facetList;

  struct FacetMD  /* facet metadata */
  {
    typedef std::stack<FacetMD> stack;
    Facet::iterator it;
    pos_t *pBuf;
    int pbeg, pend;
    FacetMD(Facet::iterator _it, pos_t *_pBuf, int _pbeg, int _pend) :
      it(_it), pBuf(_pBuf), pbeg(_pbeg), pend(_pend) {}
  };
  FacetMD::stack facetStack;

  bool partition(const FacetMD &fmd, pos_t *pBuf)
  {
    const int np = fmd.pend - fmd.pbeg;
    /* no particles left to partition */
    if (np == 0)
      return false;

    /* find max distances */
    real_t distMax = 0;
    pos_t  pMax;
    for (int i = fmd.pbeg; i < fmd.pend; i++)
    {
      const pos_t    &p = fmd.pBuf[i];
      const real_t dist = fmd.it->distance(p.pos);
      if (dist > distMax)
      {
        distMax = dist;
        pMax    = p;
      }
    }
    assert(distMax > 0);

    Facet facets[NDIM];
    for (int i = 0; i < NDIM; i++)
      facets[i] = fmd.it->makeFace(pMax, i);

    /* count particles belonging to each of the new face */

    int count[NDIM] = {0};
    std::vector<int> whichFacet(np);

    for (int i = fmd.pbeg; i < fmd.pend; i++)
    {
      const auto &p = fmd.pBuf[i].pos;
      bool used = false;
      for (int l = 0; l < NDIM; l++)
        if (facets[l].distance(p) >= 0)
        {
          assert(!used);
          count[l]++;
          whichFacet[i-fmd.pbeg] = l;
          used = true;
        }
    }

    /* compute offset */
    int pbeg[NDIM] = {0};
    int pend[NDIM] = {0};
    for (int l = 1; l < NDIM; l++)
    {
      pbeg[l] = pbeg[l-1] + count[l-1] + fmd.pbeg;
      pend[l] = pbeg[l  ];
    }
    assert(pbeg[NDIM-1] + count[NDIM-1] == fmd.pend - fmd.pbeg);

    /* sort */
    for (int i = fmd.pbeg; i < fmd.pend; i++)
    {
      const int l = whichFacet[i-fmd.pbeg];
      pBuf[pend[l]++] = fmd.pBuf[i];
    }

    /* remove old face */
    facetList.erase(fmd.it);

    /* push new faces to the stack */
    for (int l = 0; l < NDIM; l++)
    {
      auto it = fmd.it;
      facetList.insert(it, facets[l]);
      it--;
      facetStack.push(FacetMD(it, pBuf, pbeg[l], pend[l]));
    }

    return true;
  }

  void extremeSimplex(const pos_t::vector &pos)
  {
    facetList.clear();
    std::array<pos_t,NDIM+1> vtxList;

    real_t xMin = +HUGE, xMax = -HUGE;
    const int np = pos.size();
    for (int i = 0; i < np; i++)
    {
      const auto &p = pos[i];
      if (p[0] < xMin)
      {
        xMin       = p[0];
        vtxList[0] = p;
      }
      if (p[0] > xMax)
      {
        xMax       = p[0];
        vtxList[1] = p;
      }
    }
    assert(vtxList[0].idx != vtxList[1].idx);

    for (int l = 2; l < NDIM+1; l++)
    {
      real_t dMax = -HUGE;
      for (int i = 0; i < np; i++)
      {
//        bool use = true;
        real_t d = 0, n = 0;
        for (int ll = 0; ll < l; ll++)
        {
//          d += 

        }
      }


    }

  }

  void computeConvexHull(const pos_t::vector &pos)
  {
    pos_t::vector pos1(pos), pos2(pos.size());

    extremeSimplex(pos);

    while (!facetStack.empty())
      facetStack.pop();

    pos_t *pBuf1 = &pos1[0];
    pos_t *pBuf2 = &pos2[0];
    while (!facetStack.empty())
    {
      const auto fmd = facetStack.top();
      facetStack.pop();
      partition(fmd, fmd.pBuf == pBuf1 ? pBuf2 : pBuf1);
    }
  }


};




