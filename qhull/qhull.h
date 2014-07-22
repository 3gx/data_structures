#pragma once

#include <array>
#include <vector>
#include <list>
#include <stack>
#include <cassert>
#include <cmath>

struct QHull
{
  enum {NDIM = 3};
  typedef float real_t;

  typedef std::array<real_t,NDIM > vec;

  struct pos_t
  {
    vec pos;
    int idx;
    real_t& operator[](const int i)       { return pos[i]; }
    real_t  operator[](const int i) const { return pos[i]; }
  };

  typedef std::vector<pos_t> PosVector;

  struct Facet;
  typedef std::list<Facet> FacetList;

  struct Facet
  {
    std::array<real_t,NDIM+1> plane;  /* normalized plane equation */
    std::array< pos_t,NDIM  > vtx;    /* right-handed orientation */
    real_t distance(const vec &pos) const
    {
      real_t dist = plane[NDIM];
      for (int i = 0; i < NDIM; i++)
        dist += pos[i]*plane[i];
      return dist;
    }

    std::array<real_t,NDIM+1> makePlane(const std::array<pos_t,NDIM> &vtx) const
    {
      std::array<real_t,NDIM+1> plane;

      plane[0] = vtx[0][0];
      plane[1] = vtx[1][0];
      plane[2] = vtx[2][0];
      plane[3] = 0;
      //std::array<real_t,NDIM> n = cross(vtx[0], vtx[1]);

      return plane;
    }

    Facet newFace(const pos_t &p, const int facetIdx) const
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


  FacetList facetList;
  PosVector _pos1, _pos2;


  struct FacetMD
  {
    FacetList::iterator it;
    pos_t *pBuf;
    int pbeg, pend;
    FacetMD(FacetList::iterator _it, pos_t *_pBuf, int _pbeg, int _pend) :
      it(_it), pBuf(_pBuf), pbeg(_pbeg), pend(_pend) {}
  };
  std::stack<FacetMD> faceStack;

  bool partition(const FacetMD &fmd, pos_t *pBuf)
  {
    /* no particles left to partition */
    if (fmd.pend - fmd.pbeg == 0)
      return false;

    /* find max distances */
    real_t distMax =  0;
    int    idxMax  = -1;
    for (int i = fmd.pbeg; i < fmd.pend; i++)
    {
      const real_t dist = fmd.it->distance(fmd.pBuf[i].pos);
      if (dist > distMax)
      {
        distMax = dist;
        idxMax  = i;
      }
    }
    assert(distMax > 0);

    Facet facets[NDIM];
    const pos_t pMax = fmd.pBuf[idxMax];
    for (int i = 0; i < NDIM; i++)
      facets[i] = fmd.it->newFace(pMax, i);

    /* partition particle ditribution to new facets */

    int count[NDIM] = {0};
    std::vector<int> whichFacet(fmd.pend-fmd.pbeg);

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
      faceStack.push(FacetMD(it, pBuf, pbeg[l], pend[l]));
    }

    return true;
  }

  void extremeSimplex(const PosVector &pos)
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

  void computeConvexHull(const PosVector &pos)
  {
    std::vector<pos_t> pos1(pos), pos2(pos.size());

    extremeSimplex(pos);

    while (!faceStack.empty())
      faceStack.pop();

    pos_t *pBuf1 = &pos1[0];
    pos_t *pBuf2 = &pos2[0];
    while (!faceStack.empty())
    {
      const auto fmd = faceStack.top();
      faceStack.pop();
      partition(fmd, fmd.pBuf == pBuf1 ? pBuf2 : pBuf1);
    }
  }


};




