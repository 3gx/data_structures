#pragma once

#include <array>

struct QHull
{
  enum {NDIM = 3;}
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
  typedef std::List<Facet> FacetList;

  struct Facet
  {
    std::array<real_t, NDIM+1> plane;  /* normalized plane equation */
    std::array<FacetList::iterator,NDIM> ngb; /* neighbours */
    std::array<pos_t,NDIM> vtx;  /* right-handed orientation */
    real_t distance(const vec &pos) const
    {
      real_t dist = plane[NDIM];
      for (int i = 0; i < NDIM; i++)
        dist += pos[i]*plane[i];
      return dist;
    }
    Facet newFace(const pos_t &p, const int facetIdx) const
    {
      assert(idx >= 0 && idx < NDIM);

      Facet f;
      f.vtx[facetIdx] = p;
      for (int i = 0; i < NDIM; i++)
      {
        if (facetIdx != i)
          f.vtx[i] = vtx[i];
        f.ngb[i] = NULL;
      }

      f.n = makePlane(f.vtx.pos);

      return f;
    }

    std::array<real_t,NDIM+1>& makePlane(const std::array<vec,NDIM> &vtx)
    {
      std::array<real_t,NDIM+1> plane;
      std::array<real_t,NDIM> n = cross(vtx[0], vtx[1]);

      return plane;
    }

  };


  FacetList facetList;
  PosVector _pos1,_pos2;
  pos_t* pBuf, pDst;


  struct FaceMD
  {
    FaceList::iterator it;
    pos_t *pBuf;
    int pbeg, pend;
    FaceMD(FaceList::iterator _it, vec *_pBuf, int _pbeg, int _pend) :
      it(_it), pBuf(_pBuf), pbeg(_pbeg), pend(_pend) {}
  };
  stl::stack<FaceMD> faceStack;

  bool partition(const FaceMD &fmd, vec *pBuf)
  {
    /* no particles left to partition */
    if (fmd.pend - fmd.pbeg == 0)
      return false;

    /* find max distances */
    real_t distMax =  0;
    int    idxMax  = -1;
    for (int i = fmd.pbeg; i < fmd.pend; i++)
    {
      const dist = f->distance(fmd.pBuf[i].pos);
      if (dist > distMax)
      {
        distMax = dist;
        idxMax  = i;
      }
    }
    assert(distMax > 0);

    FacetList facets[NDIM];
    const pos_t pMax = fmd.pBuf[idxMax];
    for (int i = 0; i < NDIM; i++)
    {
      facets[i] = fmd.it->newFace(pMax, i);
      Facet &f = facets[i];
    
      /* fix neightbours */
    }

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
    assert(pbeg[l] + count[l] == fmd.pend - fmd.pbeg);

    /* sort */
    for (int i = fmd.pbeg; i < fmd.pend; i++)
    {
      const int l = whichFacet[i-fmd.pbeg];
      pBuf[pend[l]++] = fmd.pBuf[i];
    }

    /* remove old face */
    auto it = fmd.it;
    faceList.remove(it);

    /* insert new facets */
    facetList.insert(it, facets.begin(), facets.end());

    /* change iterator to point to the first new face */
    it -= NDIM;

    /* push new faces to the stack */
    for (int l = 0; l < NDIM; l++)
      faceStack.push(FaceMd(it++, pBuf, pbeg[l], pend[l]));

    return true;
  }

  void computeConvexHull(const PosVector &pos)
  {
    while (faceStack.empty())
      faceStack.pop();
    faceList.clear();

    _pos1.resize(pos.size());
    _pos2.resize(pos.size());

    for (auto i = 0; i < pos.size(); i++)
    {
      _pos1[i].pos = pos[i];
      _pos1[i].idx = i;
    }

    pBuf = &_pos1[0];
    pDst = &_pos2[0];
    while (faceStack.empty())
    {
      const auto f = faceStack.top();
      faceStack.pop();
      partition(f, pDst);
      std::swap(pBuf, pDst);
    }
  }


};




