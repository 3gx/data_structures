#pragma once

#include <array>
#include <vector>
#include <list>
#include <stack>
#include <cassert>
#include <cmath>
#include <iostream>

#define _EPS (1.0e-10)


#include "vector.h"
#include "linsolve.h"


template<int N>
struct QHull_t
{
  public:
    enum {NDIM = N};
    typedef double real_t;
    typedef int   id_t;

    using vec_t = Vector_t<real_t,NDIM>;

    struct Vertex
    {
      using vector = std::vector<Vertex>;
      vec_t pos;
      id_t  idx;
      real_t& operator[](const int i)       { return pos[i]; }
      real_t  operator[](const int i) const { return pos[i]; }
      operator const vec_t&() const {return pos;}
      operator vec_t&() {return pos;}
    };
    using Basis = std::array<Vertex,NDIM>;
    
    static std::pair<vec_t,real_t> planeEquation(
        const Basis &vtx,  /* first NDIM vectors for facets, NDIM+1 vector is used for orientation */
        const vec_t &posO,
        real_t  orientation)
    {
      /* compte centre of the facet */
      vec_t centre(0.0);
      for (int l = 0; l < NDIM; l++)
        centre += vtx[l];
      centre *= 1.0/NDIM;

      /* compute facet basis */
      std::array<vec_t,NDIM-1> basis;
      for (int l = 0; l < NDIM-1; l++)
      {
        basis[l] = vtx[l] - centre;
      }


      /* compute orientation vector */
      const vec_t &pos = posO - centre;
      

      std::array<std::array<real_t,NDIM-1>,NDIM-1> _m;
      std::array<real_t,NDIM-1> _b;
      for (int l = 0; l < NDIM-1; l++)
      {
        for (int ll = 0; ll < NDIM-1; ll++)
          _m[l][ll] = dot(basis[l],basis[ll]);
        _b[l] = dot(basis[l],pos);
      }

      /* solve coefficients */
      const auto& _x = linSolve<real_t,NDIM-1>(_m,_b);

      /* recontruct tangential part of the vector */
      vec_t pt(0.0);
      for (int l = 0; l < NDIM-1; l++)
        pt += _x[l]*vec_t(basis[l]);

      /* normal component of the vector */
      const vec_t &pn = pos - pt;

      const vec_t  n = pn * (1.0/norm(pn)) * orientation;
      const real_t p = -dot(n , centre);

      return std::make_pair(n,p);
    }


    struct Facet
    {
      using list = std::list<Facet>;
      using listIterator = typename list::iterator;
      using vector = std::vector<Facet>;

      std::pair<vec_t,real_t> plane;
      Basis vtx;     /* vertecies are stored in right-handed orientation */
      real_t distance(const vec_t &pos) const
      {
        return dot(plane.first,pos) + plane.second;
      }

      Facet makeFacet(const Vertex &p, const int facetIdx) const
      {
        assert(facetIdx >= 0 && facetIdx < NDIM);

        Facet f;
        f.vtx           = vtx;
        f.vtx[facetIdx] = p;
        f.plane         = planeEquation(f.vtx, vtx[facetIdx],-1.0);

        return f;
      }
    };
    typename Facet::list   facetList;
    typename Facet::vector facetVector;

    struct FacetMD  /* facet metadata */
    {
      using stack  = std::stack<FacetMD>;
      using vector = std::vector<FacetMD>;
      typename Facet::listIterator it;
      Vertex *pBuf;
      int pbeg, pend;
      FacetMD(typename Facet::listIterator _it, Vertex *_pBuf, int _pbeg, int _pend) :
        it(_it), pBuf(_pBuf), pbeg(_pbeg), pend(_pend) {}
    };
    typename FacetMD::stack facetStack;

    typename FacetMD::vector partition(const FacetMD &fmd, Vertex *pBuf)
    {
      typename FacetMD::vector newFacetsMD;
      const int np = fmd.pend - fmd.pbeg;
      /* no particles left to partition */
      if (np == 0)
        return newFacetsMD;

      /* find max distances */
      real_t distMax = 0;
      int iMax = -1;
      Vertex  pMax;
      for (int i = fmd.pbeg; i < fmd.pend; i++)
      {
        const Vertex    &p = fmd.pBuf[i];
        const real_t dist = fmd.it->distance(p);
        if (dist > distMax)
        {
          distMax = dist;
          pMax    = p;
          iMax    = i;
        }
      }
      assert(distMax > 0);

#if 0
      {
        std::cout << "splot";
        std::cout << " '-' with points notitle, '-' with lines lc 3 notitle, '-' with lines lc 2 notitle";
        // std::cout << " '-' with points notitle, '-' with lines lc 3 notitle, '-' with vector lc 2 notitle";
        std::cout << ";\n";
        fprintf(stderr, " beg= %d  end= %d\n", fmd.pbeg, fmd.pend);
        for (int i = fmd.pbeg; i < fmd.pend; i++)
        {
          const auto &p = fmd.pBuf[i];
          assert(fmd.it->distance(p) >= 0.0);
          if (i == 380)
          {
            for (int l = 0; l < NDIM; l++)
              std::cout << p[l] << " ";
            fprintf(stderr, "dist= %g \n", fmd.it->distance(p));
            fprintf(stderr, "dist0= %g \n", fmd.it->distance(pMax));
          }
          std::cout << "\n";
        }
        std::cout << "e\n";
      }
#endif

      Facet facets[NDIM];
      for (int i = 0; i < NDIM; i++)
      {
        facets[i] = fmd.it->makeFacet(pMax, i);
#if 0
        {
          const int vtx[4] =  {0,1,2,0};
          for (auto ii : vtx)
          {
            for (int l = 0; l < 3 ; l++)
              std::cout << facets[i].vtx[ii][l] << " ";
            std::cout << "\n";
          }
        }
#endif
      }
#if 0
      std::cout << "e\n";
#endif

#if 0
      {
        const int vtx[4] =  {0,1,2,0};
        for (auto i : vtx)
        {
          for (int l = 0; l < 3 ; l++)
            std::cout << fmd.it->vtx[i][l] << " ";
          std::cout << "\n";
        }
        std::cout << "e\n";
#if 0
        const auto &plane = fmd.it->plane;
        std::cout << "0";
        std::cout << " 0";
        std::cout << " 0 ";
        std::cout << plane.first[0] << " ";
        std::cout << plane.first[1] << " ";
        std::cout << plane.first[2] << " ";
        std::cout << "\ne\n";
#endif
      }
#endif

      /* count particles belonging to each of the new face */

      int count[NDIM] = {0};
      std::vector<int> whichFacet(np);

      for (int i = fmd.pbeg; i < fmd.pend; i++)
        assert(fmd.it->distance(fmd.pBuf[i]) > 0.0);

      for (int i = fmd.pbeg; i < fmd.pend; i++)
      {
        if (i == iMax)
          continue;
        const auto &p = fmd.pBuf[i];
        bool used = false;
        for (int l = 0; l < NDIM; l++)
          if (facets[l].distance(p) > 0)
          {
            if (used)
              for (int ll = 0; ll < NDIM; ll++)
                fprintf(stderr," i= %d ll= %d : dist= %g \n",
                    i, ll, facets[ll].distance(p));
            assert(!used);
            count[l]++;
            whichFacet[i-fmd.pbeg] = l;
            used = true;
            break;
          }
      }

      /* compute offset */
      int pbeg[NDIM] = {fmd.pbeg};
      int pend[NDIM] = {fmd.pbeg};
      for (int l = 1; l < NDIM; l++)
      {
        pbeg[l] = pbeg[l-1] + count[l-1];
        pend[l] = pbeg[l  ];
      }

      /* sort */
      for (int i = fmd.pbeg; i < fmd.pend; i++)
      {
        if (i == iMax)
          continue;
        const int l = whichFacet[i-fmd.pbeg];
        pBuf[pend[l]++] = fmd.pBuf[i];
      }

      /* remove old face */
      facetList.erase(fmd.it);

      newFacetsMD.reserve(NDIM);
      /* push new faces to the stack */
      for (int l = 0; l < NDIM; l++)
      {
        auto it = fmd.it;
        facetList.insert(it, facets[l]);
        it--;
        newFacetsMD.push_back(FacetMD(it, pBuf, pbeg[l], pend[l]));
      }

      return newFacetsMD;
    }

    using Simplex =  std::array<Vertex,NDIM+1>;
    Simplex extremeSimplex;

    template<int DIM>
    static real_t distance(const Simplex &simplex, vec_t pos)
    {
      assert(DIM <= NDIM);
      /* compute centre of the subspace */
      vec_t centre(0.0);
      for (int l = 0; l < DIM; l++)
        centre += simplex[l];
      centre *= 1.0/DIM;

      /* compute subspace basis */
      std::array<vec_t,DIM-1> basis;
      for (int l = 0; l < DIM; l++)
        basis[l] = simplex[l] - centre;


      /* bring pos origin to centreS */
      pos -= centre;

      /* compute basis matrix and projection of posS onto basis */
      std::array<std::array<real_t,DIM-1>, DIM-1> _m;
      std::array<real_t,DIM-1> _b;
      for (int l = 0; l < DIM-1; l++)
      {
        for (int ll = 0; ll < DIM-1; ll++)
          _m[l][ll] = dot(basis[l],basis[ll]);
        _b[l] = dot(basis[l], pos);
      }

      /* solve coefficients */
      const auto& _x = linSolve<real_t,DIM-1>(_m,_b);

      /* recontruct tangential part of the vector */
      vec_t pt(0.0);
      for (int l = 0; l < DIM-1; l++)
        pt += _x[l]*vec_t(basis[l]);

      /* normal component of the vector */
      const vec_t &pn = pos - pt;

      assert(std::abs(dot(pn,pt)) < _EPS*norm2(pos));

      return norm2(pn);
    }
    
    template<int DIM>
    static void findExtremeSimplex(const Vertex *pos, Simplex &simplex, const int np)
    {
      findExtremeSimplex<DIM-1>(pos, simplex,np);
      real_t distMax = 0;
      for (int i = 0; i < np; i++)
      {
        const Vertex &p = pos[i];
        const real_t dist = distance<DIM-1>(simplex,p);
        if (dist > distMax)
        {
          distMax = dist;
          simplex[DIM-1] = p;
        }
      }
    }


    void findExtremeSimplex(Vertex *pBuf1, Vertex *pBuf2, const int np)
    {
      Simplex simplex;
      findExtremeSimplex<NDIM+1>(pBuf1, simplex, np);

      extremeSimplex = simplex;

      /* reconstruct facets of the simplex with an outward looking normal */
      Facet facets[2];

      for (int l = 0; l < NDIM; l++)
        facets[0].vtx[l] = facets[1].vtx[l] = simplex[l];
      facets[0].plane = planeEquation(facets[0].vtx, simplex[NDIM],-1.0);
      facets[1].plane = planeEquation(facets[1].vtx, simplex[NDIM],+1.0);
      assert(std::abs(dot(facets[0].plane.first, facets[1].plane.first) + 1.0) < _EPS);

      /* sanity check */
      for (int i = 0; i < np; i++)
      {
        const auto &p = pBuf1[i];
        assert(
                (facets[0].distance(p) >= 0.0 && facets[1].distance(p) <= 0.0)
            ||  (facets[1].distance(p) >= 0.0 && facets[0].distance(p) <= 0.0)
         );
      }


      int count = 0;
      for (int i = 0; i < np; i++)
      {
        const auto &p = pBuf1[i];

        bool inUse = false;
        for (int l = 0; l < NDIM; l++)
          if (p.idx == facets[0].vtx[l].idx)
          {
            inUse = true;
            break;
          }
        if (inUse) 
          continue;
        if (facets[0].distance(p) >= 0.0)
          count++;
      }

      /*********** {outer, inner} ********/
      int pbeg[] = {0, count};
      int pend[] = {0, count};

      for (int i = 0; i < np; i++)
      {
        const auto &p = pBuf1[i];

        bool inUse = false;
        for (int l = 0; l < NDIM; l++)
          if (p.idx == facets[0].vtx[l].idx)
          {
            inUse = true;
            break;
          }
        if (inUse) 
          continue;

        if (facets[0].distance(p) >= 0.0)
          pBuf2[pend[0]++] = p;
        else
        {
          assert(facets[1].distance(p) >= 0.0);
          pBuf2[pend[1]++] = p;
        }
      }
      fprintf(stderr, " 0: pbeg= %d  pend= %d\n", pbeg[0],pend[0]);
      fprintf(stderr, " 1: pbeg= %d  pend= %d\n", pbeg[1],pend[1]);
      assert(pend[1] == np-NDIM);

      assert(facetList.empty());
      typename FacetMD::vector fmdVec;
      fmdVec.reserve(2);
      for (int l = 0; l < 2; l++)
      {
        facetList.push_back(facets[l]);
        auto it = facetList.end(); 
        it--;
        fmdVec.push_back(FacetMD(it, pBuf2, pbeg[l], pend[l]));
      }

      const auto& newFacetsMD = partition(fmdVec[1], pBuf1);
      for (auto &fmd : newFacetsMD)
        facetStack.push(fmd);
    }

    void convexHull(const typename Vertex::vector &pos)
    {
      typename Vertex::vector pos1(pos), pos2(pos.size());
      Vertex *pBuf1 = &pos1[0];
      Vertex *pBuf2 = &pos2[0];

      const int np = pos.size();
      findExtremeSimplex(pBuf1, pBuf2, np);
      return;

      while (!facetStack.empty())
        facetStack.pop();

      while (!facetStack.empty())
      {
        const auto fmd = facetStack.top();
        facetStack.pop();
        const auto& newFacetsMD = partition(fmd, fmd.pBuf == pBuf1 ? pBuf2 : pBuf1);
        for (auto &fmd : newFacetsMD)
          facetStack.push(fmd);
      }

      facetVector.clear();
      for (auto it = facetList.begin(); it != facetList.end(); it++)
        facetVector.push_back(*it);
    }

    size_t getNumFacets() const { return facetVector.size(); }
    const Facet& getFacet(const int i) const {return facetVector[i]; }
};




