#pragma once

#include <array>
#include <vector>
#include <list>
#include <stack>
#include <cassert>
#include <cmath>


#include "linsolve.h"

template<typename real_t, int N>
class Vector_t
{
  private:
    std::array<real_t,N> x;
  public:
    Vector_t() {}
    Vector_t(const real_t y) 
    {
      for (int l = 0; l < N; l++)
        x[l] = y;
    }
    Vector_t(const real_t *y) 
    {
      for (int l = 0; l < N; l++)
        x[l] = y[l];
    }
    real_t& operator[](const int i)       { return x[i]; }
    real_t  operator[](const int i) const { return x[i]; }

    /******************/
    friend real_t dot(const Vector_t &a, const Vector_t &b)
    {
      real_t sum = 0;
      for (int l = 0; l < N; l++)
        sum += a[l]*b[l];
      return sum;
    }
    friend real_t norm2(const Vector_t &a)
    {
      return dot(a,a);
    }
    friend real_t norm(const Vector_t &a)
    {
      return std::sqrt(norm2(a));
    }
    /******************/
    Vector_t& operator*=(const Vector_t &a)
    {
      for (int l = 0; l < N; l++)
        x[l] *= a[l];
      return *this;
    }
    Vector_t& operator*=(const real_t a)
    {
      for (int l = 0; l < N; l++)
        x[l] *= a;
      return *this;
    }
    friend Vector_t operator*(const Vector_t &a, const Vector_t &b) 
    {
      Vector_t res;
      for (int l = 0; l < N; l++)
        res[l] = a[l]*b[l];
      return res;
    }
    friend Vector_t operator*(const Vector_t &a, const real_t b) 
    {
      Vector_t res;
      for (int l = 0; l < N; l++)
        res[l] = a[l]*b;
      return res;
    }
    friend Vector_t operator*(const real_t a, const Vector_t &b) 
    {
      return b*a;
    }
    /******************/
    Vector_t& operator+=(const Vector_t &a)
    {
      for (int l = 0; l < N; l++)
        x[l] += a[l];
      return *this;
    }
    Vector_t& operator+=(const real_t a)
    {
      for (int l = 0; l < N; l++)
        x[l] += a;
      return *this;
    }
    friend Vector_t operator+(const Vector_t &a, const Vector_t &b) 
    {
      Vector_t res;
      for (int l = 0; l < N; l++)
        res[l] = a[l]+b[l];
      return res;
    }
    friend Vector_t operator+(const Vector_t &a, const real_t b) 
    {
      Vector_t res;
      for (int l = 0; l < N; l++)
        res[l] = a[l]+b;
      return res;
    }
    friend Vector_t operator+(const real_t a, const Vector_t &b) 
    {
      return b+a;
    }
    /******************/
    Vector_t& operator-=(const Vector_t &a)
    {
      for (int l = 0; l < N; l++)
        x[l] -= a[l];
      return *this;
    }
    Vector_t& operator-=(const real_t a)
    {
      for (int l = 0; l < N; l++)
        x[l] -= a;
      return *this;
    }
    friend Vector_t operator-(const Vector_t &a, const Vector_t &b) 
    {
      Vector_t res;
      for (int l = 0; l < N; l++)
        res[l] = a[l]-b[l];
      return res;
    }
    friend Vector_t operator-(const Vector_t &a, const real_t b) 
    {
      Vector_t res;
      for (int l = 0; l < N; l++)
        res[l] = a[l]-b;
      return res;
    }
    friend Vector_t operator+(const real_t a, const Vector_t &b) 
    {
      Vector_t res;
      for (int l = 0; l < N; l++)
        res[l] = a - b[l];
      return res;
    }
};

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

      static std::pair<vec_t,real_t> makePlane(const Basis &basis)
      {
        /* move origin to the plane */
        Basis basisP;
        for (int l = 0; l < NDIM-1; l++)
          basisP[l].pos = basis[l].pos - basis[NDIM-1].pos;

        /* find a unit vector that is not parallel to the plane */
        vec_t unitVec(0.0);
        unitVec[0] = 1.0;

        int el = 0;
        vec_t planeVec = basisP[el++];
        while (dot(planeVec,unitVec) == 0)
        {
          planeVec = basisP[el++];
          assert(el < NDIM-1);
        }

        /* compute plane equation */
        vec_t n = unitVec - unitVec*(planeVec * (1.0/sqrt(norm2(planeVec))));
        n *= 1.0/sqrt(norm2(n));
        real_t p = -dot(n,basis[0]);

        return std::make_pair(n,p);
      }

      Facet makeFace(const Vertex &p, const int facetIdx) const
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
    typename Facet::list facetList;
    typename Facet::vector facetVector;

    struct FacetMD  /* facet metadata */
    {
      using stack = std::stack<FacetMD>;
      typename Facet::listIterator it;
      Vertex *pBuf;
      int pbeg, pend;
      FacetMD(typename Facet::listIterator _it, Vertex *_pBuf, int _pbeg, int _pend) :
        it(_it), pBuf(_pBuf), pbeg(_pbeg), pend(_pend) {}
    };
    typename FacetMD::stack facetStack;

    bool partition(const FacetMD &fmd, Vertex *pBuf)
    {
      const int np = fmd.pend - fmd.pbeg;
      /* no particles left to partition */
      if (np == 0)
        return false;

      /* find max distances */
      real_t distMax = 0;
      Vertex  pMax;
      for (int i = fmd.pbeg; i < fmd.pend; i++)
      {
        const Vertex    &p = fmd.pBuf[i];
        const real_t dist = fmd.it->distance(p);
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
        const auto &p = fmd.pBuf[i];
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

    using Simplex =  std::array<Vertex,NDIM+1>;


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
      const auto _x = linSolve<real_t,DIM-1>(_m,_b);

      /* recontruct tangential part of the vector */
      vec_t pt(0.0);
      for (int l = 0; l < DIM-1; l++)
        pt += _x[l]*vec_t(basis[l]);

      /* normal component of the vector */
      const vec_t pn = pos - pt;

      assert(std::abs(dot(pn,pt)) < 1.0e-10*norm2(pos));

      return norm2(pn);
    }
    
    Simplex extremeSimplex;
    template<int DIM>
    static void findExtremeSimplex(const typename Vertex::vector &pos, Simplex &simplex)
    {
      findExtremeSimplex<DIM-1>(pos, simplex);
      const int np = pos.size();
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

    void findExtremeSimplex(const typename Vertex::vector &pos)
    {
      findExtremeSimplex<NDIM+1>(pos,extremeSimplex);
    }

    void convexHull(const typename Vertex::vector &pos)
    {
      typename Vertex::vector pos1(pos), pos2(pos.size());

      findExtremeSimplex(pos);

      while (!facetStack.empty())
        facetStack.pop();

      Vertex *pBuf1 = &pos1[0];
      Vertex *pBuf2 = &pos2[0];
      while (!facetStack.empty())
      {
        const auto fmd = facetStack.top();
        facetStack.pop();
        partition(fmd, fmd.pBuf == pBuf1 ? pBuf2 : pBuf1);
      }

      facetVector.clear();
      for (auto it = facetList.begin(); it != facetList.end(); it++)
        facetVector.push_back(*it);
    }

    size_t getNumFacets() const { return facetVector.size(); }
    const Facet& getFacet(const int i) const {return facetVector[i]; }
};




