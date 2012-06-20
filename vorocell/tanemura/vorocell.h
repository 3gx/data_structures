#ifndef __VOROCELL_H__
#define __VOROCELL_H__

/* Tanemura's algorithm */

#include <cassert>
#include <vector>
#include <deque>
#include "vector3.h"

namespace Voronoi
{
  typedef double real;
  typedef vector3<real>  vec3;

  struct Site
  {
    typedef std::vector<Site> Vector;
    long long idx;
    vec3 pos;
  };

  struct Plane
  {
    vec3 n;
    Plane(const vec3 &ipos, const vec3 &jpos)
    {
       n = ipos%jpos;
    }
    real operator()(const vec3 pos) const
    {
      return n*pos;
    }
  };


  template<class T, const int N>
    class Array
    {
      private:
        int n;
        T data[N];

      public:
        Array(const int _n = 0) : n(_n) {}
        Array(const std::vector<T> &_data)
        {
          n = _data.size();
          assert(n <= N);
          for (int i = 0; i < n; i++)
            data[i] = _data[i];
        }
        void clear()
        {
          n = 0;
        }
        const T& operator[](const int i) const {return data[i];}
        T& operator[](const int i)       {return data[i];}

        void push_back(const T &t) {assert(n<N); data[n++] = t;}
        int size() const { return n; }
        void resize(const int size) 
        {
          assert (size <= N);
          n = size;
        }

        int capacity() const {return N;}

        bool erase(const T &t) 
        {
          for (int i = 0; i < n; i++)
            if (data[i] == t)
            {
              n--;
              std::swap(data[i], data[n]);
              return true;
            }
          return false;
        }

        T erase(const int i)
        {
          assert(i >= 0);
          assert(i < n);
          n--;
          std::swap(data[i], data[n]);
          return data[n];
        }
    };

  template<const int N, const int NCRIT>
    struct ConnectivityMatrix2D
    {
      private:
        int matrix[N*(N+1)/2];
        int _flags[N];

      public:

        ConnectivityMatrix2D() 
        {
          clear();
        }

        void clear()
        {
          for (int i = 0; i < N*(N+1)/2; i++)
            matrix[i] = 0;

          for (int i = 0; i < N; i++)
            _flags[i] = 0;
        }

        /* i: x
         * j: y */
        unsigned int  operator()(const int i, const int j) const { return matrix[map(i,j)]; }
        unsigned int  inc       (const int i, const int j) 
        {
          unsigned int &val = matrix[map(i,j)];
          val++;
          const int inc = (int)(val == NCRIT) - (int)(val == NCRIT+1);
          _flags[i] += inc;
          _flags[j] += inc;

          return val;
        }

        /************/

        bool isFullyConnected(const int i) const {return _flags[i] == 0;}

      private:

        /* i: x
         * j: y */
        int map(const int i, const int j) const
        {
#if 1 /* LOWER_PACKED */
          return i + ((2*N-j)*(j-1)>>1);
#else /* UPPER_PACKED */
          return i + (j*(j-1)>>1);
#endif
        }
    };

  struct Tetrahedron
  {
    private:
      int v1, v2, v3;
    public:
      Tetrahedron(const int _v1, const int _v2, const int _v3) : v1(_v1), v2(_v2), v3(_v3) {}

      int vertex1() const {return v1;}
      int vertex2() const {return v2;}
      int vertex3() const {return v3;}

      std::pair<int,int> pair(const int vX) const
      {
        if      (v1 == vX) return std::make_pair(v2,v3);
        else if (v2 == vX) return std::make_pair(v1,v3);
        else if (v3 == vX) return std::make_pair(v1,v2);
        else
          assert(0);
      }
  };


  template<const int N>
    struct Cell
    {
      typedef ConnectivityMatrix2D<N, 1>      Edges;
      typedef ConnectivityMatrix2D<N, 2>  Triangles;

      private:
      Edges     edges;
      Triangles triangles;

      std::deque < int>   vertexQueue    ;
      std::vector<bool> isVertexQueued   ;
      std::vector<bool>   vertexCompleted;

      std::vector<Tetrahedron> tetrahedra;
      std::vector<   int     >     nbList;

      Array<int, N> faceVtx[N];

      public:
      Cell(const Site::Vector &sites)
      {
        build(sites);
        
      }

      void build(const Site::Vector &sites)
      {
        assert((int)sites.size() <= N);
        clear(sites.size());
      }

      private:

      void clear(const int nSites)
      {
        edges      .clear();
        triangles  .clear();

        vertexCompleted.clear();
        isVertexQueued .clear();
        vertexQueue    .clear();

        vertexCompleted.resize(nSites);
        for (int i = 0; i < nSites; i++)
          vertexCompleted[i] = false;


        tetrahedra.clear();
        nbList    .clear();

        for (int i = 0; i < N; i++)
          faceVtx[i].clear();
      }

      private:

      void buildCell()
      {
      }

      void completeCell(const Site::Vector &siteList)
      {
        const int nSites = siteList.size();
        while (!vertexQueue.empty())
        {
          /* step 4.2:
           *  extract vertex from the list
           */
          const int iVertex = vertexQueue.front();
          vertexQueue.pop_front();

#if 1  /* sanity check: the vertex haven't yet registered in nbList */
          for (std::vector<int>::const_iterator it = nbList.begin(); it != nbList.end(); it++)
            assert(*it !=  iVertex);
#endif
          nbList.push_back(iVertex);

          if (edges.isFullyConnected(iVertex))
            vertexCompleted[iVertex] = 1;

          /* step 4.3:
           *  the vertex is completed, proceed to the next one
           */
          if (vertexCompleted[iVertex]) 
          {
            assert(edges.isFullyConnected(iVertex));
            continue;
          }

          /* step 4.4:
           *  find a tetrahedron (i, iVertex, jVertex, kVertex) with at least one
           *  incomplete face
           */

          const Tetrahedron &t = tetrahedra[faceVtx[iVertex].back()];
          const std::pair<int,int> vpair = t.pair(iVertex);
          int jVertex = vpair.first;
          int kVertex = vpair.second;
          if (edges(iVertex, jVertex) > 1)
            std::swap(jVertex, kVertex);
          assert(edges(iVertex, jVertex) == 1);
          assert(edges(iVertex, kVertex)  > 1);

          /* step 4.5-4.7 */

          while(triangles(iVertex, jVertex) < 2)
          {
            /* step 4.5 - 4.6: 
             *  search a vertex on the opposite side of the kVertex 
             *  (in the half-space bounded by ijFace that does not contain kVertex)
             */
            const vec3 &ipos = siteList[iVertex];
            const vec3 &jpos = siteList[jVertex];
            const vec3 &kpos = siteList[kVertex];
            const Plane plane(ipos, jpos);
            const int  sideK = plane(kpos) > 0.0;
            real largeNum = +1e10;
            vec3  cpos(0.0);
            int lVertex = -1;

            /* hot-spot: finding 4th vertex of the new tetrahedron */

            for (int i = 0; i < nSites; i++)
            {
              const vec3 &pos = siteList[i];
              const int  side = plane(pos) > 0.0;
              const real dist = pos*(pos + cpos) + largeNum;
              const bool chck = !vertexCompleted[i] &&
                (i != iVertex) && (i != jVertex) && (i != kVertex);

              if (dist < 0.0 && side^sideK = 1 && chck)
              {
                real radius = 0.0;
                cpos = sphere(ipos, jpos, pos, radius);
                largeNum = 0.0;
                lVertex = i;
              }
            }

            assert(lVertex >= 0);

            /* step 4.7:
             *  register new tetrahedron (iVertex, jVertex, lVertex) 
             */
            tetrahedra.push_back(Tetrahedron(iVertex, jVertex, lVertex));
            faceVtx[iVertex].push_back(tetrahedra.size()-1);
            faceVtx[jVertex].push_back(tetrahedra.size()-1);
            faceVtx[lVertex].push_back(tetrahedra.size()-1);

            if (!isVertexQueued[lVertex])
            {
              vertexQueue.push_back(lVertex);
              isVertexQueued[lVertex] = 1;
            }

            edges.inc(iVertex, jVertex);
            edges.inc(iVertex, lVertex);
            edges.inc(jVertex, lVertex);

            triangles.inc(iVertex, jVertex);
            triangles.inc(iVertex, lVertex);
            triangles.inc(jVertex, lVertex);
            kVertex = jVertex;
            jVertex = lVertex;
          }

          /* step 4.8: */
          assert(edges.isFullyConnected(iVertex));
          vertexCompleted[iVertex] = 1;
        }
      }


      vec3 sphere(const vec3 &ip, const vec3 &jp, const vec3 &kp, real &radius)
      {
        const real di = ip.norm2();
        const real dj = jp.norm2();
        const real dk = kp.norm2();

        real XI, YI, ZI, DI;
        real XJ, YJ, ZJ, DJ;
        real XK, YK, ZK, DK;
        real XC, YC, ZC, RR, D;

        D=XI*YJ*ZK+XJ*YK*ZI+XK*YI*ZJ-XK*YJ*ZI-XI*YK*ZJ-XJ*YI*ZK;
        assert (D != 0.0);
        D=0.50/D;
        XC=DI*YJ*ZK+DJ*YK*ZI+DK*YI*ZJ-DK*YJ*ZI-DI*YK*ZJ-DJ*YI*ZK;
        XC=XC*D;
        RR=XC*XC;
        YC=XI*DJ*ZK+XJ*DK*ZI+XK*DI*ZJ-XK*DJ*ZI-XI*DK*ZJ-XJ*DI*ZK;
        YC=YC*D;
        RR=RR+YC*YC;
        ZC=XI*YJ*DK+XJ*YK*DI+XK*YI*DJ-XK*YJ*DI-XI*YK*DJ-XJ*YI*DK;
        ZC=ZC*D;
        RR=RR+ZC*ZC;

        radius = RR;
        return vec3(XC, YC, ZC);
      }



    };
}

#endif /* __VOROCELL_H__ */
