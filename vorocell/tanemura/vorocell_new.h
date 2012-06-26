#ifndef __VOROCELL_H__
#define __VOROCELL_H__

/* Tanemura's algorithm */

#include <cassert>
#include <vector>
#include <stack>
#include <deque>
#include <set>
#include <algorithm>
#include <cmath>
#include "vector3.h"

template<class T>
inline T __min(const T a, const T b) {return a < b ? a : b;}

template<class T>
inline T __max(const T a, const T b) {return a > b ? a : b;}

template<class T>
inline T __abs(const T a) {return a < T(0.0) ? -a : a;}

template<class T>
inline T __sign(const T a) {return a < T(0.0) ? (T)-1.0 : (T)+1.0;}

struct PackedInt2
{
  private:
  unsigned int value;

  public:
  PackedInt2() {}
  PackedInt2(const unsigned int x)                       : value((x<<16) + x) {}
  PackedInt2(const unsigned int x, const unsigned int y) : value((x<<16) + y) {}

  unsigned int first () const {return (value >> 16) & 0x0000FFFFU;}
  unsigned int second() const {return (value      ) & 0x0000FFFFU;}
};

struct SortedInt2
{
  private:
  unsigned int value;

  public:
  SortedInt2() : value(0) {}
  SortedInt2(const unsigned int x, const unsigned int y) :
    value(x < y ? (y<<16) + x : (x<<16) + y)
    {
      assert(lo() < hi());
    }

  unsigned int hi() const {return (value >> 16) & 0x0000FFFFU;}
  unsigned int lo() const {return (value      ) & 0x0000FFFFU;}
  friend bool operator==(const SortedInt2 &lhs, const SortedInt2 &rhs) {return lhs.value == rhs.value;}
  friend bool operator< (const SortedInt2 &lhs, const SortedInt2 &rhs) {return lhs.value <  rhs.value;}
  friend bool operator!=(const SortedInt2 &lhs, const SortedInt2 &rhs) {return lhs.value != rhs.value;}
  friend bool operator> (const SortedInt2 &lhs, const SortedInt2 &rhs) {return lhs.value >  rhs.value;}
  friend bool operator>=(const SortedInt2 &lhs, const SortedInt2 &rhs) {return lhs.value >= rhs.value;}
  friend bool operator<=(const SortedInt2 &lhs, const SortedInt2 &rhs) {return lhs.value <= rhs.value;}
};

template<class T1, class T2>
struct cmp_data
{
  bool operator()(const std::pair<T1, T2> &lhs, const std::pair<T1, T2> &rhs) const
  {
    return lhs.first < rhs.first;
  }
};


namespace Voronoi
{
  typedef double real;
  typedef vector3<real>  vec3;

  struct Site
  {
    typedef std::vector<Site> Vector;
    vec3 pos;
    int idx;
    float r;
    Site() {}
    Site(const vec3 &_pos, const long long _idx) : pos(_pos), idx(_idx) {}
    bool operator()(const Site &s1, const Site &s2) const
    {
      return s1.r < s2.r;
    }
  };

  struct Plane
  {
    vec3 n;
    Plane(const vec3 &ipos, const vec3 &jpos, const bool normalize = false)
    {
      n = ipos%jpos;
      if (normalize)
        n *= 1.0/n.abs();
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
        const T& back() const {assert(n > 0); return data[n-1];};
        T& back() {assert(n > 0); return data[n-1];};

        void push_back(const T &t) {assert(n<N); data[n++] = t;}
        int size() const { return n; }
        bool empty() const { return n == 0;}
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

  typedef std::vector<PackedInt2> List;
  template<const int N>
    struct ConnectivityMatrix2D
    {
      private:
        List matrix[N*(N+1)/2];
        int  first[N];
        std::stack<int> list;

      public:

        ConnectivityMatrix2D() 
        {
          for (int i = 0; i < N*(N+1)/2; i++)
          {
            matrix[i].reserve(8);
          }
          for (int i = 0; i < N; i++)
            first[i] = -1;
        }

        void clear()
        {
          for (int i = 0; i < N; i++)
            first[i] = -1;
          while (!list.empty())
          {
            matrix[list.top()].clear();
            list.pop();
          }
        }

        /* i: x
         * j: y */
        const List& operator()(const int i, const int j) const { return matrix[map(i,j)]; }
        void add(const int i, const int j, const PackedInt2 p)
        {
          const int addr = map(i,j);
          if (matrix[addr].empty())
            list.push(addr);
          matrix[addr].push_back(p);
          first[i] = j;
          first[j] = i;
        }
        const int operator[](const int i) const { return first[i]; }

      private:
        int map(const int _i, const int _j) const
        {
#if 0   /* LOWER_PACKED */
          const int i = __max(_i, _j);
          const int j = __min(_i, _j) + 1;
          return i + ((2*N-j)*(j-1)>>1);
#else   /* UPPER_PACKED */
          const int i = __min(_i, _j);
          const int j = __max(_i, _j);
          return i + ((j*(j+1))>>1);
#endif
        }
    };

  struct Tetrahedron
  {
    private:
      int v1, v2, v3;
      vec3 c;
    public:
      Tetrahedron(const int _v1, const int _v2, const int _v3, const vec3 &_c) : 
        v1(_v1), v2(_v2), v3(_v3), c(_c) {}

      int vertex1() const {return v1;}
      int vertex2() const {return v2;}
      int vertex3() const {return v3;}
      const vec3& centre() const { return c; }

      std::pair<int,int> pair(const int vX) const
      {
        if      (v1 == vX) return std::make_pair(v2,v3);
        else if (v2 == vX) return std::make_pair(v1,v3);
        else if (v3 == vX) return std::make_pair(v1,v2);
        else
        {
          assert(0);
          return std::make_pair(-1,-1);
        }
      }

      bool operator==(const Tetrahedron &t) const
      {
        int list1[3] = {  v1,   v2,   v3};
        int list2[3] = {t.v1, t.v2, t.v3};
        std::sort(list1, list1+3);
        std::sort(list2, list2+3);
        return list1[0]==list2[0] && list1[1]==list2[1] && list1[2]==list2[2];
      };
      bool operator!=(const Tetrahedron &t) const { return !(*this == t); }
  };

  struct Face
  {
    vec3 norm;
    real area;
    Face() {}
    Face(const vec3 &n, const real A) : norm(n), area(A) {}
  };

  template<const int N>
    struct Cell
    {

      private:
      typedef ConnectivityMatrix2D<N>  Triangles;
//      typedef Array<int, N/2> FaceArray;
      Triangles triangles;
      real eps, eps2;

      std::deque < int>   vertexQueue    ;
      std::vector<bool> isVertexQueued   ;
      std::vector<bool>   vertexCompleted;

      std::vector<Tetrahedron> tetraList;
      std::vector<   int     >     nbList;
      std::vector<  Face     >   faceList;

      std::vector< int> vtxUse;
      std::stack < int> vtxList;
      std::vector<vec3> vtxPos;
      Site::Vector     sites;

      std::vector< std::pair<real, vec3> > angle_vec_pair;
#if 0
      std::deque<int> incompleteTetra;
      FaceArray faceVtx[N];
#endif

      real cellVolume;

      public:
      Cell(const real _eps = 1.0e-13) : eps(_eps), eps2(_eps*_eps)
      {
        vertexCompleted.reserve(N);
        isVertexQueued .reserve(N);

        tetraList      .reserve(N);
        nbList         .reserve(N);
        faceList       .reserve(N);
        vtxUse         .reserve(N);
        vtxPos         .reserve(N);
        sites          .reserve(N);
        angle_vec_pair .reserve(N);
      }
      
      private:
      void clear(const int nSites)
      {
        vertexQueue.clear();
        triangles  .clear();

        vertexCompleted.resize(nSites);
        isVertexQueued .resize(nSites);
        vtxUse         .resize(nSites);
        sites          .resize(nSites);
        for (int i = 0; i < nSites; i++)
          isVertexQueued[i] = vertexCompleted[i] = false;

        tetraList .clear();
        nbList    .clear();
        faceList  .clear();

        cellVolume = 0.0;
      }
      
      public:

      int nb() const {return nbList.size();}
      real volume() const {return cellVolume;}

      bool build(const Site::Vector &siteList)
      {
        const double t00 = get_wtime();

        const int nSite = siteList.size();
        assert(nSite <= N);
        clear(nSite);

        /* step 1: 
         *  find site i, nearest to the origin
         */
        double t1 = get_wtime();
        int i = -1;
        real r2min = HUGE;
        for (int ix = 0; ix < nSite; ix++)
        {
          const vec3 &pos = siteList[ix].pos;
          const real   r2 = pos.norm2();
          assert(r2 > 0.0);
          if (r2 < r2min)
          {
            i     = ix;
            r2min = r2;
          }
        }
        assert(i >= 0);

        const vec3 &ipos = siteList[i].pos;
        assert(ipos.norm2() > 0.0);

        double t2 = get_wtime();
        dt_40 += t2 - t1;
        //t1 = t2;

        /* step 2:
         *  find site j, so that the triangle (origin, i, j)
         *  has minimal circumradius
         */
        r2min = HUGE;
        int j = -1;
        const real ipos2 = ipos.norm2();
        for (int ix = 0; ix < nSite; ix++)
          if (ix != i)
          {
            const vec3 &pos = siteList[ix].pos;
            const real r2ij = (ipos - pos).norm2();
            if (r2ij >= r2min) continue;

            const real area = (ipos%pos).norm2();
            const real pos2 =       pos .norm2();
            if (area < eps2*pos2) continue;
            const real r2 = ipos2*pos2*r2ij/area;
            if (r2 < r2min)
            {
              j     = ix;
              r2min = r2;
            }
          }
        assert(j >= 0);
        assert(j != i);
        const vec3 &jpos = siteList[j].pos;

        t2 = get_wtime();
        dt_44 += t2 - t1;
        //        t1 = t2;

        /* step 3:
         *  find site k, so that tetrahedron (origin, i,j,k)
         *  has minimal circumradius
         *  also site l that forms adjacent tetrahedron (origin, i, j, l)
         */
        int k = -1;
        int l = -1;
        real largek = -1e10, largel = -1e10;
        real     rk =   0.0,     rl =   0.0;
        vec3  cposk(0.0),     cposl(0.0);
        const Plane plane(ipos, jpos, true);
        for (int ix = 0; ix < nSite; ix++)
          if (ix != i && ix != j)
          {
            const vec3 &pos = siteList[ix].pos;
            const real ploc = plane(pos);
            const bool side  = ploc > 0.0;
            const real dist1 = pos*(pos + cposk) + largek;
            const real dist2 = pos*(pos + cposl) + largel;
            if (side && dist1 < 0.0)
            {
              if (ploc*ploc < eps2*pos.norm2()) continue;
              if (__abs(dist1) < eps) continue; 
              cposk = sphere(ipos, jpos, pos, rk)*(real)(-2.0);
              largek = 0.0;
              k = ix;
            }
            else if (!side && dist2 < 0.0)
            {
              if (ploc*ploc < eps2*pos.norm2()) continue;
              if (__abs(dist2) < eps) continue; 
              cposl = sphere(ipos, jpos, pos, rl)*(real)(-2.0);
              largel = 0.0;
              l = ix;
            }
          }
        assert(k >= 0);
        assert(l >= 0);
        const vec3 &kpos = siteList[k].pos;
        const vec3 &lpos = siteList[l].pos;
        assert(k != l);
        assert(k != i);
        assert(k != j);
        assert(l != i);
        assert(l != j);
        assert(__abs(plane(kpos)) > eps*kpos.abs());
        assert(__abs(plane(lpos)) > eps*lpos.abs());
        assert(plane(siteList[k].pos)*plane(siteList[l].pos) < 0.0);

        assert(rl > 0.0);
        assert(rk > 0.0);

        triangles.add(i,j, PackedInt2(k, tetraList.size()));
        triangles.add(j,k, PackedInt2(i, tetraList.size()));
        triangles.add(k,i, PackedInt2(j, tetraList.size()));
        tetraList.push_back(Tetrahedron(i,j,k, cposk*(real)(-0.5)));
        
        triangles.add(i,j, PackedInt2(l, tetraList.size()));
        triangles.add(j,l, PackedInt2(i, tetraList.size()));
        triangles.add(l,i, PackedInt2(j, tetraList.size()));
        tetraList.push_back(Tetrahedron(i,j,l, cposl*(real)(-0.5)));

#if 0
        fprintf(stderr," i= %d j= %d k= %d l= %d\n",
            i,j,k,l);
#endif

        vertexQueue.push_back(i);
        isVertexQueued[i] = true;

        vertexQueue.push_back(j);
        isVertexQueued[j] = true;

        vertexQueue.push_back(k);
        isVertexQueued[k] = true;

        vertexQueue.push_back(l);
        isVertexQueued[l] = true;

        dt_50 += get_wtime() - t1;

        if (!completeCell(siteList))
          return false;
        dt_60 += get_wtime() - t00;

        /* compute volume and area of each faces */

#if 0
        const int nnb = nbList.size();
        faceList.reserve(nnb);
        cellVolume = 0.0;
        Face face(0.0, 0.0);
        for (int i = 0; i < nnb; i++)
        {
          const int j = nbList[i];
#if 0   /* sanity check: pass over all tetrahedra to make sure they are all complete */
#error
          for (int i= 0; i < faceVtx[j].size(); i++)
            assert(isComplete(tetrahedra[faceVtx[j][i]], j));
#endif
#if 0
          if (!buildFace(faceVtx[j], cellVolume, face)) return false;
#endif
          faceList.push_back(face);
        }
#endif

        dt_70 += get_wtime() - t00;

#if 0  /* sanity check */
        const int nt = tetrahedra.size();
        for (int i = 0 ; i < nt-1; i++)
          for (int j = i+1; j < nt; j++)
            if (tetrahedra[i] == tetrahedra[j])
            {
              fprintf(stderr, "i= %d  (%d %d %d)  j= %d (%d %d %d) \n",
                  i, tetrahedra[i].vertex1(), tetrahedra[i].vertex2(), tetrahedra[i].vertex3(),
                  j, tetrahedra[j].vertex1(), tetrahedra[j].vertex2(), tetrahedra[j].vertex3());
              assert(0);
            }
#endif
        return true;
      }


      private:

      bool completeCell(const Site::Vector &siteList)
      {
        const double tX = get_wtime();
        const int nSites = siteList.size();
        for (int i = 0; i < nSites; i++)
        {
          sites [i] = Site(siteList[i].pos, i);
          vtxUse[i] = 1;
        }

#if 0
#define SWAP(i,n) {std::swap(sites[i], sites[--n]);  iMap[sites[i].idx] = i; iMap[sites[n].idx] = n; }
        std ::vector<int> iMap(nSites);
        for (int i = 0; i < nSites; i++)
          sites[i].idx = i;

#else
#define SWAP(i,n) {}
#endif
        while (!vertexQueue.empty())
        {
          /* step 4.2:
           *  extract vertex from the list
           */
          const int iVertex = vertexQueue.front();
          vertexQueue.pop_front();

          /* step 4.3-4.4:
           *  find a tetrahedron (i, iVertex, jVertex, kVertex) with at least one
           *  incomplete face
           */

          int jVertex = triangles[iVertex];
          assert(jVertex >= 0);
          assert(!triangles(iVertex, jVertex).empty());
          int kVertex = triangles(iVertex, jVertex).begin()->first();

          const int endVertex = kVertex;
          const double tY = get_wtime();
          assert(vtxList.empty());

          vtxUse[iVertex] = 0;
          vtxUse[jVertex] = 0;
          vtxList.push(iVertex);
          vtxList.push(jVertex);


          int nSites_loc = siteList.size();
          SWAP(iMap[iVertex], nSites_loc);
          SWAP(iMap[jVertex], nSites_loc);

          vtxPos.clear();
          vtxPos.push_back(tetraList[triangles(iVertex, jVertex).begin()->second()].centre());

          int cnt = 0;
          while (jVertex != endVertex)
          {
            if (cnt++ > 100) assert(0);
            const vec3 &ipos = siteList[iVertex].pos;
            const vec3 &jpos = siteList[jVertex].pos;
            const vec3 &kpos = siteList[kVertex].pos;

            const Plane plane(ipos, jpos);
            const real kloc  = plane(kpos);
            const int  sideK = kloc > 0.0;
            assert(kloc*kloc > eps2*kpos.norm2());
            real largeNum = -1e10;
            vec3  cpos(0.0);
            int lVertex = -1;

            const List &list = triangles(iVertex, jVertex);
            for (List::const_iterator it = list.begin(); it != list.end(); it++)
            {
              const int i = it->first();
              const vec3 &pos = siteList[i].pos;
              const int  side = plane(pos) > 0.0;
              if (side^sideK && vtxUse[i])
              {
                lVertex = i;
                cpos    = tetraList[it->second()].centre();
                break;
              }
            }

            if (lVertex < 0)
            {
              /* vertex is not found in tetrahedra list, search the siteList*/
              myNMX += nSites;

              for (int i = 0; i < nSites_loc; i++)
              {
                const Site   &s = sites[i];
                const real  loc = plane(s.pos);
                const int  side = loc > 0.0;
                const real dist = s.pos*(s.pos + cpos) + largeNum;
                flop += 20;
                if (dist < -eps && side^sideK && vtxUse[s.idx])
                {
                  if (loc*loc < eps2*s.pos.norm2()) continue;
                  flop += 92;
                  real radius = 0.0;
                  const vec3 _cpos = sphere(ipos, jpos, s.pos, radius)*(real)(-2.0);
                  if (radius > 0.0) 
                  {
                    cpos = _cpos;
                    largeNum = 0.0;
                    lVertex  = s.idx;
                  }
                }
              }
              assert(lVertex >= 0);
#if 0      /* sanity check */
              {
                List::const_iterator begin = triangles(iVertex, jVertex).begin();
                List::const_iterator end   = triangles(iVertex, jVertex).end  ();
                for (List::const_iterator it = begin; it != end; it++)
                  assert(it->first() != (unsigned int)lVertex); 
                begin = triangles(jVertex, lVertex).begin();
                end   = triangles(jVertex, lVertex).end  ();
                for (List::const_iterator it = begin; it != end; it++)
                  assert(it->first() != (unsigned int)iVertex); 
                begin = triangles(lVertex, iVertex).begin();
                end   = triangles(lVertex, iVertex).end  ();
                for (List::const_iterator it = begin; it != end; it++)
                  assert(it->first() != (unsigned int)jVertex); 
              }
#endif
#if 0
              assert(triangles(iVertex,jVertex).size() < 2);
              assert(triangles(jVertex,lVertex).size() < 2);
              assert(triangles(lVertex,iVertex).size() < 2);
#endif
              triangles.add(iVertex, jVertex, PackedInt2(lVertex, tetraList.size()));
              triangles.add(jVertex, lVertex, PackedInt2(iVertex, tetraList.size()));
              triangles.add(lVertex, iVertex, PackedInt2(jVertex, tetraList.size()));
              cpos *= (real)(-0.5);
              tetraList.push_back(Tetrahedron(iVertex, jVertex, lVertex, cpos));
            }
            else
            {
#if 0       /* sanity check */
              bool flag = false;
              List::const_iterator begin = triangles(iVertex, jVertex).begin();
              List::const_iterator end   = triangles(iVertex, jVertex).end  ();
              for (List::const_iterator it = begin; it != end; it++)
                flag |= it->first() == (unsigned int)lVertex; 
              assert(flag);
              flag = false;
              begin = triangles(jVertex, lVertex).begin();
              end   = triangles(jVertex, lVertex).end  ();
              for (List::const_iterator it = begin; it != end; it++)
                flag |= it->first() == (unsigned int)iVertex;
              assert(flag);
              flag = false;
              begin = triangles(lVertex, iVertex).begin();
              end   = triangles(lVertex, iVertex).end  ();
              for (List::const_iterator it = begin; it != end; it++)
                flag |= it->first() == (unsigned int)jVertex;
              assert(flag);
#endif
            }

            /* vertex must be found at this point */

            if (!isVertexQueued[lVertex])
            {
              vertexQueue.push_back(lVertex);
              isVertexQueued[lVertex] = true;
            }

            assert(lVertex >= 0);
            if (vtxUse[lVertex] == 1)
            {
              vtxList.push(lVertex);
              vtxUse[lVertex] = 0;
              SWAP(iMap[lVertex], nSites_loc);
            }

            vtxPos.push_back(cpos);

            kVertex = jVertex;
            jVertex = lVertex;
          }
          dt_20 += get_wtime() - tY;

          while(!vtxList.empty())
          {
            vtxUse[vtxList.top()] = 1;
            vtxList.pop();
          }

          if (!vtxPos.empty())
          {
#if 1
            const int nvtx = vtxPos.size();
            vec3 cpos(0.0);
            for (int i = 0; i < nvtx; i++)
              cpos += vtxPos[i];
            cpos *= 1.0/(real)nvtx;

            vec3 area(0.0);
            real vol (0.0);
            vtxPos.push_back(vtxPos[0]);
            for (int i = 0; i < nvtx; i++)
            {
              const vec3 &v1 = vtxPos[i  ];
              const vec3 &v2 = vtxPos[i+1];
              area  += (v1-cpos)%(v2-cpos);
              vol   += __abs(cpos*(v1%v2));
            }
            cellVolume += vol*(1.0/6.0);
            const real A = area.abs();
            if (A > 0.0)
            {
              faceList.push_back(Face(area/A, 0.5*A));
              nbList.push_back(iVertex);
            }
#else
            Face f;
            if (traceFace(vtxPos, cellVolume, f))
            {
              faceList.push_back(f);
              nbList.push_back(iVertex);
            }
#endif
          }


          vertexCompleted[iVertex] = true;
        }
        dt_30 += get_wtime() - tX;
        return true;
      }


      vec3 sphere(const vec3 &ip, const vec3 &jp, const vec3 &kp, real &radius) const
      {

        real XI, YI, ZI, DI;
        real XJ, YJ, ZJ, DJ;
        real XK, YK, ZK, DK;
        real XC, YC, ZC, RR, D;

        XI = ip.x;
        YI = ip.y;
        ZI = ip.z;
        DI = ip.norm2();

        XJ = jp.x;
        YJ = jp.y;
        ZJ = jp.z;
        DJ = jp.norm2();

        XK = kp.x;
        YK = kp.y;
        ZK = kp.z;
        DK = kp.norm2();

        D=XI*YJ*ZK+XJ*YK*ZI+XK*YI*ZJ-XK*YJ*ZI-XI*YK*ZJ-XJ*YI*ZK;

#if 0
        if (__abs(D*D) < eps*DI*DJ*DK) {radius = -1.0; return vec3(0.0);}
#else
        if (D == 0.0) {radius = -1.0; return vec3(0.0);}
#endif

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

      real triangle(const vec3 &jpos, const vec3 &kpos) const
      {
        real radius;
        sphere(jpos, kpos, Plane(jpos, kpos).n, radius);
        return radius;
      }

      bool traceFace(const std::vector<vec3> &vtxPos, real &volume, Face &face)
      {
        const int n = vtxPos.size();
        angle_vec_pair.resize(n);
        vec3 cpos(0.0);
        for (int i= 0; i < n; i++)
        {
          const vec3 &jpos = vtxPos[i];
          cpos += jpos;
          angle_vec_pair[i].second = jpos;
        }
        cpos *= 1.0/(real)n;

        const vec3 &posA = angle_vec_pair[0].second - cpos;
        //assert(posA.norm2() > 0.0);
        if (posA.norm2() == 0)
          return false;
        const vec3 unitA = posA * (1.0/posA.abs());

        int iv = 1;
#if 0
        const real eps4a = eps2*eps2*posA.norm2();
        for (iv = 1; iv < n; iv++)
        {
          const vec3 &posB = angle_vec_pair[iv].second - cpos;
          const real q = (posA%posB).norm2();
          if (q*q > eps4a*posB.norm2())
            break;
        }
        assert(iv < n);
        if (iv == n) 
          return false;
#endif
        const vec3 &posB = angle_vec_pair[iv].second - cpos;
        vec3  unitN  = posA%posB;
//        assert(unitN.norm2() > 0.0);
        if (unitN.norm2() == 0.0)
          return false;
        unitN *= 1.0/unitN.abs();
        for (int j = 0; j < n; j++)
        {
          const vec3  jpos = angle_vec_pair[j].second - cpos;
          const real    r2 = jpos.norm2();
          if (r2 == 0.0) return false;
          assert(r2 > 0.0);
          const real cos    =  unitA * jpos;
          const real sin    = (unitA % jpos) * unitN;
#if 1
          const real cos2   =  cos*__abs(cos) * (1.0/r2);
          angle_vec_pair[j] = std::make_pair(sin >  0.0 ? -cos2 : 2.0+cos2, jpos);
#else
          angle_vec_pair[j] = std::make_pair(std::atan2(sin, cos), jpos);
#endif
        }
        std::sort(angle_vec_pair.begin(), angle_vec_pair.end(), cmp_data<real, vec3>());
        angle_vec_pair.push_back(angle_vec_pair[0]);

        /* comput area & volume */
        real vol  = 0.0;
        vec3 area(0.0);
        const vec3 &v0 = cpos;
        for (int i = 0; i < n; i++)
        {
          const vec3 &v1 = angle_vec_pair[i  ].second;
          const vec3 &v2 = angle_vec_pair[i+1].second;
          area  += v1%v2;
          vol   += __abs(v0*((v1+v0)%(v2+v0)));
        }
        volume += vol*(1.0/6.0);

        face = Face(unitN, 0.5*area.abs());
        return true;
      }

    };




}

#endif /* __VOROCELL_H__ */
