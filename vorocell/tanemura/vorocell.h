#ifndef __VOROCELL_H__
#define __VOROCELL_H__

/* Tanemura's algorithm */

#include <cassert>
#include <vector>
#include <deque>
#include <algorithm>
#include "vector3.h"

template<class T>
inline T __min(const T a, const T b) {return a < b ? a : b;}

template<class T>
inline T __max(const T a, const T b) {return a > b ? a : b;}

template<class T>
inline T __abs(const T a) {return a < T(0.0) ? -a : a;}

inline float __atan2(float y, float x)
{
  float t0, t1, t3, t4;

  t3 = __abs(x);
  t1 = __abs(y);
  t0 = __max(t3, t1);
  t1 = __min(t3, t1);
  t3 = float(1) / t0;
  t3 = t1 * t3;

  t4 = t3 * t3;
  t0 =         - float(0.013480470);
  t0 = t0 * t4 + float(0.057477314);
  t0 = t0 * t4 - float(0.121239071);
  t0 = t0 * t4 + float(0.195635925);
  t0 = t0 * t4 - float(0.332994597);
  t0 = t0 * t4 + float(0.999995630);
  t3 = t0 * t3;

  t3 = (abs(y) > abs(x)) ? float(1.570796327) - t3 : t3;
  t3 = (x < 0) ?  float(3.141592654) - t3 : t3;
  t3 = (y < 0) ? -t3 : t3;

  return t3;
}


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
    long long idx;
    Site() {}
    Site(const vec3 &_pos, const long long _idx) : pos(_pos), idx(_idx) {}
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

  template<const int N>
    struct ConnectivityMatrix2D
    {
      private:
        int matrix[N][N];

      public:

        ConnectivityMatrix2D() 
        {
          clear();
        }

        void clear()
        {
          for (int j = 0; j < N; j++)
            for (int i = 0; i < N; i++)
              matrix[j][i] = 0;
        }

        /* i: x
         * j: y */
        int   operator()(const int i, const int j) const { return matrix[__min(i,j)][__max(i,j)]; }
        int&  operator()(const int i, const int j)       { return matrix[__min(i,j)][__max(i,j)]; }
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
    Face(const vec3 &n, const real A) : norm(n), area(A) {}
  };

  template<const int N>
    struct Cell
    {
      typedef ConnectivityMatrix2D<N>  Triangles;

      private:
      Triangles triangles;
      real eps;

      std::deque < int>   vertexQueue    ;
      std::vector<bool> isVertexQueued   ;
      std::vector<bool>   vertexCompleted;

      std::vector<Tetrahedron> tetrahedra;
      std::vector<   int     >     nbList;
      std::vector<  Face     >   faceList;
      std::vector< std::pair<real, vec3> > angle_vec_pair;
      std::deque<int> incompleteTetra;

      Array<int,  N> faceVtx[N];

      real cellVolume;

      public:
      Cell(const real _eps = 1.0e-11) : eps(_eps)
      {
        angle_vec_pair.reserve(N);
      }

      int nb() const {return nbList.size();}
      real volume() const {return cellVolume;}

      void build(const Site::Vector &siteList)
      {
        const double t00 = get_wtime();

        assert((int)siteList.size() <= N);
        clear(siteList.size());
        const int nSite = siteList.size();

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
        t1 = t2;

        /* step 2:
         *  find site j, so that the triangle (origin, i, j)
         *  has minimal circumradius
         */
        r2min = HUGE;
        int j = -1;
        for (int ix = 0; ix < nSite; ix++)
          if (ix != i)
          {
            const vec3 &pos = siteList[ix].pos;
            const real   r2 = triangle(ipos, pos);
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
        t1 = t2;

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
            if (__abs(ploc) < eps*pos.abs()) continue;
            const bool side  = ploc > 0.0;
            const real dist1 = pos*(pos + cposk) + largek;
            const real dist2 = pos*(pos + cposl) + largel;
            if (side && dist1 < 0.0)
            {
              cposk = sphere(ipos, jpos, pos, rk)*(real)(-2.0);
              largek = 0.0;
              k = ix;
            }
            else if (!side && dist2 < 0.0)
            {
              cposl = sphere(ipos, jpos, pos, rl)*(real)(-2.0);
              largel = 0.0;
              l = ix;
            }
          }
#if 0
        fprintf(stderr, "i= %d j= %d k= %d l= %d\n", i,j,k,l);
        fprintf(stderr, "i= %d: %g %g %g\n", i, ipos.x, ipos.y, ipos.z);
        fprintf(stderr, "j= %d: %g %g %g\n", j, jpos.x, jpos.y, jpos.z);
#endif
        assert(k >= 0);
        assert(l >= 0);
        const vec3 &kpos = siteList[k].pos;
        const vec3 &lpos = siteList[l].pos;
#if 0
        fprintf(stderr, "k= %d: %g %g %g  rk=%g\n", k, kpos.x, kpos.y, kpos.z, rk);
        fprintf(stderr, "l= %d: %g %g %g  rl=%g\n", l, lpos.x, lpos.y, lpos.z, rl);
#endif
        assert(k != l);
        assert(k != i);
        assert(k != j);
        assert(l != i);
        assert(l != j);
        assert(__abs(plane(kpos)) > eps*kpos.abs());
        assert(__abs(plane(lpos)) > eps*lpos.abs());
        assert(plane(siteList[k].pos)*plane(siteList[l].pos) < 0.0);

#if 0
        fprintf(stderr, "cposl= %g %g %g \n", cposl.x, cposl.y, cposl.z);
        fprintf(stderr, "cposl= %g %g %g \n", cposl.x, cposl.y, cposl.z);
#endif
        assert(rl > 0.0);
        assert(rk > 0.0);

        faceVtx[i].push_back(tetrahedra.size());
        faceVtx[j].push_back(tetrahedra.size());
        faceVtx[k].push_back(tetrahedra.size());
        tetrahedra.push_back(Tetrahedron(i,j,k, cposk*(real)(-0.5)));

        assert(triangles(i,j) < 2);
        assert(triangles(j,k) < 2);
        assert(triangles(i,k) < 2);

        triangles(i,j)++;
        triangles(i,k)++;
        triangles(j,k)++;

        faceVtx[i].push_back(tetrahedra.size());
        faceVtx[j].push_back(tetrahedra.size());
        faceVtx[l].push_back(tetrahedra.size());
        tetrahedra.push_back(Tetrahedron(i,j,l, cposl*(real)(-0.5)));

        assert(triangles(i,j) < 2);
        assert(triangles(i,l) < 2);
        assert(triangles(j,l) < 2);
        triangles(i,j)++;
        triangles(i,l)++;
        triangles(j,l)++;

        vertexQueue.push_back(i);
        isVertexQueued[i] = true;

        vertexQueue.push_back(j);
        isVertexQueued[j] = true;

        vertexQueue.push_back(k);
        isVertexQueued[k] = true;

        vertexQueue.push_back(l);
        isVertexQueued[l] = true;

        dt_50 += get_wtime() - t1;

        completeCell(siteList);
        dt_60 += get_wtime() - t00;

        /* compute volume and area of each faces */

        const int nnb = nbList.size();
        faceList.reserve(nnb);
        cellVolume = 0.0;
        for (int i = 0; i < nnb; i++)
        {
          const int j = nbList[i];
#if 0   /* sanity check: pass over all tetrahedra to make sure they are all complete */
          for (int i= 0; i < faceVtx[j].size(); i++)
            assert(isComplete(tetrahedra[faceVtx[j][i]], j));
#endif
          faceList.push_back(buildFace(faceVtx[j], cellVolume));
        }

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
      }

      private:

      bool isComplete(const Tetrahedron &t, const int vtx) const
      {
        const std::pair<int,int> vpair = t.pair(vtx);
        return 
          triangles(vtx, vpair.first ) == 2 &&
          triangles(vtx, vpair.second) == 2;
      }

      void clear(const int nSites)
      {
        triangles  .clear();
        vertexQueue.clear();

        vertexCompleted.resize(nSites);
        isVertexQueued .resize(nSites);
        for (int i = 0; i < nSites; i++)
          isVertexQueued[i] = vertexCompleted[i] = false;

        tetrahedra.clear();
        nbList    .clear();
        faceList  .clear();

        for (int i = 0; i < N; i++)
          faceVtx[i].clear();
      }

      private:

      void completeCell(const Site::Vector &siteList)
      {
        const double tX = get_wtime();
        assert(incompleteTetra.empty());
        const int nSites = siteList.size();

        while (!vertexQueue.empty())
        {
          /* step 4.2:
           *  extract vertex from the list
           */
          const int iVertex = vertexQueue.front();
          vertexQueue.pop_front();

#if 0  /* sanity check: the vertex haven't yet registered in nbList */
          for (std::vector<int>::const_iterator it = nbList.begin(); it != nbList.end(); it++)
            assert(*it !=  iVertex);
#endif
          nbList.push_back(iVertex);

          /* step 4.3:
           *  the vertex is completed, proceed to the next one
           */
          assert(incompleteTetra.empty());
          for (int i = 0; i < faceVtx[iVertex].size(); i++)
            if (!isComplete(tetrahedra[faceVtx[iVertex][i]], iVertex))
              incompleteTetra.push_back(faceVtx[iVertex][i]);

          if (!incompleteTetra.empty())
            assert(!vertexCompleted[iVertex]);

          /* step 4.4:
           *  find a tetrahedron (i, iVertex, jVertex, kVertex) with at least one
           *  incomplete face
           */

          /* otherwise, find the face that lacks adjacent tetrahedron */

          //          const double tAA = get_wtime();
          if (!incompleteTetra.empty())
            incomplT++;

          while(!incompleteTetra.empty())
          {
            const Tetrahedron &t = tetrahedra[incompleteTetra.front()];
            incompleteTetra.pop_front();
            if (isComplete(t, iVertex)) continue;
            incompl++;

            const std::pair<int,int> vpair = t.pair(iVertex);
            int jVertex = vpair.first;
            int kVertex = vpair.second;
            if (triangles(iVertex, jVertex) != 1)
              std::swap(jVertex, kVertex);
            assert(triangles(iVertex, jVertex) == 1);
            assert(triangles(iVertex, kVertex) == 2);

            /* step 4.5-4.7 */

            //            const double t00 = get_wtime();
            while(triangles(iVertex, jVertex) != 2)
            {
              assert(triangles(iVertex, jVertex) < 2);
              /* step 4.5 - 4.6: 
               *  search a vertex on the opposite side of the kVertex 
               *  (in the half-space bounded by ijFace that does not contain kVertex)
               */
              const vec3 &ipos = siteList[iVertex].pos;
              const vec3 &jpos = siteList[jVertex].pos;
              const vec3 &kpos = siteList[kVertex].pos;
              const Plane plane(ipos, jpos);
              const int  sideK = plane(kpos) > 0.0;
              real largeNum = -1e10;
              vec3  cpos(0.0);
              int lVertex = -1;

              /* hot-spot: finding 4th vertex of the new tetrahedron */

              //              const double tA = get_wtime();
              for (int i = 0; i < nSites; i++)
              {
                const vec3 &pos = siteList[i].pos;
                const int  side = plane(pos) > 0.0;
                const real dist = pos*(pos + cpos) + largeNum;
                const bool skip = vertexCompleted[i] || 
                  (i == iVertex) || (i == jVertex) || (i == kVertex);
             
                 
                if (dist < 0.0 && side^sideK && !skip)
                {
                 // if (triangles(iVertex,i)==2 || triangles(jVertex,i)==2) continue;

                  real radius = 0.0;
                  const vec3 _cpos = sphere(ipos, jpos, pos, radius)*(real)(-2.0);
                  //assert(radius > 0.0);
                  if (radius < 0.0) continue;
                  cpos = _cpos;
                  largeNum = 0.0;
                  lVertex = i;
                }
              }
              //              dt_00 += get_wtime() - tA;
              if (lVertex < 0)
              {
                fprintf(stderr, " nb= %d\n", (int)nbList.size());
              }
              assert(lVertex >= 0);

              /* step 4.7:
               *  register new tetrahedron (iVertex, jVertex, lVertex) 
               */
              faceVtx[iVertex].push_back(tetrahedra.size());
              faceVtx[jVertex].push_back(tetrahedra.size());
              faceVtx[lVertex].push_back(tetrahedra.size());
              tetrahedra.push_back(Tetrahedron(iVertex, jVertex, lVertex, cpos*(real)(-0.5)));

              if (!isVertexQueued[lVertex])
              {
                vertexQueue.push_back(lVertex);
                isVertexQueued[lVertex] = true;
              }

#if 0
              fprintf(stderr, "(%d,%d)= %d\n", __min(iVertex,jVertex), __max(iVertex,jVertex), triangles(iVertex,jVertex));
              fprintf(stderr, "(%d,%d)= %d\n", __min(iVertex,lVertex), __max(iVertex,lVertex), triangles(iVertex,lVertex));
              fprintf(stderr, "(%d,%d)= %d\n", __min(jVertex,lVertex), __max(jVertex,lVertex), triangles(jVertex,lVertex));
              fprintf(stderr, " --- \n");
#endif
              assert(triangles(iVertex, jVertex) < 2);
              assert(triangles(iVertex, lVertex) < 2);
              assert(triangles(jVertex, lVertex) < 2);
              triangles(iVertex, jVertex)++;
              triangles(iVertex, lVertex)++;
              triangles(jVertex, lVertex)++;

#if 0         /* sanity check */
              bool complete = true; 
              for (int i = 0; i < faceVtx[kVertex].size(); i++)
                complete &= isComplete(tetrahedra[faceVtx[kVertex][i]], kVertex);
              if (complete) vertexCompleted[kVertex] = true;

              complete = true; 
              for (int i = 0; i < faceVtx[lVertex].size(); i++)
                complete &= isComplete(tetrahedra[faceVtx[lVertex][i]], lVertex);
              if (complete) vertexCompleted[lVertex] = true;

              complete = true; 
              for (int i = 0; i < faceVtx[jVertex].size(); i++)
                complete &= isComplete(tetrahedra[faceVtx[jVertex][i]], jVertex);
              if (complete) vertexCompleted[jVertex] = true;

              complete = true; 
              for (int i = 0; i < faceVtx[jVertex].size(); i++)
                complete &= isComplete(tetrahedra[faceVtx[iVertex][i]], iVertex);
              if (complete) break;

#endif

              kVertex = jVertex;
              jVertex = lVertex;
              //              dtA += get_wtime() - tB;
            }
            //            dt_10 += get_wtime() - t00;
          }
          //          dt_20 += get_wtime() - tAA;

          /* step 4.8: */
#if 0   /* sanity check: pass over all tetrahedra to make sure they are all complete */
          for (int i= 0; i < faceVtx[iVertex].size(); i++)
            assert(isComplete(tetrahedra[faceVtx[iVertex][i]], iVertex));
#endif
          vertexCompleted[iVertex] = true;
        }
        dt_30 += get_wtime() - tX;
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

      Face buildFace(const Array<int, N> vtxList, real &volume)
      {
        const int n = vtxList.size();
        vec3 cpos(0.0);
        for (int i= 0; i < n; i++)
          cpos += tetrahedra[vtxList[i]].centre();
        cpos *= 1.0/(real)n;

        const vec3 &posA = tetrahedra[vtxList[0]].centre() - cpos;
        const vec3 &posB = tetrahedra[vtxList[1]].centre() - cpos;
        const vec3 &unitA = posA * (1.0/posA.abs());
        const vec3 &unitB = posB * (1.0/posB.abs());
        const vec3  unitN  = unitA%unitB;
        angle_vec_pair.clear();
#if 0
        angle_vec_pair.push_back(std::make_pair(unitA*unitA, posA));
        angle_vec_pair.push_back(std::make_pair(unitA*unitB, posB));
        for (int j = 2; j < n; j++)
        {
          const vec3 &jpos = tetrahedra[vtxList[j]].centre() - cpos;
          const vec3 unitj = jpos * (1.0/jpos.abs());
          const float cos =  unitA * unitj;
          const float sin = (unitA % unitj) * unitN;
          if (sin >= 0.0f)  angle_vec_pair.push_back(std::make_pair(  -cos, jpos));
          else              angle_vec_pair.push_back(std::make_pair(10+cos, jpos));
        }
#error
#else
        angle_vec_pair.push_back(std::make_pair(0.0, posA));
        for (int j = 1; j < n; j++)
        {
          const vec3 &jpos = tetrahedra[vtxList[j]].centre() - cpos;
          const real x =  unitA * jpos;
          const real y = (unitA % jpos) * unitN;
          const real dot = jpos * unitN;
          assert(__abs(dot) < jpos.abs()*1.0e-10);
          angle_vec_pair.push_back(std::make_pair(std::atan2(y, x), jpos));
        }
#endif
        assert((int)angle_vec_pair.size() == n);
        std::sort(angle_vec_pair.begin(), angle_vec_pair.end(), cmp_data<real, vec3>());
        angle_vec_pair.push_back(angle_vec_pair[0]);

        /* comput area & volume */
        real area = 0.0;
        real vol  = 0.0;
        const vec3 &v0 = cpos;
        for (int i = 0; i < n; i++)
        {
          const vec3 &v1 = angle_vec_pair[i  ].second;
          const vec3 &v2 = angle_vec_pair[i+1].second;
          area   += (v1%v2).abs();
          vol    += __abs(v0*((v1+v0)%(v2+v0)));
        }
        area   *= 0.5;
        volume += vol*(1.0/6.0);

        return Face(unitN, area);
      }

    };




}

#endif /* __VOROCELL_H__ */
