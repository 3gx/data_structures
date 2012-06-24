#ifndef __VOROCELL_DEGERATE_H__
#define __VOROCELL_DEGERATE_H__

/* Tanemura's algorithm */

namespace VoronoiDegenerate
{
  typedef double real;
  typedef vector3<real>  vec3;

  typedef Voronoi::Site Site;
#if 0
  struct Site
  {
    typedef std::vector<Site> Vector;
    vec3 pos;
    long long idx;
    Site() {}
    Site(const vec3 &_pos, const long long _idx) : pos(_pos), idx(_idx) {}
  };
#endif

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
      real r;
    public:
      Tetrahedron(const int _v1, const int _v2, const int _v3, const vec3 &_c, const real _r) : 
        v1(_v1), v2(_v2), v3(_v3), c(_c), r(_r) {}

      int vertex1() const {return v1;}
      int vertex2() const {return v2;}
      int vertex3() const {return v3;}
      const vec3& centre() const { return c; }
      real radius() const { return r;}

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
      int nface() const {return faceList.size();}
      real volume() const {return cellVolume;}

      bool build(const Site::Vector &siteList)
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
            if ((ipos%pos).norm2() == 0.0) continue;
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
        const real eps2 = eps*eps;
        for (int ix = 0; ix < nSite; ix++)
          if (ix != i && ix != j)
          {
            const vec3 &pos = siteList[ix].pos;
            const real ploc = plane(pos);
            if (ploc*ploc < eps2*pos.norm2()) continue;
            const bool side  = ploc > 0.0;
            const real dist1 = pos*(pos + cposk) + largek;
            const real dist2 = pos*(pos + cposl) + largel;
            if (side && dist1 < 0.0)
            {
              if (__abs(dist1) < eps) continue; 
              cposk = sphere(ipos, jpos, pos, rk)*(real)(-2.0);
              largek = 0.0;
              k = ix;
            }
            else if (!side && dist2 < 0.0)
            {
              if (__abs(dist2) < eps) continue; 
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
        assert(l >= 0);
        const vec3 &lpos = siteList[l].pos;
//        fprintf(stderr, "l= %d: %g %g %g  rl=%g\n", l, lpos.x, lpos.y, lpos.z, rl);

        const vec3 &kpos = siteList[k].pos;
        assert(k >= 0);
//        fprintf(stderr, "k= %d: %g %g %g  rk=%g\n", k, kpos.x, kpos.y, kpos.z, rk);

#if 0
        for (int m = 0; m < 10; m++)
        {
          real radius;
          sphere(ipos, jpos, siteList[m].pos, radius);
          fprintf(stderr, "m= %d: r= %g side= %g \n",
              m, radius, plane(siteList[m].pos));
        }
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
        fprintf(stderr, "cposk= %g %g %g \n", cposk.x, cposk.y, cposk.z);
#endif
        assert(rl > 0.0);
        assert(rk > 0.0);
        if (rl < rk)
        {
          std::swap(k,l);
          std::swap(rk,rl);
          std::swap(cposk, cposl);
        }

        faceVtx[i].push_back(tetrahedra.size());
        faceVtx[j].push_back(tetrahedra.size());
        faceVtx[k].push_back(tetrahedra.size());
        tetrahedra.push_back(Tetrahedron(i,j,k, cposk*(real)(-0.5), rk));

#if 0
        assert(triangles(i,j) < 2);
        assert(triangles(j,k) < 2);
        assert(triangles(i,k) < 2);

        triangles(i,j)++;
        triangles(i,k)++;
        triangles(j,k)++;

        faceVtx[i].push_back(tetrahedra.size());
        faceVtx[j].push_back(tetrahedra.size());
        faceVtx[l].push_back(tetrahedra.size());
        tetrahedra.push_back(Tetrahedron(i,j,l, cposl*(real)(-0.5), rl));

        assert(triangles(i,j) < 2);
        assert(triangles(i,l) < 2);
        assert(triangles(j,l) < 2);
        triangles(i,j)++;
        triangles(i,l)++;
        triangles(j,l)++;
#endif

        vertexQueue.push_back(i);
        isVertexQueued[i] = true;

        vertexQueue.push_back(j);
        isVertexQueued[j] = true;

        vertexQueue.push_back(k);
        isVertexQueued[k] = true;

#if 0
        vertexQueue.push_back(l);
        isVertexQueued[l] = true;
#endif

        dt_50 += get_wtime() - t1;

        if (!completeCell(siteList))
          return false;
        dt_60 += get_wtime() - t00;

        /* compute volume and area of each faces */

        const int nnb = nbList.size();
        faceList.reserve(nnb);
        cellVolume = 0.0;
//        fprintf(stderr, " -- nnb= %d\n", nnb);
        Face face(0.0, 0.0);
        for (int i = 0; i < nnb; i++)
        {
          const int j = nbList[i];
#if 0   /* sanity check: pass over all tetrahedra to make sure they are all complete */
          for (int k= 0; k < faceVtx[j].size()-1; k++)
            for (int l = k+1; l < faceVtx[j].size(); l++)
            assert(tetrahedra[faceVtx[j][k]] != tetrahedra[faceVtx[j][l]]);
#endif
          if (buildFace(faceVtx[j], cellVolume, face))
            faceList.push_back(face);
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
        return true;
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

      bool completeCell(const Site::Vector &siteList)
      {
        const double tX = get_wtime();
        const int nSites = siteList.size();
        std::vector<int> vtxCount(nSites, 0);
        std::stack<int>  vtxUsed;

        while (!vertexQueue.empty())
        {
          /* step 4.2:
           *  extract vertex from the list
           */
          const int iVertex = vertexQueue.front();
          vertexQueue.pop_front();

          nbList.push_back(iVertex);

          /* step 4.3:
           *  the vertex is completed, proceed to the next one
           */

          assert(!faceVtx[iVertex].empty());

          const int it = faceVtx[iVertex][0];
          const Tetrahedron &t = tetrahedra[it];

          faceVtx[iVertex].clear();
          faceVtx[iVertex].push_back(it);

          const std::pair<int,int> vpair = t.pair(iVertex);
          int jVertex = vpair.first;
          int kVertex = vpair.second;

          /* step 4.5-4.7 */

          vtxCount[iVertex]++; vtxUsed.push(iVertex);
          vtxCount[jVertex]++; vtxUsed.push(jVertex);
          const double tY = get_wtime();
          while(vtxCount[jVertex] != 2)
          {
            assert(vtxCount[jVertex] < 2);
            /* step 4.5 - 4.6: 
             *  search a vertex on the opposite side of the kVertex 
             *  (in the half-space bounded by ijFace that does not contain kVertex)
             */
            const vec3 &ipos = siteList[iVertex].pos;
            const vec3 &jpos = siteList[jVertex].pos;
            const vec3 &kpos = siteList[kVertex].pos;
            const Plane plane(ipos, jpos, true);
            const int  sideK = plane(kpos) > 0.0;
            assert(__abs(plane(kpos)) >= eps*kpos.abs());
            real largeNum = -1e10;
            vec3  cpos(0.0);
            int lVertex = -1;
            real radius = -1;

            /* hot-spot: finding 4th vertex of the new tetrahedron */

            vtxCount[iVertex] = -1-vtxCount[iVertex];
            vtxCount[jVertex] = -1-vtxCount[jVertex];
            vtxCount[kVertex] = -1-vtxCount[kVertex];
            for (int i = 0; i < nSites; i++)
            {
              const vec3 &pos = siteList[i].pos;
              const int  side = plane(pos) > 0.0;
              const real dist = pos*(pos + cpos) + largeNum;
              const bool  use = vtxCount[i] >= 0;

              if (dist < 0.0 && side^sideK && use)
              {
                if (__abs(dist) < eps) continue; 
                if (__abs(plane(pos)) < eps*pos.abs()) continue;

                const vec3 _cpos = sphere(ipos, jpos, pos, radius)*(real)(-2.0);
                if (radius < 0.0) continue;
                cpos = _cpos;
                largeNum = 0.0;
                lVertex = i;
              }
            }
            vtxCount[iVertex] = -1-vtxCount[iVertex];
            vtxCount[jVertex] = -1-vtxCount[jVertex];
            vtxCount[kVertex] = -1-vtxCount[kVertex];
            assert(lVertex >= 0);
            assert(radius > 0.0);

            /* step 4.7:
             *  register new tetrahedron (iVertex, jVertex, lVertex) 
             */
            if (!isVertexQueued[lVertex])
            {
              vertexQueue.push_back(lVertex);
              isVertexQueued[lVertex] = true;
            }
           
            faceVtx[iVertex].push_back(tetrahedra.size());
            if (!vertexCompleted[jVertex]) 
              faceVtx[jVertex].push_back(tetrahedra.size());
            if (!vertexCompleted[lVertex]) 
              faceVtx[lVertex].push_back(tetrahedra.size());
            tetrahedra.push_back(Tetrahedron(iVertex, jVertex, lVertex, cpos*(real)(-0.5), radius));
            assert(vtxCount[jVertex] < 2);
            if (vtxCount[lVertex] == 2)
            {
              kVertex = jVertex;
              jVertex = lVertex;
              continue;
            }
            assert(vtxCount[lVertex] < 2);

            if (vtxCount[lVertex] == 0) 
              vtxUsed.push(lVertex);

            vtxCount[jVertex]++;
            vtxCount[lVertex]++;

            kVertex = jVertex;
            jVertex = lVertex;
          }
          while (!vtxUsed.empty())
          {
            vtxCount[vtxUsed.top()] = 0;
            vtxUsed.pop();
          }
          dt_20 += get_wtime() - tY;

          /* step 4.8: */
#if 0   /* sanity check: pass over all tetrahedra to make sure they are all complete */
          for (int i= 0; i < faceVtx[iVertex].size(); i++)
            assert(isComplete(tetrahedra[faceVtx[iVertex][i]], iVertex));
#endif
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

      bool buildFace(const Array<int, N> vtxList, real &volume, Face &face)
      {
        const int n = vtxList.size();
        vec3 cpos(0.0);
        for (int i= 0; i < n; i++)
          cpos += tetrahedra[vtxList[i]].centre();
        cpos *= 1.0/(real)n;

        const vec3 &posA = tetrahedra[vtxList[0]].centre() - cpos;
        int iv = 0;
        for (iv = 1; iv < n; iv++)
        {
          const vec3 &posB = tetrahedra[vtxList[iv]].centre() - cpos;
          const real q = (posA%posB).norm2();
          if (q*q > eps*eps*eps*eps*posA.norm2()*posB.norm2())
            break;
        }
        if (iv == n) return false;
        const vec3 &posB = tetrahedra[vtxList[iv]].centre() - cpos;
#if 0
        fprintf(stderr, "c= %g %g %g | %g\n", cpos.x, cpos.y, cpos.z, cpos.abs());
        fprintf(stderr, "a= %g %g %g | %g\n", posA.x, posA.y, posA.z, posA.abs());
        fprintf(stderr, "b= %g %g %g | %g\n", posB.x, posB.y, posB.z, posB.abs());
#endif
        const vec3 &unitA = posA * (1.0/posA.abs());
        const vec3 &unitB = posB * (1.0/posB.abs());
        vec3  unitN  = unitA%unitB;
#if 0
        fprintf(stderr, "unitA= %g %g %g | %g\n", unitA.x, unitA.y, unitA.z, unitA.abs());
        fprintf(stderr, "unitB= %g %g %g | %g\n", unitB.x, unitB.y, unitB.z, unitB.abs());
        fprintf(stderr, "unitN= %g %g %g | %g\n", unitN.x, unitN.y, unitN.z, unitN.abs());
#endif
        if (unitN.abs() == 0.0) return false;
        unitN *= 1.0/unitN.abs();

        angle_vec_pair.clear();
        angle_vec_pair.push_back(std::make_pair(0.0, posA));
#if 0
        fprintf(stderr, " -- nvtx= %d\n", n);
#endif
        for (int j = 1; j < n; j++)
        {
          const vec3 &jpos = tetrahedra[vtxList[j]].centre() - cpos;
          const real x =  posA * jpos;
          const real y = (posA % jpos) * unitN;
          const real dot = jpos * unitN;
#if 1
          if (!(__abs(dot) < jpos.abs()*1.0e-10))
            return false;
#endif
          assert(__abs(dot) < jpos.abs()*1.0e-10);
          angle_vec_pair.push_back(std::make_pair(std::atan2(y, x), jpos));
#if 0
          fprintf(stderr, " -- i= %d: pos= %g %g %g  phi= %g\n", 
              j, jpos.x, jpos.y, jpos.z, std::atan2(y,x));
#endif
        }

        assert((int)angle_vec_pair.size() == n);
        std::sort(angle_vec_pair.begin(), angle_vec_pair.end(), cmp_data<real, vec3>());
        angle_vec_pair.push_back(angle_vec_pair[0]);
#if 0
        for (int i = 0; i < n+1; i++)
        {
          const vec3 &jpos = angle_vec_pair[i].second;
          fprintf(stderr, " -- i= %d: pos= %g %g %g  phi= %g\n", 
              i, jpos.x, jpos.y, jpos.z, angle_vec_pair[i].first);
        }
#endif

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

#if 0
        fprintf(stderr, " vol= %g faceA= %g  n= %g %g %g | %g\n",
            vol, area, unitN.x, unitN.y, unitN.z, unitN.abs());
#endif
        face = Face(unitN, area);
        return true;
      }

    };




}

#endif /* __VOROCELL_DEGERATE_H__ */
