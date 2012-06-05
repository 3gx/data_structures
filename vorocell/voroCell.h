#ifndef __VOROCELL_H__
#define __VOROCELL_H__

#include "Plane3D.h"
#include <map>

#define SMALLDIST (1.0e-10)

template<const int M, const int N>
class boolMatrix
{
  private:
  unsigned int matrix[M][N>>3];

  public:

  boolMatrix() 
  {
    assert((N & 7) == 0);
    for (int j = 0; j < M; j++)
      for (int i = 0; i < (N>>3); i++)
        matrix[j][i] = 0;
  }

  unsigned int bit(const int ch) const {return 1 << ch;}

  bool operator()(const int row, const int col) const
  {
    return matrix[row][col>>3] & bit(col&7);
  }

  void   set(const int row, const int col) {matrix[row][col>>3] |=  bit(col&7);}
  void unset(const int row, const int col) {matrix[row][col>>3] &= ~bit(col&7);}
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

namespace Voro
{

  struct Site
  {
    typedef std::vector<Site> Vector;
    Point3D p;
    int id;
    Site() {}
    Site(const Point3D &_p, int _id) : p(_p), id(_id) {}

    operator  Point3D() const {return  p;}

    bool operator==(const Site &s) const {return s.id == id;}
    bool operator!=(const Site &s) const {return s.id != id;}
  };

  struct Vertex
  {
    typedef std::vector<Vertex> Vector;
    int s[3];  /* sites */
    Point3D p;
    Vertex() {}
    Vertex(const Point3D &_p) : p(_p) {}
  };

  struct Edge
  {
    typedef std::vector<Edge> Vector;
    int beg, end;
    Edge() {}
    Edge(const int _beg, const int _end) : beg(_beg), end(_end) {};
  };

  struct Face
  {
    enum {NVTX=1024};
//    typedef std::deque<Face> Vector;

    private:
    Site    site[3];       /*  generating site   */
    int     nvtx;          /* number of vertices */
    Point3D vtxList[NVTX]; /* list of vertices */

    public:
    Face() {}
    Face(const Site &s1, const Site &s2, const Site &s3) : nvtx(0)
    {
      site[0] = s1;
      site[1] = s2;
      site[2] = s3;
    }

    void clear() 
    {
      nvtx = 0;
    }

    const Site& getSite(const int i = 0) const {return site[i];}

    Point3D operator[](const int i) const {return vtxList[i];}
    int nVtx() const {return nvtx;}
    void insert(const Point3D &vtx) 
    {
      assert(nvtx <= NVTX);
      vtxList[nvtx++] = vtx;
    }

  };

  struct Cell
  {
    enum {NFACEMAX=128, NSITEMAX=1024};

    typedef Array<Site, NSITEMAX> siteArray;
    typedef Array<bool, NSITEMAX> boolArray;
    typedef Array<Face, NFACEMAX> faceArray;

    private:
    siteArray _siteList;     /*  input site list */
    boolArray _siteFlag;
    faceArray _faceList;     /* list of    faces */
    double    _distVtx_s2;   /* distance to the most distant vertex */

    public:
    int nFace() const {return _faceList.size();}
    double r2() const {return _distVtx_s2;}
    Cell() {}
    void build(const Site &site0, const Site::Vector &sites) 
    {
      _siteList = sites;
      _siteFlag.resize(sites.size());
      _faceList.clear();

      for (int i = 0; i < _siteList.size(); i++)
      {
        _siteList[i].id = i;
        _siteFlag[i]    = 0;
        _siteList[i].p -= site0;
      }
      

      const int nnb = nearestNeighbour();
      const Site s1 = _siteList.erase(nnb); /* Step 1 */

      Point3D footP;
      const Plane3D face01(s1, 0.5);  /* plane in which face 1 lies */
      const Site s2 = _siteList.erase(nearest(face01, s1.p*0.5, footP)); /* Step 2*/

      const Plane3D face02(s2, 0.5);  /* plane in which face 1 lies */
      const Line3D  edge12 = intersect(face01, face02);
      const Site s3 = _siteList.erase(nearest(face02, edge12, footP));   /* Step 3 */

      _siteList.push_back(s1);
      _siteList.push_back(s2);
      _siteList.push_back(s3);

      _faceList.push_back(Face(s1,s2,s3));
      _siteFlag[s1.id] = true;


      /* now build faces */

      int nface = 0;
//      while (nface < _faceList.size())
      {
#if 0
        fprintf(stderr, " --------------- \n");
        fprintf(stderr, "nface= %d  nfaces= %d \n", nface, _faceList.size());
#endif
        buildFace(_siteList, _faceList[nface++]);
      }

      for (int i = 0; i < _siteFlag.size(); i++)
        _siteFlag[i] = 0;
    }

    bool buildFace(
        siteArray  sites,
        Face      &face)
    {

      const Site &s1 = face.getSite(0);
      const Site &s2 = face.getSite(1);
      const Site &s3 = face.getSite(2);
      
      sites.erase(s1);
      sites.erase(s3);

      const Plane3D face01(s1, 0.5);
      const Plane3D face02(s2, 0.5); 
      const Plane3D face03(s3, 0.5);  
      
      Point3D newVtx;
      Point3D    vtx = intersect(face03, intersect(face01, face02));
      Plane3D  faceA = face02;
      Plane3D  faceB = face03;
      Site        sA = s2;
      Site        sB = s3;

      store_distVtx(vtx);
      face.insert(vtx);

      if (!_siteFlag[sA.id])
      {
        _faceList.push_back(Face(sA, sB, s1));
        _siteFlag[sA.id] = true;
      }
      if (!_siteFlag[sB.id])
      {
        _faceList.push_back(Face(sB, s1, sA));
        _siteFlag[sB.id] = true;
      }

      int cnt = 0;
      while(1)
      {
        cnt++;
        const int nb = nearest(sites, faceA, intersect(face01, faceB), vtx, newVtx);
        const Site sBnew = sites.erase(nb);
        if (!_siteFlag[sBnew.id]) 
        {
          _faceList.push_back(Face(sBnew, s1, sB));
          _siteFlag[sBnew.id] = true;
        }

        store_distVtx(vtx);
        face.insert(vtx);
        if (sBnew == s2)
          break;

        sA  = sB;
        sB  = sBnew;
        faceA = faceB;
        faceB = Plane3D(sB, 0.5);

        vtx = newVtx;
      }

      return true;
    }

    /* keeps track of the farthers vertex */
    bool store_distVtx(const Point3D &vtx)
    {
      const double s2 = vtx.norm2();
      if (s2 <= _distVtx_s2) return false;

      _distVtx_s2 = s2;
      return true;
    }

    /* find the nearest site in "sites" behind the "plane" whose face
     * intersection with the "line" gives the "vertex" closes to the "point"
     */
    int nearest(
        const siteArray &sites,
        const Plane3D   &plane,
        const Line3D    &line,
        const Point3D   &point,
        Point3D         &vertex) const
    {
      const int n  = sites.size();
      int       nb = -1;
      double s2min = HUGE;
      for (int i = 0; i < n; i++)
      {
        const Point3D ptmp = intersect(Plane3D(sites[i], 0.5), line);
        if (distance(plane, ptmp) >= -SMALLDIST) continue;
        const double    s2 = (ptmp - point).norm2();
        assert(s2 > 0.0);
        if (s2 < s2min)
        {
          nb     = i;
          s2min  = s2;
          vertex = ptmp;
        }
      }

      assert(nb >= 0);
      return nb;
    }


    /* Find the nearest site behind the "plane" whose face itersection with the
     * "line" is the closest to the "point"
     */
    int nearest(
        const  Plane3D &plane,
        const   Line3D &line,
        const  Point3D &point)
    {
      const int n = _siteList.size();
      int       nb = -1;
      double s2min = HUGE;
      for (int i = 0; i < n; i++)
      {
        const Point3D ptmp = intersect(Plane3D(_siteList[i], 0.5), line);
        if (distance(plane, ptmp) > SMALLDIST) continue;
        const double    s2 = (ptmp - point).norm2();
        assert(s2 > 0.0);
        if (s2 < s2min)
        {
          nb    = i;
          s2min = s2;
        }
      }

      assert(nb >= 0);
      return nb;
    }


    /* Find the nearest site whose face intersection with the "plane" is  
     * closest line to the "point"  
     */ 
    int nearest(
        const Plane3D &plane, 
        const Point3D &point, 
        Point3D &footP) const
    {
      const int n  = _siteList.size();
      int       nb = -1;
      double s2min = HUGE;
      for (int i = 0; i < n; i++)
      {
        const  Line3D line = intersect(Plane3D(_siteList[i].p, 0.5), plane);
        const Point3D foot = perpendicularFoot(line, point);
        const double    s2 = (foot - point).norm2();
        assert (s2 > 0.0);
        if (s2 < s2min)
        {
          nb    = i;
          s2min = s2;
          footP = foot;
        }
      }

      assert(nb >= 0);
      return nb;
    }

    /* Find the nearest neighbour site in _siteList */
    int nearestNeighbour() const
    {
      const int n  = _siteList.size();
      int       nnb = -1;
      double s2min = HUGE;
      for (int i = 0; i < n; i++)
      {
        const double s2 = _siteList[i].p.norm2();
        assert(s2 > 0.0);
        if (s2 < s2min)
        {
          nnb   = i;
          s2min = s2;
        }
      }

      assert(nnb >= 0);
      return nnb;
    }
  };

}

#endif /* __VOROCELL_H__ */
