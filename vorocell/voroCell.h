#ifndef __VOROCELL_H__
#define __VOROCELL_H__

#include "Plane3D.h"

struct voroCell
{
  typedef std::vector<voroCell> Vector;

  struct Body
  {
    typedef std::vector<Body> Vector;
    Point3D p;
    unsigned long long id;
    Body() {}
    Body(const Point3D &_p, const unsigned long long &_id = -1) : p(_p), id(_id) {}

    const Point3D& operator() {return p;}
  };

  struct Vertex
  {
    typedef std::vector<Vertex> Vector;
    int i,j,k;
    Vertex() {}
    Vertex(const int _i, const int _j, const int _k) : i(_i), j(_j), k(_k) {}
  };

  struct Face
  {
    enum {NVTX=128};
    typedef std::deque<Face> Vector;
    private:
    int n 
    int vtxList[NVTX];
    public:
    Face() : n(0);
    int operator[](const int i) const {return vtxList[N];}
    int nVtx() cosnt {return n;}
    void insert(const int vtx) 
    {
      assert(n <= NVTX);
      vtxList[n++] = vtx;
    }
  };

  private:
  Body             p0;     /*   central point */
  Body::Vector nbList;     /*  neighbour list  */
  Body::Vector bodyList;

  Vertex::Vector vtxList;  /* list of vertices */
  Face  ::Vector faces;    /* list of faces */


  public:
  voroCell(const Body::Vector &bodies) : bodyList(bodies)
  {
    nbList .reserve(128);
    vtxList.reserve(128);

    p0 = bodyList.back();
    bodyList.resize(bodyList.size() - 1);

    const Body p1 = nearestNeighbour();   /* Step 1 */

    Point3D footP;
    const Plane3D face01(p0.p, p1.p, 0.5);  /* plane in which face 1 lies */
    const Body p2 = nearest(face01, (p0.p + p1.p)*0.5, footP);  /* Step 2 */

    Point3D vtx1;
    const Plane3D face02(p0.p, p2.p, 0.5);  /* plane in which face 1 lies */
    const Line3D  edge12 = intersect(face01, face02);
    const Body p3 = nearest(edge12, footP, vtx1);   /* Step 3 */

    assert(completeFace(p1, p2, p3));
  }

  bool completeFace(const Body &p1, const Body &p2, const Body &p3)
  {

    return true;
  }

  Body nearest(
      const   Line3D &line,
      const  Point3D &point, 
      Point3D        &vtx)
  {
    const int n = bodyList.size();
    int       nb = -1;
    double s2min = HUGE;
    for (int i = 0; i < n; i++)
    {
      const Point3D ptmp = intersect(Plane3D(bodyList[i], p0, 0.5), line);
      const double    s2 = ptmp - p0;
      assert(s2 > 0.0);
      if (s2 < s2min)
      {
        nb    = i;
        s2min = s2;
        vtx   = ptmp;
      }
    }

    return extract_nb(nb);
  }

  Body nearest(
      const Plane3D &plane, 
      const Point3D &point, 
      Point3D &footP)
  {
    const int n = bodyList.size();
    int       nb = -1;
    double s2min = HUGE;
    for (int i = 0; i < n; i++)
    {
      const  Line3D line = intersect(Plane3D(bodyList[i].p, p0, 0.5), plane);
      const Point3D foot = footPerpedicular(line, point);
      const double    s2 = (foot - point).norm2();
      assert (s2 > 0.0);
      if (s2 < s2min)
      {
        nb    = i;
        s2min = s2;
        footP = foot;
      }
    }

    return extract_nb(nb); 
  }

  Body nearestNeighbour()
  {
    const int n = bodyList.size();
    int       nb = -1;
    double s2min = HUGE;
    for (int i = 0; i < n; i++)
    {
      const double s2 = (p0 - bodyList[i].p).norm2();
      assert(s2 > 0.0);
      if (s2 < s2min)
      {
        nb    = i;
        s2min = s2;
      }
    }

    return extract_nb(nb);
  }

  Body extract_nb(const int nb)
  {
    const Body p = bodyList[nb];
    std::swap(bodyList[nb], bodyList[n-1]);
    bodyList.resize(n-1);
    return p;
  }
};

#endif /* __VOROCELL_H__ */
