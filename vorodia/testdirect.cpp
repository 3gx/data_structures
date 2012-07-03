#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include "DirectPoly.h"
#include "vorocell.h"
#include "mytimer.h"

typedef double real;
typedef vector3<real> vec3;

struct Particle
{
  typedef std::vector<Particle> Vector;
  vec3 pos;
  int  idx;
  float r2;
  Particle() {}
  Particle(const vec3 &p, const int i) : pos(p), idx(i) {}
  bool operator() (const Particle &p1 ,const Particle &p2) const
  {
    return p1.r2 < p2.r2;
  }
};

int main (int argc, char * argv [])
{
  assert(argc > 1);
  const int np = atoi(argv[1]);
  fprintf(stderr, " np= %d\n", np);


  Particle::Vector ptcl(np);
  for (int i = 0; i < np; i++)
    ptcl[i] = Particle(vec3(drand48(), drand48(), drand48()), i);

  const int NFACEMAX = 128;
  DirectPolyhedron<NFACEMAX> d;
  Voronoi::Cell<NFACEMAX> cell;
  Voronoi::Site::Vector list;
  list.reserve(NFACEMAX);

  const double t0 = get_wtime();
  const int cnt = 15;
  std::vector<int> nblist;
  std::swap(ptcl[0], ptcl[42255]);
  ptcl[0].idx = 0;
  ptcl[42255].idx = 42255;
  for (int icnt = 0; icnt < cnt; icnt++)
  {
    const real f = 0.05;
    const vec3 pos(
        0.5 + (1.0 - 2.0*drand48()) * f,
        0.5 + (1.0 - 2.0*drand48()) * f,
        0.5 + (1.0 - 2.0*drand48()) * f);

#if 0
    if (icnt < 1) continue;
    if (icnt > 1) break;
#endif
    d.clear();
    int nsuccess = 0;
    for (int i= 0; i < np; i++)
      nsuccess += d.push((ptcl[i].pos - pos)*0.5, i, 1.25);
#if 0
    for (int i= 0; i < np; i++)
      nsuccess += d.push((ptcl[i].pos - pos)*0.5, i);
#endif
    fprintf(stderr, " --- \n");
    fprintf(stderr, " icnt= %d: nsuccess= %d  nface= %d\n", icnt, nsuccess, d.nface());

#if 1
    list.clear();
    for (int i = 0; i < d.nface(); i++)
      list.push_back(Voronoi::Site(ptcl[d[i]].pos-pos, d[i]));
#if 0
    list.push_back(Voronoi::Site(vec3(-2.0*pos.x,0.0,0.0), -1));
    list.push_back(Voronoi::Site(vec3(2.0*(1.0-pos.x),0.0,0.0), -2));
    list.push_back(Voronoi::Site(vec3(0.0, -2.0*pos.y,0.0), -3));
    list.push_back(Voronoi::Site(vec3(0.0, 2.0*(1.0-pos.y),0.0), -4));
    list.push_back(Voronoi::Site(vec3(0.0, 0.0, -2.0*pos.z), -5));
    list.push_back(Voronoi::Site(vec3(0.0, 0.0, 2.0*(1.0-pos.z)), -6));
#endif
    cell.build(list);
    fprintf(stderr, "  nf= %d vol= %g\n", cell.nb(), cell.volume());
#endif
    
#if 0
    nblist.clear();
    for (int i = 0; i < d.nface(); i++)
      nblist.push_back(ptcl[d[i]].idx);

    d.clear();
    Particle::Vector ptcl1(ptcl);
    for (int i= 0; i < np; i++)
      ptcl1[i].r2 = (ptcl1[i].pos - pos).norm2()*0.25;
    std::sort(ptcl1.begin(), ptcl1.end(), Particle());
    nsuccess = 0;
#if 1
    for (int i= 0; i < np; i++)
      nsuccess += d.push((ptcl1[i].pos - pos)*0.5, i);
#endif
    fprintf(stderr, "  sorted: nsuccess= %d  nface= %d\n", nsuccess, d.nface());

#if 0
    std::sort(nblist.begin(), nblist.end());
    for (int i = 0; i < (int)nblist.size(); i++)
      fprintf(stderr, " %d", nblist[i]);
    fprintf(stderr, " \n");

    nblist.clear();
    for (int i = 0; i < d.nface(); i++)
      nblist.push_back(ptcl1[d[i]].idx);
    std::sort(nblist.begin(), nblist.end());
    for (int i = 0; i < d.nface(); i++)
      fprintf(stderr, " %d", nblist[i]);
    fprintf(stderr, " \n");
#endif
#endif
  }
  const double dt = get_wtime() - t0;

  fprintf(stderr, " done in %g sec : %g ptcl/sec \n", dt/cnt, np*cnt/dt);


  return 0;
}
