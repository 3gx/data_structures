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

  const int NFACEMAX = 32;
  DirectPolyhedron<NFACEMAX> d;
  Voronoi::Cell<NFACEMAX> cell;
  Voronoi::Site::Vector list;
  list.reserve(NFACEMAX);

  const double t0 = get_wtime();
  const int cnt = 10;
  for (int icnt = 0; icnt < cnt; icnt++)
  {
    const real f = 0.05;
    const vec3 pos(
        0.5 + (1.0 - 2.0*drand48()) * f,
        0.5 + (1.0 - 2.0*drand48()) * f,
        0.5 + (1.0 - 2.0*drand48()) * f);
    d.clear();
    int nsuccess = 0;
    for (int i= 0; i < np; i++)
      nsuccess += d.push(ptcl[i].pos - pos, i);
    fprintf(stderr, " --- \n");
    fprintf(stderr, " icnt= %d: nsuccess= %d  nface= %d\n", icnt, nsuccess, d.nface());

#if 1
    list.clear();
    for (int i = 0; i < d.nface(); i++)
      list.push_back(Voronoi::Site(ptcl[d[i]].pos-pos, d[i]));
    cell.build(list);
    fprintf(stderr, "  nf= %d vol= %g\n", cell.nb(), cell.volume());
#endif
    
#if 0
    d.clear();
    for (int i= 0; i < np; i++)
      ptcl[i].r2 = (ptcl[i].pos - pos).norm2();
    std::sort(ptcl.begin(), ptcl.end(), Particle());
    nsuccess = 0;
#if 1
    for (int i= 0; i < np; i++)
      nsuccess += d.push(ptcl[i].pos - pos);
#endif
    fprintf(stderr, "  sorted: nsuccess= %d  nface= %d\n", nsuccess, d.nface());
#endif
  }
  const double dt = get_wtime() - t0;

  fprintf(stderr, " done in %g sec : %g ptcl/sec \n", dt/cnt, np*cnt/dt);


  return 0;
}
