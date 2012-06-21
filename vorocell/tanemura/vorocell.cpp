#include "vorocell.h"
#include "mytimer.h"
#include <iostream>

typedef Voronoi::real real;
typedef Voronoi::vec3 vec3;

struct cmp_dist
{
  bool operator()(const std::pair<float, int> &lhs, const std::pair<float, int> &rhs) const
  {
    return lhs.first < rhs.first;
  }
};

int main(int argc, char * argv[])
{
  double eps;
  int ns, nrg;
  double lx, ly, lz;
  int np;

  int idum;
  double fdum;
  std::string sdum;

  std::cin >> eps >> sdum;
  std::cin >> ns >> nrg >> sdum >> sdum;
  std::cin >> fdum >> sdum;
  std::cin >> idum >> sdum >> idum;
  std::cin >> idum >> sdum;
  std::cin >> idum >> sdum;
  std::cin >> idum >> sdum;
  std::cin >> idum >> sdum;
  std::cin >> idum >> np >> sdum >> sdum;
  std::cin >> lx >> ly >> lz >> sdum >> sdum >> sdum;

  fprintf(stderr, " eps= %g \n", eps);
  fprintf(stderr, " ns= %d  nrg= %d\n", ns, nrg);
  fprintf(stderr, " np= %d\n", np);
  fprintf(stderr, " l= %g %g %g \n", lx, ly, lz);

  Voronoi::Site::Vector sites(np);
  vec3 min(+HUGE), max(-HUGE);
  for (int i = 0; i < np; i++)
  {
    std::cin >> idum >> 
      sites[i].pos.x >> 
      sites[i].pos.y >> 
      sites[i].pos.z;
    sites[i].idx = i;
    min = mineach(min, sites[i].pos);
    max = maxeach(max, sites[i].pos);
  }
  for (int i = 0; i < np; i++)
    sites[i].pos -= min;
  min -= min;
  max -= min;
  fprintf(stderr, " min= %g %g %g \n", min.x, min.y, min.z);
  fprintf(stderr, " max= %g %g %g \n", max.x, max.y, max.z);
  for (int i = 0; i < np; i++)
  {
    assert(sites[i].pos.x < lx);
    assert(sites[i].pos.y < ly);
    assert(sites[i].pos.z < lz);
  }

  Voronoi::Site::Vector sitesP;
  sitesP.reserve(4*np);
#if 1 /* periodic */
  const real  f = 0.5;
  const real dx = f * lx;
  const real dy = f * ly;
  const real dz = f * lz;
  for (int i = 0; i < np; i++)
  {
    Voronoi::Site &s = sites[i];
    sitesP.push_back(s);

    if      (s.pos.x      < dx) sitesP.push_back(Voronoi::Site(vec3(s.pos.x+lx, s.pos.y, s.pos.z), -1-s.idx));
    else if (lx - s.pos.x < dx) sitesP.push_back(Voronoi::Site(vec3(s.pos.x-lx, s.pos.y, s.pos.z), -1-s.idx));

    if      (s.pos.y      < dy) sitesP.push_back(Voronoi::Site(vec3(s.pos.x, s.pos.y+ly, s.pos.z), -1-s.idx));
    else if (ly - s.pos.y < dy) sitesP.push_back(Voronoi::Site(vec3(s.pos.x, s.pos.y-ly, s.pos.z), -1-s.idx));

    if      (s.pos.z      < dz) sitesP.push_back(Voronoi::Site(vec3(s.pos.x, s.pos.y, s.pos.z+lz), -1-s.idx));
    else if (lz - s.pos.z < dz) sitesP.push_back(Voronoi::Site(vec3(s.pos.x, s.pos.y, s.pos.z-lz), -1-s.idx));
  }
#else /* reflecting */
#endif
  fprintf(stderr, " np= %d  Pnp= %d\n", (int)sites.size(), (int)sitesP.size());
  assert(!sitesP.empty());

  double dt_search = 0.0;
  double dt_voro   = 0.0;
  const double tbeg = get_wtime();
  Voronoi::Site::Vector list;
  list.reserve(ns);

  std::vector< std::pair<float, int> > dist(sitesP.size());
  for (int i = 0; i < np; i++)
  {
    const Voronoi::Site &s = sites[i];

    double t0 = get_wtime();
    for (int i = 0; i < (const int)sitesP.size(); i++)
      dist[i] = std::make_pair((sitesP[i].pos - s.pos).norm2(), i);
    std::nth_element(dist.begin(), dist.begin() + ns, dist.end(), cmp_dist());
    list.clear();
    for (int i = 0; i < ns+1; i++)
      if (dist[i].second > 0.0f)
        list.push_back(Voronoi::Site(sitesP[dist[i].second].pos - s.pos, sitesP[dist[i].second].idx));
    fprintf(stderr, " list.size()=  %d  ns= %d \n", (int)list.size(), ns);
    assert((int)list.size() == ns);
    double t1 = get_wtime();
    dt_search += t1 - t0;
  }
  const double tend = get_wtime();

  fprintf(stderr , " dt_search=  %g sec \n", dt_search);
  fprintf(stderr , " dt_voro  =  %g sec \n", dt_voro);
  fprintf(stderr , " dt_total =  %g sec [ sum= %g ]\n", tend - tbeg,
      dt_search + dt_voro);



  return 0;
};
