#include "mytimer.h"
double dt_00, dt_10, dt_20, dt_30, dt_40, dt_44, dt_50, dt_60, dt_70;
double dtA;
unsigned long long flop = 0, myNMX = 0;
#if 0
#include "vorocell.h"
#else
#include "vorocell_new.h"
#endif
#include "vorocell_degenerate.h"
#include <iostream>

typedef Voronoi::real real;
typedef Voronoi::vec3 vec3;

int main(int argc, char * argv[])
{
  dt_00=dt_10=dt_20=dt_30=dt_40=dt_44=dt_50=dt_60=dt_70=0.0;
  dtA = 0.0;
  flop = 0;
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

#if 0
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
#else
  const int N = 10;
  Voronoi::Site::Vector sites;
  sites.reserve(np);
  vec3 min(+HUGE), max(-HUGE);
  {
    const real dx = lx / N;
    const real dy = lx / N;
    const real dz = lx / N;
    for (int k = 0; k < N; k++)
      for (int j = 0; j < N; j++)
        for (int i = 0; i < N; i++)
        {
          sites.push_back(Voronoi::Site(vec3(dx*(i+0.5), dy*(j+0.5), dz*(k+0.5)), sites.size()));
          vec3 &pos = sites.back().pos;
          const real fac = 0.0e-12;
          pos.x += fac*(1.0 - 2.0*drand48())*dx;
          pos.y += fac*(1.0 - 2.0*drand48())*dy;
          pos.z += fac*(1.0 - 2.0*drand48())*dz;
          min = mineach(min, sites.back().pos);
          max = maxeach(max, sites.back().pos);
        }
  }
  std::random_shuffle(sites.begin(), sites.end());
  assert((int)sites.size() == np);

#endif
  fprintf(stderr, " min= %g %g %g \n", min.x, min.y, min.z);
  fprintf(stderr, " max= %g %g %g \n", max.x, max.y, max.z);
  assert(min.x > 0.0);
  assert(min.y > 0.0);
  assert(min.z > 0.0);
  assert(max.x < lx);
  assert(max.y < ly);
  assert(max.z < lz);

  Voronoi::Site::Vector sitesP;
  sitesP.reserve(8*np);
#if 1 /* periodic */
  const real  f = 0.5;
  assert(f <= 0.5);
  const real dx = (0.5 - f) * lx;
  const real dy = (0.5 - f) * ly;
  const real dz = (0.5 - f) * lz;
  const vec3 cpos(0.5*lx, 0.5*ly, 0.5*lz);
  for (int i = 0; i < np; i++)
  {
    const Voronoi::Site &s0 = sites[i];
    sitesP.push_back(s0);
    //   fprintf(stderr, "%g %g %g \n", s0.pos.x, s0.pos.y, s0.pos.z);
    for (int oct = 1; oct < 8; oct++)
    {
      Voronoi::Site s = s0;
      int ioct = 0;
      if ((oct&1) && __abs(s.pos.x-cpos.x) > dx) {s.pos.x += lx * (s.pos.x <= cpos.x ? +1.0 : -1.0); ioct += 1;}
      if ((oct&2) && __abs(s.pos.y-cpos.y) > dy) {s.pos.y += ly * (s.pos.y <= cpos.y ? +1.0 : -1.0); ioct += 2;}
      if ((oct&4) && __abs(s.pos.z-cpos.z) > dz) {s.pos.z += lz * (s.pos.z <= cpos.z ? +1.0 : -1.0); ioct += 4;}
      if (oct == ioct)
      {
        s.idx = -1-s.idx;
        //        fprintf(stderr, "%g %g %g \n", s.pos.x, s.pos.y, s.pos.z);
        sitesP.push_back(s);
      }
    }
  }
#else  /* reflecting */
#define REFLECTING
  sitesP = sites;
#endif

  fprintf(stderr, " np= %d  Pnp= %d\n", (int)sites.size(), (int)sitesP.size());
  assert(!sitesP.empty());

  double dt_search = 0.0;
  double dt_sort   = 0.0;
  double dt_voro   = 0.0;
  double nface     = 0.0;
  const double tbeg = get_wtime();
  double volume    = 0.0;
  int nfailed = 0;

 
 const int CNT = 10; 
  for (int cnt = 0; cnt < CNT; cnt++)
  {
    volume    = 0.0;
    nfailed = 0;

#pragma omp parallel reduction(+:dt_search, dt_sort, dt_voro, nface, volume, nfailed)
    {
      Voronoi          ::Cell<128> cell;
      VoronoiDegenerate::Cell<128> cell_degenerate;
      std::vector< std::pair<real, int> > dist(sitesP.size());

      Voronoi::Site::Vector list;
      list.reserve(ns);


#pragma omp for
      for (int i = 0; i < np; i++)
      {
        const Voronoi::Site &s = sites[i];
#if 0
        fprintf(stderr, "i= %d  %g %g %g \n", i, s.pos.x, s.pos.y, s.pos.z);
#endif

        double t0 = get_wtime();
        for (int j = 0; j < (const int)sitesP.size(); j++)
          dist[j] = std::make_pair((sitesP[j].pos - s.pos).norm2(), j);
        
        list.clear();
#ifdef REFLECTING
        std::nth_element(dist.begin(), dist.begin() + ns-3, dist.end(), cmp_data<real, int>());
        for (int j = 0; j < ns+1; j++)
          if (dist[j].first > 0.0)
            list.push_back(Voronoi::Site(sitesP[dist[j].second].pos - s.pos, sitesP[dist[j].second].idx));
        assert((int)list.size() == ns);
#if 1
        list.resize(ns-3);
        if (s.pos.x < 0.5*lx)  list.push_back(Voronoi::Site(vec3(-2.0*      s.pos.x,  0.0, 0.0), -1-s.idx));
        else                   list.push_back(Voronoi::Site(vec3( 2.0*(lx - s.pos.x), 0.0, 0.0), -1-s.idx));
        if (s.pos.y < 0.5*ly)  list.push_back(Voronoi::Site(vec3(0.0, -2.0*      s.pos.y,  0.0), -1-s.idx));
        else                   list.push_back(Voronoi::Site(vec3(0.0,  2.0*(ly - s.pos.y), 0.0), -1-s.idx));
        if (s.pos.z < 0.5*lz)  list.push_back(Voronoi::Site(vec3(0.0, 0.0, -2.0*      s.pos.z ), -1-s.idx));
        else                   list.push_back(Voronoi::Site(vec3(0.0, 0.0,  2.0*(lz - s.pos.z)), -1-s.idx));
#else
        list.resize(ns-6);
        list.push_back(Voronoi::Site(vec3(-2.0*      s.pos.x,  0.0, 0.0), -1-s.idx));
        list.push_back(Voronoi::Site(vec3( 2.0*(lx - s.pos.x), 0.0, 0.0), -1-s.idx));
        list.push_back(Voronoi::Site(vec3(0.0, -2.0*      s.pos.y,  0.0), -1-s.idx));
        list.push_back(Voronoi::Site(vec3(0.0,  2.0*(ly - s.pos.y), 0.0), -1-s.idx));
        list.push_back(Voronoi::Site(vec3(0.0, 0.0, -2.0*      s.pos.z ), -1-s.idx));
        list.push_back(Voronoi::Site(vec3(0.0, 0.0,  2.0*(lz - s.pos.z)), -1-s.idx));
#endif
#else
        std::nth_element(dist.begin(), dist.begin() + ns, dist.end(), cmp_data<real, int>());
        for (int j = 0; j < ns+1; j++)
          if (dist[j].first > 0.0f)
            list.push_back(Voronoi::Site(sitesP[dist[j].second].pos - s.pos, sitesP[dist[j].second].idx));
#endif
        assert((int)list.size() == ns);
        for (int j = 0; j < ns; j++)
          list[j].r = list[j].pos.norm2();
        double t1 = get_wtime();
        dt_search += t1 - t0;
        t0 = t1;

//        std::sort(list.begin(), list.end(), Voronoi::Site());
//        std::random_shuffle(list.begin(), list.end());
        t1 = get_wtime();
        dt_sort += t1 - t0;
        t0 = t1;
        if (!cell.build(list))
        {
          fprintf(stderr, "i =%d\n", i);
          assert(0);
#if 1
          assert(cell_degenerate.build(list));
          volume += cell_degenerate.volume();
          nface  += cell_degenerate.nb();
#endif
          nfailed++;
        }
        else
        {
          volume += cell.volume();
          nface  += cell.nb();
        }
        t1 = get_wtime();
        dt_voro += t1 - t0;
      }
    }
  }
  const double tend = get_wtime();

  fprintf(stderr , " dt_search=  %g sec \n", dt_search);
  fprintf(stderr , " dt_sort  =  %g sec \n", dt_sort);
  fprintf(stderr , " dt_voro  =  %g sec \n", dt_voro);
  fprintf(stderr,  "   dt_00=  %g \n" ,dt_00);
  fprintf(stderr,  "   dt_10=  %g \n" ,dt_10);
  fprintf(stderr,  "   dt_20=  %g GFLOP/s: %g  flop= %llu M\n" ,dt_20, flop/dt_20/1e9, flop/1000000 );
  fprintf(stderr,  "   dt_30=  %g GFLOP/s: %g rate= %g cell/sec\n" ,dt_30, flop/dt_30/1e9, np*CNT/dt_30);
  fprintf(stderr,  "   dt_40=  %g \n" ,dt_40);
  fprintf(stderr,  "   dt_44=  %g \n" ,dt_44);
  fprintf(stderr,  "   dt_50=  %g \n" ,dt_50);
  fprintf(stderr,  "   dt_60=  %g \n" ,dt_60);
  fprintf(stderr,  "   dt_70=  %g \n" ,dt_70);
  fprintf(stderr,  "   dtA=  %g \n" ,dtA);
  //  fprintf(stderr,  "    sum= %g\n", dt_00+dt_10+dt_20+dt_30+dt_40+dt_50+dt_60);
  fprintf(stderr , " dt_total =  %g sec [ sum= %g ]\n", tend - tbeg,
      dt_search + dt_voro);
  fprintf(stderr, "   nface= %g \n", nface/np);

  fprintf(stderr, " volume= %g   exact= %g  diff= %g \n",
      volume, lx*ly*lz, (volume-lx*ly*lz)/(lx*ly*lz));

  fprintf(stderr,  "ncell= %d  nfailed= %d : %g \n",      np, nfailed, (real)nfailed/np);

  fprintf(stderr, "myNMX/pp= %g\n", myNMX*1.0/(CNT*np));

  return 0;
};
