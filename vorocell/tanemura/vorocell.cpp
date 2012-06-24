#include "mytimer.h"
double dt_00, dt_10, dt_20, dt_30, dt_40, dt_44, dt_50, dt_60, dt_70;
double dtA;
long long incomplT = 0;
long long incompl  = 0;
#include "vorocell.h"
#include "vorocell_degenerate.h"
#include <iostream>

typedef Voronoi::real real;
typedef Voronoi::vec3 vec3;

int main(int argc, char * argv[])
{
  dt_00=dt_10=dt_20=dt_30=dt_40=dt_44=dt_50=dt_60=dt_70=0.0;
  dtA = 0.0;
  double eps;
  int ns, nrg;
  double lx, ly, lz;
  int np;

  int idum;
  double fdum;
  std::string sdum;

  incomplT = incompl = 0;

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

#if 1
  Voronoi::Site::Vector sites(np);
  vec3 min(+HUGE), max(-HUGE);
  for (int i = 0; i < np; i++)
  {
    std::cin >> idum >> 
      sites[i].pos.x >> 
      sites[i].pos.y >> 
      sites[i].pos.z;
    vec3 &pos = sites.back().pos;
    sites[i].idx = i;
    min = mineach(min, sites[i].pos);
    max = maxeach(max, sites[i].pos);
  }
#else
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
          min = mineach(min, sites.back().pos);
          max = maxeach(max, sites.back().pos);
        }
  }
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
#else /* reflecting */
  const real f = 0.5;
  const real dx = f*lx;
  const real dy = f*ly;
  const real dz = f*lz;
  const real ff = 1.0e-9;
  const real flx = ff*lx;
  const real fly = ff*ly;
  const real flz = ff*lz;
  for (int i = 0; i < np; i++)
  {
    const Voronoi::Site &s0 = sites[i];
    sitesP.push_back(s0);
    for (int oct = 1; oct < 8; oct++)
    {
      Voronoi::Site s = s0;
#if 1
      if (oct > 3) continue;
      if (oct == 1)
      {
        if (s.pos.x <  0.5*lx) s.pos.x = -s.pos.x;
        else                   s.pos.x = 2.0*lx - s.pos.x;
      }
      if (oct == 2)
      {
        if (s.pos.y <  0.5*ly) s.pos.y = -s.pos.y;
        else                   s.pos.y = 2.0*ly - s.pos.y;
      }
      if (oct == 3)
      {
        if (s.pos.z <  0.5*lz) s.pos.z = -s.pos.z;
        else                   s.pos.z = 2.0*lz - s.pos.z;
      }
      if (1)
      {
        s.idx = -1-s.idx;
#if 0
        const real dx = flx*(1.0-2.0*drand48());
        const real dy = fly*(1.0-2.0*drand48());
        const real dz = flz*(1.0-2.0*drand48());
        if (s.pos.x + dx > 0.0 && s.pos.x + dx < lx) s.pos.x += dx;
        if (s.pos.y + dy > 0.0 && s.pos.y + dy < ly) s.pos.y += dy;
        if (s.pos.z + dz > 0.0 && s.pos.z + dz < lz) s.pos.z += dz;
#endif
        sitesP.push_back(s);
      }
#else
      if (oct&1)
      {
        if  (s.pos.x <  0.5*lx) s.pos.x =        - s.pos.x;
        else                   s.pos.x = 2.0*lx - s.pos.x;
      }
      if (oct&2)
      {
        if   (s.pos.y <  0.5*ly) s.pos.y =        - s.pos.y;
        else                     s.pos.y = 2.0*ly - s.pos.y;
      }
      if (oct&4)
      {
        if  (s.pos.z < 0.5*lz) s.pos.z =        - s.pos.z;
        else                   s.pos.z = 2.0*lz - s.pos.z;
      }
      if (1)
      {
        s.idx = -1-s.idx;
        sitesP.push_back(s);
      }
#endif
    }
  }
#endif
#if 0
  for (int i = 0; i < (int)sitesP.size(); i++)
  {
    fprintf(stdout, "%g %g %g \n",
        sitesP[i].pos.x,
        sitesP[i].pos.y,
        sitesP[i].pos.z);
  }
  assert(0);
#endif
  fprintf(stderr, " np= %d  Pnp= %d\n", (int)sites.size(), (int)sitesP.size());
  assert(!sitesP.empty());

  double dt_search = 0.0;
  double dt_voro   = 0.0;
  double nface     = 0.0;
  const double tbeg = get_wtime();
  double volume    = 0.0;
  int nfailed = 0;

  for (int cnt = 0; cnt < 10; cnt++)
  {
    volume    = 0.0;
    nfailed = 0;

#pragma omp parallel reduction(+:dt_search, dt_voro, nface, volume, nfailed)
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
#if 1
        std::nth_element(dist.begin(), dist.begin() + ns, dist.end(), cmp_data<real, int>());
#else
        std::sort(dist.begin(), dist.end(), cmp_data<real, int>());
#endif
        list.clear();
        for (int j = 0; j < ns+1; j++)
          if (dist[j].first > 0.0f)
          {
            assert(dist[j].first <= dist[ns].first);
            list.push_back(Voronoi::Site(sitesP[dist[j].second].pos - s.pos, sitesP[dist[j].second].idx));
#if 0
            fprintf(stdout, " %g %g %g   %g dist= %g  j= %d\n", 
                list.back().pos.x,
                list.back().pos.y,
                list.back().pos.z,
                list.back().pos.norm2(),
                dist[j].first, j
                );
#endif
          }
        assert((int)list.size() == ns);
        double t1 = get_wtime();
        dt_search += t1 - t0;

        t0 = t1;
        if (!cell.build(list))
        {
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
  fprintf(stderr , " dt_voro  =  %g sec \n", dt_voro);
  fprintf(stderr,  "   dt_00=  %g \n" ,dt_00);
  fprintf(stderr,  "   dt_10=  %g \n" ,dt_10);
  fprintf(stderr,  "   dt_20=  %g \n" ,dt_20);
  fprintf(stderr,  "   dt_30=  %g \n" ,dt_30);
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
  fprintf(stderr, "  ratio= %g \n", incompl*1.0/incomplT);

  fprintf(stderr, " volume= %g   exact= %g  diff= %g \n",
      volume, lx*ly*lz, (volume-lx*ly*lz)/(lx*ly*lz));

  fprintf(stderr,  "ncell= %d  nfailed= %d : %g \n",      np, nfailed, (real)nfailed/np);


  return 0;
};
