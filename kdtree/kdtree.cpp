#include "kdtree.h"
#include <cstdio>
#include <cstdlib>
#include "plummer.h"
#include "mytimer.h"

int main(int argc, char * argv[])
{
  int n_bodies = 10240;
  if (argc > 1) n_bodies = atoi(argv[1]);
  assert(n_bodies > 0);
  fprintf(stderr, "n_bodies= %d\n", n_bodies);

  const double t00 = get_wtime();
#if 0
  const Plummer data(n_bodies);
#endif
  const double t10 = get_wtime();

  Particle::Vector ptcl;
  ptcl.reserve(n_bodies);
  for (int i = 0; i < n_bodies; i++)
  {
#if 0
    ptcl.push_back(Particle(data.pos[i], data.mass[i]));
#else
    ptcl.push_back(Particle(
          vec3(drand48(), drand48(), drand48()),
          1.0/n_bodies));
#endif
  }


  const double t20 = get_wtime();
  fprintf(stderr, " -- Buidling kdTree -- \n");
  kdTree tree(ptcl);
  const double t30 = get_wtime();

  fprintf(stderr, " -- Searching range ngb -- \n");
  std::vector<int> nblist;
  nblist.reserve(1024);
  int nb = 0;
  const int nb_mean = 32;
  const real s = std::pow(3.0/(4.0*M_PI)*(double)nb_mean/(double)n_bodies, 1.0/3.0);
#pragma omp parallel for reduction(+:nb)
  for (int i = 0; i < n_bodies; i++)
  {
    nb += tree.find_range_nb(ptcl[i].pos, s);
  }
  const double t40 = get_wtime();

  const int K = 8;
#if 0
  fprintf(stderr, " -- Searching k-nearest ngb -- \n");
  for (int i = 0; i < n_bodies; i++)
  {
    int klist[K];
    tree.find_knb<K>(ptcl[i].pos, klist);
  }
#endif
  const double t50 = get_wtime();

  fprintf(stderr, " -- Searching nearest ngb -- \n");
  {
    real s = 0.0;
    real s2 = 0.0;
#pragma omp parallel for reduction(+:s,s2)
    for (int i = 0; i < n_bodies; i++)
    {
      const int node_id = tree.find_nnb(ptcl[i].pos);
      const int  j = tree[node_id].body_idx();
      const real r = (ptcl[i].pos - ptcl[j].pos).abs();
#if 0 /* correctness check */
      real s2min = HUGE;
      int  jmin  = -1;
      for (int jx = 0; jx < n_bodies; jx++)
      {
        const real r2 = (ptcl[i].pos - ptcl[jx].pos).norm2();
        if (r2 < s2min && r2 != 0.0f)
        {
          s2min = r2;
          jmin = jx;
        } 
      }
      assert(jmin == j);
#endif
      s  += r;
      s2 += r*r;
    }
    s  *= 1.0/n_bodies;
    s2 *= 1.0/n_bodies;
    const real ds = std::sqrt(s2 - s*s);
    fprintf(stderr, "<r> = %g  sigma= %g \n", s, ds);
  }
  const double t60 = get_wtime();

  fprintf(stderr, " Timing info: \n");
  fprintf(stderr, " -------------\n");
  fprintf(stderr, "   Plummer:  %g sec \n", t10 -t00);
  fprintf(stderr, "   Copy:     %g sec \n", t20 -t10);
  fprintf(stderr, "   kdTree:   %g sec \n", t30 -t20);
  fprintf(stderr, "   ngb :     %g sec, nb= %g \n", t40 -t30, (real)nb/(real)n_bodies);
  fprintf(stderr, "   kgb :     %g sec, K= %d \n", t50 -t40, K);
  fprintf(stderr, "   ngb :     %g sec \n", t60 -t50);



  return 0;
}
