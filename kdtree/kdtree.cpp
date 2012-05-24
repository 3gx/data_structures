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
  const Plummer data(n_bodies);
  const double t10 = get_wtime();

  Particle::Vector ptcl;
  ptcl.reserve(n_bodies);
  for (int i = 0; i < n_bodies; i++)
    ptcl.push_back(Particle(data.pos[i], data.mass[i]));


  const double t20 = get_wtime();
  fprintf(stderr, " -- Buidling kdTree -- \n");
  kdTree tree(ptcl);
  const double t30 = get_wtime();

  fprintf(stderr, " -- Searching range ngb -- \n");
  std::vector<int> nblist;
  nblist.reserve(1024);
  int nb = 0;
  const real s = 0.1;
  const real s2 = s*s;
#pragma omp parallel for reduction(+:nb)
  for (int i = 0; i < n_bodies; i++)
  {
    std::vector<int> nblist;
    nblist.reserve(1024);
    tree.find_range_nb(ptcl[i].pos, s2, nblist);
    nb += nblist.size();
  }
  const double t40 = get_wtime();

  const int K = 16;
#if 0
  fprintf(stderr, " -- Searching k-nearest ngb -- \n");
  int klist[K];
  for (int i = 0; i < n_bodies; i++)
  {
    tree.find_knb<K>(ptcl[i].pos, klist);
  }
#endif
  const double t50 = get_wtime();
  fprintf(stderr, " Timing info: \n");
  fprintf(stderr, " -------------\n");
  fprintf(stderr, "   Plummer:  %g sec \n", t10 -t00);
  fprintf(stderr, "   Copy:     %g sec \n", t20 -t10);
  fprintf(stderr, "   kdTree:   %g sec \n", t30 -t20);
  fprintf(stderr, "   ngb :     %g sec, nb= %g \n", t40 -t30, (real)nb/(real)n_bodies);
  fprintf(stderr, "   kgb :     %g sec, K= %d \n", t50 -t40, K);



  return 0;
}
