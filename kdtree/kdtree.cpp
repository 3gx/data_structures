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

  fprintf(stderr, " Timing info: \n");
  fprintf(stderr, " -------------\n");
  fprintf(stderr, "   Plummer:  %g sec \n", t10 -t00);
  fprintf(stderr, "   Copy:     %g sec \n", t20 -t10);
  fprintf(stderr, "   kdTree:   %g sec \n", t30 -t20);



  return 0;
}
