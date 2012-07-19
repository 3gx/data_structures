#include "octree.h"
#include "DirectPoly.h"

#include "plummer.h"
#include "mytimer.h"

int main(int argc, char * argv[])
{
  int n_bodies = 10240;
  if (argc > 1) n_bodies = atoi(argv[1]);
  assert(n_bodies > 0);
  fprintf(stderr, "n_bodies= %d\n", n_bodies);

  Particle::Vector ptcl;
  ptcl.reserve(n_bodies);

#if 1
#define PLUMMER
  const Plummer data(n_bodies);
  for (int i = 0; i < n_bodies; i++)
  {
    ptcl.push_back(Particle(data.pos[i], data.mass[i]));
  }
#elif 0  /* reads IC */
  int dummy;
  std::cin >> dummy >> n_bodies;
  fprintf(stderr, " -- input file: nbodies= %d\n", n_bodies);
  for (int i = 0; i < n_bodies; i++)
  {
    Particle p;
    std::cin >> p.pos.x >> p.pos.y >> p.pos.z >> p.h >> p.nb;
    p.h *= 2.0;
    ptcl.push_back(p);
  }
#elif 1
  {
    const int nb_mean = 128;
    const real s = std::pow(3.0/(4.0*M_PI)*(double)nb_mean/(double)n_bodies, 1.0/3.0);
    for (int i = 0; i < n_bodies; i++)
    {
      ptcl.push_back(Particle(
            vec3(drand48(), drand48(), drand48()),
            1.0/n_bodies, s));
    }
  }
#endif

  return 0;
}
