#include "vorocell.h"
#include <iostream>

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
  for (int i = 0; i < np; i++)
  {
    std::cin >> idum >> 
      sites[i].pos.x >> 
      sites[i].pos.y >> 
      sites[i].pos.z;
    sites[i].idx = i;
  }





  return 0;
};
