#include "voroCell.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>

int main(int argc, char * argv[])
{
  assert(argc > 1);
  const int n = atoi(argv[1]);

  fprintf(stderr, "nsites= %d\n", n);

  Voro::Site::Vector sites(n);
  for (int i = 0; i < n; i++)
    sites[i] = Voro::Site(dvec3(drand48(), drand48(), drand48()), i);

  const Voro::Site s0(dvec3(0.5, 0.5, 0.5), -1);
  const Voro::Cell cell(s0, sites);

  fprintf(stderr , "nface= %d  r= %g\n", cell.nFace(), std::sqrt(cell.r2()));


  return 0;
  
};

