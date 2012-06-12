#include "voroCell.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "mytimer.h"

int main(int argc, char * argv[])
{
  assert(argc > 2);
  const int Nsite   = atoi(argv[1]);
  const int Nrepeat = atoi(argv[2]);

//  srand48(3456);

  fprintf(stderr, "nsites= %d nrepeat= %d\n", Nsite, Nrepeat);

  Voro::Site::Vector sites(Nsite);
  for (int i = 0; i < Nsite; i++)
    sites[i] = Voro::Site(dvec3(drand48(), drand48(), drand48()), i);

  Voro::Site::Vector searchSites(Nrepeat);
  for (int i = 0; i < Nrepeat; i++)
  {
    const double x = 0.5 + (1.0 - 2.0*drand48())*0.1;
    const double y = 0.5 + (1.0 - 2.0*drand48())*0.1;
    const double z = 0.5 + (1.0 - 2.0*drand48())*0.1;
    searchSites[i] = Voro::Site(dvec3(x, y, z), -1);
  }

  Voro::Cell cell;
  const double t00 = get_wtime();
  for (int i = 0; i < Nrepeat; i++)
  {
    cell.build(searchSites[i], sites);
    if (i == Nrepeat-1)
      fprintf(stderr , "nface= %d  r= %g\n", cell.nFace(), std::sqrt(cell.r2()));
  }
  const double t10 = get_wtime();
  fprintf(stderr, " %d done in %g sec  [rate: %g cell/sec]\n",
      Nrepeat, t10 - t00, Nrepeat/(t10-t00));



  return 0;

};

