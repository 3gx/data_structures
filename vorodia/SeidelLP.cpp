#include "SeidelLP.h"
#include "mytimer.h"

int main(int argc, char * argv [])
{

  const int N = 100000;
  SeidelLP<N> lp;

#if 0
  srand48(12345);
  int n = 1000;
  for (int i = 0; i < n; i++)
  {
    const real lx = 0.1;
    const real ly = 0.1;
    const real lz = 0.1;
    const real nx = (1.0-2.0*drand48())*lx;
    const real ny = (1.0-2.0*drand48())*ly;
    const real nz = (1.0-2.0*drand48())*lz;
    const real x = -nx + 0.5;
    const real y = -ny + 0.5;
    const real z = -nz + 0.5;
    lp.push(HalfSpace(vec3(nx,ny,nz), vec3(x, y, z)));
  }
#else
    lp.push(HalfSpace(vec3(+1.0, 0.0, 0.0), vec3(0.25, 0.5, 0.5)));
    lp.push(HalfSpace(vec3(-1.0, 0.0, 0.0), vec3(0.75, 0.5, 0.5)));
    lp.push(HalfSpace(vec3(0.0, +1.0, 0.0), vec3(0.5, 0.25, 0.5)));
    lp.push(HalfSpace(vec3(0.0, -1.0, 0.0), vec3(0.5, 0.75, 0.5)));
    lp.push(HalfSpace(vec3(0.0, 0.0, +1.0), vec3(0.5, 0.5, 0.25)));
    lp.push(HalfSpace(vec3(0.0, 0.0, -1.0), vec3(0.5, 0.5, 0.75)));
    lp.push(HalfSpace(vec3(0.0, 0.0, +1.0), vec3(0.5, 0.5, 0.40)));
    lp.push(HalfSpace(vec3(1.0/sqrt(2.0), 1.0/sqrt(2.0), 0.0), vec3(0.3, 0.3, 0.0)));
#endif

  fprintf(stderr, " nspace= %d\n", lp.nspace());


  const double t0 = get_wtime();
  const vec3 pos = lp.solve();
  const double dt = get_wtime() - t0;
  fprintf(stderr, " pos= %g %g %g   done in %g sec\n", pos.x, pos.y, pos.z, dt);

  return 0;
}

