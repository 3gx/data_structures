#include "SeidelLP.h"
#include "mytimer.h"

int main(int argc, char * argv [])
{
  int n = 1000;
  srand48(120);
//  srand48(123);

  if (argc > 1)
    srand48(atoi(argv[1]));

  if (argc > 2)
    n = atoi(argv[2]);
  fprintf(stderr, "n= %d\n", n);

  const int N = 100000;
  SeidelLP<N> lp(1.0);

#if 1
  for (int i = 0; i < n; i++)
  {
    const real lx = 0.1;
    const real ly = 0.1;
    const real lz = 0.1;
    const real nx = (1.0-2.0*drand48())*lx;
    const real ny = (1.0-2.0*drand48())*ly;
    const real nz = (1.0-2.0*drand48())*lz;
    const real x = -nx;
    const real y = -ny;
    const real z = -nz;
    lp.push(HalfSpace(vec3(nx,ny,nz), vec3(x, y, z)));
  }
#else
    lp.push(HalfSpace(vec3(+1.0, 0.0, 0.0), vec3(0.25, 0.5, 0.5)));
    lp.push(HalfSpace(vec3(-1.0, 0.0, 0.0), vec3(0.75, 0.5, 0.5)));
    lp.push(HalfSpace(vec3(0.0, +1.0, 0.0), vec3(0.5, 0.25, 0.5)));
    lp.push(HalfSpace(vec3(0.0, -1.0, 0.0), vec3(0.5, 0.75, 0.5)));
    lp.push(HalfSpace(vec3(0.0, 0.0, +1.0), vec3(0.5, 0.5, 0.25)));
    lp.push(HalfSpace(vec3(0.0, 0.0, -1.0), vec3(0.5, 0.5, 0.75)));
 //   lp.push(HalfSpace(vec3(0.0, 0.0, +1.0), vec3(0.5, 0.5, 0.40)));
//    lp.push(HalfSpace(vec3(0.0, 0.0, -1.0), vec3(0.5, 0.5, 0.30)));
//    lp.push(HalfSpace(vec3(0.0, +1.0, 0.0), vec3(0.5, 0.3, 0.30)));
    lp.push(HalfSpace(vec3(-1.0, -1.0,  0.0), vec3(0.3, 0.3, 0.0)));
#endif

  fprintf(stderr, " nspace= %d\n", lp.nspace());
  const int nrep = 100;


  vec3 cvec(-0.0, 0.0, +1.0);
  cvec = vec3(drand48(), drand48(), drand48());
  vec3 pos = lp.solve(cvec, false);
  {
    const double t0 = get_wtime();
    fprintf(stderr, " pos= %g %g %g \n", pos.x, pos.y, pos.z);
    for (int i = 0; i < nrep; i++)
    {
      vec3 pos1 = lp.solve(cvec, false);
      if ((pos1 - pos).norm2() > 1.0e-20*pos.norm2())
        fprintf(stderr, " pos= %g %g %g \n", pos1.x, pos1.y, pos1.z);
    }
    const double dt = get_wtime() - t0;
    fprintf(stderr, " norm  done in %g sec\n", dt/nrep);
  }
  {
    const double t0 = get_wtime();
    for (int i = 0; i < nrep; i++)
    {
      vec3 pos1 = lp.solve(cvec, true);
      if ((pos1 - pos).norm2() > 1.0e-20*pos.norm2())
        fprintf(stderr, " rpos= %g %g %g \n", pos1.x, pos1.y, pos1.z);
    }
    const double dt = get_wtime() - t0;
    fprintf(stderr, " rand  done in %g sec\n", dt/nrep);
  }

  return 0;
}

