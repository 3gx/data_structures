#include <cstdio>
#include <cstdlib>
#include "plummer.h"
#include "octree.h"
#include "mytimer.h"

int main(int argc, char * argv[])
{
  int n_bodies = 10240;
  if (argc > 1) n_bodies = atoi(argv[1]);
  assert(n_bodies > 0);
  fprintf(stderr, "n_bodies= %d\n", n_bodies);

  const double t00 = get_wtime();
  Particle::Vector ptcl;
  ptcl.reserve(n_bodies);
#if 0
  const Plummer data(n_bodies);
  for (int i = 0; i < n_bodies; i++)
  {
    ptcl.push_back(Particle(data.pos[i], data.mass[i]));
  }
#else
  for (int i = 0; i < n_bodies; i++)
  {
    ptcl.push_back(Particle(
          vec3(drand48(), drand48(), drand48()),
          1.0/n_bodies));
  }
#endif
  
  Octree::Body::Vector octBodies;
  octBodies.reserve(n_bodies);

  const double t10 = get_wtime();
  fprintf(stderr, " -- Compute bounding box -- \n");

  vec3 rmin(+HUGE);
  vec3 rmax(-HUGE);

  for (int i = 0; i < n_bodies; i++)
  {
    octBodies.push_back(Octree::Body(ptcl[i], i));
    rmin = mineach(rmin, ptcl[i].pos);
    rmax = maxeach(rmax, ptcl[i].pos);
  }
  const vec3 centre = (rmax + rmin)*0.5;
  const vec3 vsize  =  rmax - rmin;
  const real  size  = __max(__max(vsize.x, vsize.y), vsize.z);
  real size2 = 1.0;
  while (size2 > size) size2 *= 0.5;
  while (size2 < size) size2 *= 2.0;

  const int n_nodes = n_bodies;
  Octree tree(centre, size2, n_nodes);

  const double t20 = get_wtime();
  fprintf(stderr, " -- Buidling octTree -- \n");
  for (int i = 0; i < n_bodies; i++)
  {
    tree.push(octBodies[i], i, octBodies);
  }
  fprintf(stderr, "ncell= %d n_nodes= %d  depth= %d\n",
      tree.ncell, tree.n_nodes, tree.depth);
  double t30 = get_wtime();
  assert(tree.sanity_check(octBodies) == n_bodies);
  t30 = get_wtime();
  
  fprintf(stderr, " -- Dump morton -- \n");
  std::vector<int> morton_list;
  morton_list.reserve(n_bodies);
  tree.morton_dump(morton_list);
  fprintf(stderr, "morton_list.size()= %d  n_bodies= %d\n",
      (int)morton_list.size(), n_bodies);
  for (std::vector<int>::iterator it = morton_list.begin(); it != morton_list.end(); it++)
    assert(*it < n_bodies);
  assert((int)morton_list.size() == n_bodies);
  const double t40 = get_wtime();
  
  fprintf(stderr, " -- Shuffle octBodies -- \n");
  Octree::Body::Vector octBodiesSorted;
  octBodiesSorted.reserve(n_bodies);
  for (std::vector<int>::iterator it = morton_list.begin(); it != morton_list.end(); it++)
  {
    assert(*it < n_bodies);
    octBodiesSorted.push_back(octBodies[*it]);
#if 0
    fprintf(stdout, "%g %g %g \n", 
        octBodiesSorted.back().pos.x,
        octBodiesSorted.back().pos.y,
        octBodiesSorted.back().pos.z);
#endif
  }
  const double t50 = get_wtime();

  fprintf(stderr, " -- Buidling octTreeSorted -- \n");
  tree.clear();
  for (int i = 0; i < n_bodies; i++)
  {
    tree.push(octBodiesSorted[i], i, octBodiesSorted);
  }
  fprintf(stderr, "ncell= %d n_nodes= %d  depth= %d\n",
      tree.ncell, tree.n_nodes, tree.depth);
  const double t60 = get_wtime();
  fprintf(stderr, " -- Sanity check -- \n");
  assert(tree.sanity_check(octBodiesSorted) == n_bodies);
  const double t63 = get_wtime();
  fprintf(stderr, " -- Inner boundary -- \n");
  const boundary rootBnd = tree.inner_boundary(octBodiesSorted);
  fprintf(stderr, " bnd= %g %g %g  size= %g %g %g \n",
      rootBnd.center().x,
      rootBnd.center().y,
      rootBnd.center().z,
      rootBnd.hlen().x,
      rootBnd.hlen().y,
      rootBnd.hlen().z);
  fprintf(stderr, " c= %g %g %g size= %g\n",
      tree.root_centre.x,
      tree.root_centre.y,
      tree.root_centre.z,
      tree.root_size*0.5);
  
  const double t67 = get_wtime();
  assert(tree.sanity_check1(octBodiesSorted) == n_bodies);
  const double t68 = get_wtime();

  fprintf(stderr, " -- Range search -- \n");
  int nb = 0;
#if 0
  const int nb_mean = 32;
  const real s = std::pow(3.0/(4.0*M_PI)*(double)nb_mean/(double)n_bodies, 1.0/3.0);
#pragma omp parallel for reduction(+:nb)
  for (int i = 0; i < n_bodies; i++)
  {
    nb += tree.range_search(ptcl[i].pos, s, octBodiesSorted);
  }
#endif
  const double t70 = get_wtime();



  fprintf(stderr, " Timing info: \n");
  fprintf(stderr, " -------------\n");
  fprintf(stderr, "   Plummer:  %g sec \n", t10 -t00);
  fprintf(stderr, "   BBox:     %g sec \n", t20 -t10);
  fprintf(stderr, "   Tree:     %g sec \n", t30 -t20);
  fprintf(stderr, "   Morton:   %g sec \n", t40 -t30);
  fprintf(stderr, "   Shuffle:  %g sec \n", t50 -t40);
  fprintf(stderr, "   TreeSort: %g sec \n", t60 -t50);
  fprintf(stderr, "   Sanity:   %g sec \n", t63 -t60);
  fprintf(stderr, "   Boundary: %g sec \n", t67 -t63);
  fprintf(stderr, "   Sanity1:   %g sec \n", t68 -t67);
  fprintf(stderr, "   RangeS:   %g sec <nb>= %g \n", t70 -t68, (real)nb/n_bodies);

  return 0;
}
