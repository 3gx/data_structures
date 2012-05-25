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
#if 1
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
  assert(tree.sanity_check() == n_bodies);
  const double t30 = get_wtime();
  
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
  assert(tree.sanity_check() == n_bodies);
  const double t60 = get_wtime();

  fprintf(stderr, " Timing info: \n");
  fprintf(stderr, " -------------\n");
  fprintf(stderr, "   Plummer:  %g sec \n", t10 -t00);
  fprintf(stderr, "   BBox:     %g sec \n", t20 -t10);
  fprintf(stderr, "   Tree:     %g sec \n", t30 -t20);
  fprintf(stderr, "   Morton:   %g sec \n", t40 -t30);
  fprintf(stderr, "   Shuffle:  %g sec \n", t50 -t40);
  fprintf(stderr, "   TreeSort: %g sec \n", t60 -t50);

  return 0;
}
