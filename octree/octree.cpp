#include <cstdio>
#include <cstdlib>
#include "octree.h"
#include "plummer.h"
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
  {
    const int nb_mean = 32;
    const real s = std::pow(3.0/(4.0*M_PI)*(double)nb_mean/(double)n_bodies, 1.0/3.0);
    for (int i = 0; i < n_bodies; i++)
    {
      ptcl.push_back(Particle(
            vec3(drand48(), drand48(), drand48()),
            1.0/n_bodies, s));
    }
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
    tree.insert(octBodies[i]);
  }
  fprintf(stderr, "ncell= %d nnode= %d nleaf= %d n_nodes= %d  depth= %d\n",
      tree.get_ncell(), tree.get_nnode(), tree.get_nleaf(), n_nodes, tree.get_depth());
  const double t30 = get_wtime();

  fprintf(stderr, " -- Dump morton -- \n");
  std::vector<int> morton_list;
  morton_list.reserve(n_bodies);
  tree.tree_dump<true>(morton_list);
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
        octBodiesSorted.back().pos().x,
        octBodiesSorted.back().pos().y,
        octBodiesSorted.back().pos().z);
#endif
  }

  const double t50 = get_wtime();
  fprintf(stderr, " -- Buidling octTreeSorted -- \n");
  tree.clear();
  for (int i = 0; i < n_bodies; i++)
  {
    tree.insert(octBodiesSorted[i]);
  }
  fprintf(stderr, "ncell= %d nnode= %d nleaf= %d n_nodes= %d  depth= %d\n",
      tree.get_ncell(), tree.get_nnode(), tree.get_nleaf(), n_nodes, tree.get_depth());

  const double t60 = get_wtime();
  fprintf(stderr, " -- Compute boundaries -- \n");
  tree.computeBoundaries<true>();
  const boundary root_innerBnd = tree.root_innerBoundary();
  fprintf(stderr, " rootBnd_inner= %g %g %g  size= %g %g %g \n",
      root_innerBnd.center().x,
      root_innerBnd.center().y,
      root_innerBnd.center().z,
      root_innerBnd.hlen().x,
      root_innerBnd.hlen().y,
      root_innerBnd.hlen().z);
  const boundary root_outerBnd = tree.root_outerBoundary();
  fprintf(stderr, " rootBnd_outer= %g %g %g  size= %g %g %g \n",
      root_outerBnd.center().x,
      root_outerBnd.center().y,
      root_outerBnd.center().z,
      root_outerBnd.hlen().x,
      root_outerBnd.hlen().y,
      root_outerBnd.hlen().z);
  fprintf(stderr, "rootCentre:= %g %g %g  rootSize= %g \n",
      tree.get_rootCentre().x,
      tree.get_rootCentre().y,
      tree.get_rootCentre().z,
      tree.get_rootSize());

  const double t63 = get_wtime();
  assert(tree.sanity_check<true>() == n_bodies);
  const double t65 = get_wtime();
  tree.buildLeafList<true>();
  const double t68 = get_wtime();

  fprintf(stderr, " -- Range search -- \n");
  int nb = 0;
#if 1
#pragma omp parallel for reduction(+:nb)
  for (int i = 0; i < n_bodies; i++)
  {
#if 0
    nb += tree.range_search<true>(octBodies[i]);
#else
    nb += tree.range_search<true>(octBodiesSorted[i]);
#endif
  }
#endif
  const double t70 = get_wtime();
  const int nleaf = tree.nLeaf();
  fprintf(stderr, " -- Range search Leaf-Leaf : nleaf=%d  nbody= %d-- \n", nleaf, n_bodies);
  int nbL = 0;
#if 1
#pragma omp parallel for reduction(+:nbL)
  for (int i = 0; i < nleaf; i++)
  {
    int nb[Octree::NLEAF];
    const Octree::Leaf& leaf = tree.getLeaf(i);
    tree.range_search<true>(nb, leaf);
    for (int j = 0; j < leaf.nb(); j++)
      nbL += nb[j];
  }
#endif
  const double t75 = get_wtime();

  fprintf(stderr, " -- Remove ptcl -- \n");
  int nrm = 0;
#if 1
  for (int i = 0; i < n_bodies; i++)
  {
    tree.remove(Octree::Body(ptcl[i].pos, i));
    nrm++;
  }
  fprintf(stderr, "ncell= %d nnode= %d nleaf= %d n_nodes= %d  depth= %d\n",
      tree.get_ncell(), tree.get_nnode(), tree.get_nleaf(), n_nodes, tree.get_depth());
#endif

  const double t80 = get_wtime();

  fprintf(stderr, " -- Insert ptcl -- \n");
  int nins = 0;
#if 1
  for (int i = 0; i < n_bodies; i++)
  {
    tree.insert(octBodiesSorted[i]);
    nins++;
  }
  fprintf(stderr, "ncell= %d nnode= %d nleaf= %d n_nodes= %d  depth= %d\n",
      tree.get_ncell(), tree.get_nnode(), tree.get_nleaf(), n_nodes, tree.get_depth());
#endif

  const double t90 = get_wtime();

  fprintf(stderr, " -- Range search 1 -- \n");
  int nb1 = 0;
#if 1
#pragma omp parallel for reduction(+:nb)
  for (int i = 0; i < n_bodies; i++)
  {
#if 0
    nb1 += tree.range_search<true>(octBodies[i]);
#else
    nb1 += tree.range_search<true>(octBodiesSorted[i]);
#endif
  }
#endif
  const double t100 = get_wtime();


  fprintf(stderr, " Timing info: \n");
  fprintf(stderr, " -------------\n");
  fprintf(stderr, "   Plummer:  %g sec \n", t10 -t00);
  fprintf(stderr, "   BBox:     %g sec \n", t20 -t10);
  fprintf(stderr, "   Tree:     %g sec \n", t30 -t20);
  fprintf(stderr, "   Morton:   %g sec \n", t40 -t30);
  fprintf(stderr, "   Shuffle:  %g sec \n", t50 -t40);
  fprintf(stderr, "   TreeSort: %g sec \n", t60 -t50);
  fprintf(stderr, "   Boundary: %g sec \n", t63 -t60);
  fprintf(stderr, "   Sanity:   %g sec \n", t65 -t63);
  fprintf(stderr, "   LeafList: %g sec \n", t68 -t65);
  fprintf(stderr, "   RangeS:   %g sec <nb>= %g \n", t70 -t68, (real)nb /n_bodies);
  fprintf(stderr, "   RangeL:   %g sec <nb>= %g \n", t75 -t70, (real)nbL/n_bodies);
  fprintf(stderr, "   Remove:   %g sec nrm= %d \n", t80 - t70, nrm);
  fprintf(stderr, "   Insert:   %g sec nins= %d \n", t90 - t80, nins);
  fprintf(stderr, "   RangeS1:   %g sec <nb>= %g \n", t100 -t90, (real)nb1/n_bodies);

  return 0;
}
