#include "octree.h"
#include "DirectPoly.h"

#include "plummer.h"
#include "mytimer.h"

enum {NGROUP = 256};
typedef Octree::GroupT<NGROUP> octGroup;

int main(int argc, char * argv[])
{
  int n_bodies = 10240;
  if (argc > 1) n_bodies = atoi(argv[1]);
  assert(n_bodies > 0);
  fprintf(stderr, "n_bodies= %d\n", n_bodies);

  Particle::Vector ptcl;
  ptcl.reserve(n_bodies);

  /*  generate or import particles */

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
  
  /* construct the tree */

  Octree::Body::Vector octBodies;
  octBodies.reserve(n_bodies);

  fprintf(stderr, " -- Compute bounding box -- \n");

  vec3 rmin(+HUGE);
  vec3 rmax(-HUGE);

  for (int i = 0; i < n_bodies; i++)
  {
    octBodies.push_back(Octree::Body(ptcl[i], i));
    rmin = mineach(rmin, ptcl[i].pos);
    rmax = maxeach(rmax, ptcl[i].pos);
  }
  vec3 centre = (rmax + rmin)*0.5;
  const vec3 vsize  =  rmax - rmin;
  const real  size  = __max(__max(vsize.x, vsize.y), vsize.z);
  real size2 = 1.0;
  while (size2 > size) size2 *= 0.5;
  while (size2 < size) size2 *= 2.0;
  centre = 0.0;
  rmin = vec3(-size2);
  rmax = vec3(+size2);
  for (int i = 0; i < n_bodies; i++)
  {
    assert(ptcl[i].pos.x >= rmin.x);
    assert(ptcl[i].pos.y >= rmin.y);
    assert(ptcl[i].pos.z >= rmin.z);
    assert(ptcl[i].pos.x <= rmax.x);
    assert(ptcl[i].pos.y <= rmax.y);
    assert(ptcl[i].pos.z <= rmax.z);
  }


  const int n_nodes = n_bodies;
  Octree tree(centre, size2, n_nodes);
  fprintf(stderr, " centre= %g %g %g   size= %g\n",
      centre.x, centre.y, centre.z, size2);
  
  fprintf(stderr, " -- Buidling octTree -- \n");
  for (int i = 0; i < n_bodies; i++)
  {
    tree.insert(octBodies[i]);
  }
  fprintf(stderr, "ncell= %d nnode= %d nleaf= %d n_nodes= %d  depth= %d np= %d\n",
      tree.get_ncell(), tree.get_nnode(), tree.get_nleaf(), n_nodes, tree.get_depth(), tree.get_np());
  
  fprintf(stderr, " -- Dump morton -- \n");
  std::vector<int> morton_list;
  morton_list.reserve(n_bodies);
  tree.tree_dump<true>(morton_list);
  fprintf(stderr, "morton_list.size()= %d  n_bodies= %d\n",
      (int)morton_list.size(), n_bodies);
  for (std::vector<int>::iterator it = morton_list.begin(); it != morton_list.end(); it++)
    assert(*it < n_bodies);
  
  fprintf(stderr, " -- Shuffle octBodies -- \n");
  {
    Octree::Body::Vector octBodiesSorted;
    octBodiesSorted.reserve(n_bodies);
    for (std::vector<int>::iterator it = morton_list.begin(); it != morton_list.end(); it++)
    {
      assert(*it < n_bodies);
      octBodiesSorted.push_back(octBodies[*it]);
    }
    octBodies = octBodiesSorted;
  }

  fprintf(stderr, " -- Buidling octTreeSorted -- \n");
  tree.clear();
  for (int i = 0; i < n_bodies; i++)
    tree.insert(octBodies[i]);
  fprintf(stderr, "ncell= %d nnode= %d nleaf= %d n_nodes= %d  depth= %d\n",
      tree.get_ncell(), tree.get_nnode(), tree.get_nleaf(), n_nodes, tree.get_depth());
  
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
  fprintf(stderr, "rootCentre:= %g %g %g  rootSize= %g \n",
      tree.get_rootCentre().x,
      tree.get_rootCentre().y,
      tree.get_rootCentre().z,
      tree.get_rootSize());

  assert(tree.sanity_check<true>() == n_bodies);
  tree.buildLeafList<true>();

  /* ready to construct vorocells ... */
  
  const bool SORT = 0 ? true : false;  /* use peano-sort inside the group */
  octGroup::Vector groupList;
  tree.buildGroupList<SORT, true>(groupList);
 

  const int NFMAX = 1024;
  DirectPolyhedron<NFMAX> direct;
  const double t0 = get_wtime();
  for (int i = 0; i < n_bodies; i++)
  {
    direct.clear();
    const vec3 ipos = octBodies[i].vector_pos();

    const real f = 1.0;
    direct.push(vec3(2.0*(rmin.x-ipos.x), 0.0, 0.0), -1, f);
    direct.push(vec3(2.0*(rmax.x-ipos.x), 0.0, 0.0), -2, f);
    direct.push(vec3(0.0, 2.0*(rmin.y-ipos.y), 0.0), -3, f);
    direct.push(vec3(0.0, 2.0*(rmax.y-ipos.y), 0.0), -4, f);
    direct.push(vec3(0.0, 0.0, 2.0*(rmin.z-ipos.z)), -5, f);
    direct.push(vec3(0.0, 0.0, 2.0*(rmax.z-ipos.z)), -6, f);
#if 0
    tree.buildDirectPolyhedron(octBodies[i], direct, f);
#else
    for (int j= 0; j < n_bodies; j++)
      if (i != j)
        direct.push(octBodies[j].vector_pos() - ipos, j, f);
#endif
    fprintf(stderr, "i= %d  nface= %d\n", i, direct.nface());
  }
  const double t1 = get_wtime();
  fprintf(stderr, " -- done in %g sec -- \n", t1 - t0);

  return 0;

}
