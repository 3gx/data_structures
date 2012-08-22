#include <vector>
#include "vector3.h"
typedef double real;
typedef vector3<real> vec3;

template<class T>
inline T __min(const T a, const T b) {return a < b ? a : b;}

template<class T>
inline T __max(const T &a, const T &b) {return a > b ? a : b;}

template<class T>
inline T __sign(const T a) {return a < T(0.0) ? (T)-1.0 : (T)+1.0;}

template<class T>
inline T __abs(const T a) {return a < T(0.0) ? -a : a;}

struct Particle
{
  typedef std::vector<Particle> Vector;
  vec3 pos;
  int  id;

  Particle() {}
  Particle(const vec3 &_pos, const int _id) :
    pos(_pos), id(_id) {}
};
struct CmpDist
{
  vec3 ipos;
  CmpDist(const vec3 &vec) : ipos(vec) {}
  bool operator()(const Particle &p1, const Particle &p2) const
  {
    return (p1.pos-ipos).norm2() < (p2.pos - ipos).norm2();
  }
};


#include "MSW.h"
#include "SeidelLP.h"
#include "octree.h"
#include "vorocell.h"

#include "plummer.h"
#include "mytimer.h"

enum {NGROUP = 128};
typedef Tree::Octree::GroupT<NGROUP> octGroup;

struct CmpList
{
  vec3 vec;
  CmpList(const vec3 &v) : vec(v) {}
  bool operator()(const Voronoi::Site &s1, const Voronoi::Site &s2) const
  {
    return (s1.pos).norm2() < (s2.pos).norm2();
  }
};

int main(int argc, char * argv[])
{
  int n_bodies = 512;
  if (argc > 1) n_bodies = atoi(argv[1]);
  assert(n_bodies > 0);
  fprintf(stderr, "n_bodies= %d\n", n_bodies);

  Particle::Vector ptclO;
  ptclO.reserve(n_bodies);

  /*  generate or import particles */

#if 0
#define PLUMMER
  const Plummer data(n_bodies);
  for (int i = 0; i < n_bodies; i++)
  {
    ptclO.push_back(Particle(data.pos[i], data.mass[i]));
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
    ptclO.push_back(p);
  }
#elif 1
  {
    for (int i = 0; i < n_bodies; i++)
      ptclO.push_back(Particle(vec3(
              1.0-2.0*drand48(), 
              1.0-2.0*drand48(), 
              1.0-2.0*drand48())*0.5,
            i));
  }
#endif

  Particle::Vector ptclP;
  ptclP.reserve(2*n_bodies);
  const vec3 rminD(-0.5);
  const vec3 rmaxD(+0.5);
  {
#if 0 /* periodic */
    const real  f = 0.5;
    assert(f <= 0.5);
    const real lx = rmaxD.x - rminD.x;
    const real ly = rmaxD.y - rminD.y;
    const real lz = rmaxD.z - rminD.z;
    const real dx = (0.5 - f) * lx;
    const real dy = (0.5 - f) * ly;
    const real dz = (0.5 - f) * lz;
    const vec3 cpos(0.5*lx, 0.5*ly, 0.5*lz);
    for (int i = 0; i < n_bodies; i++)
    {
      const Particle &p0 = ptclO[i];
      ptclP.push_back(p0);
      for (int oct = 1; oct < 8; oct++)
      {
        Particle s = p0;
        s.pos -= rminD;
        int ioct = 0;
        if ((oct&1) && __abs(s.pos.x-cpos.x) > dx) {s.pos.x += lx * (s.pos.x <= cpos.x ? +1.0 : -1.0); ioct += 1;}
        if ((oct&2) && __abs(s.pos.y-cpos.y) > dy) {s.pos.y += ly * (s.pos.y <= cpos.y ? +1.0 : -1.0); ioct += 2;}
        if ((oct&4) && __abs(s.pos.z-cpos.z) > dz) {s.pos.z += lz * (s.pos.z <= cpos.z ? +1.0 : -1.0); ioct += 4;}
        s.pos += rminD;
        if (oct == ioct)
          ptclP.push_back(s);
      }
    }
#else  /* reflecting */
#define REFLECTING
    ptclP = ptclO;
#endif
  }


  /* construct the tree */

  Tree::Octree::Body::Vector octBodies;
  octBodies.reserve(n_bodies);

  fprintf(stderr, " -- Compute bounding box -- \n");

  const int n_bodies0 = n_bodies;
  n_bodies = ptclP.size();
  fprintf(stderr, "nbodiesP= %d\n", n_bodies);

  vec3 rminT(+HUGE);
  vec3 rmaxT(-HUGE);

  for (int i = 0; i < n_bodies; i++)
  {
    octBodies.push_back(Tree::Octree::Body(ptclP[i], i));
    rminT = mineach(rminT, ptclP[i].pos);
    rmaxT = maxeach(rmaxT, ptclP[i].pos);
  }
  const vec3 centreT = (rmaxT + rminT)*0.5;
  const vec3 vsizeT  =  rmaxT - rminT;
  const real  sizeT  = __max(__max(vsizeT.x, vsizeT.y), vsizeT.z);
  fprintf(stderr, " centre= %g %g %g  size= %g \n", centreT.x, centreT.y, centreT.z, sizeT);
  real size2T = 1.0;
  while (size2T > sizeT) size2T *= 0.5;
  while (size2T < sizeT) size2T *= 2.0;

  const int n_nodes = n_bodies;
  Tree::Octree tree(centreT, size2T, n_nodes);
  fprintf(stderr, " >>> centre= %g %g %g   size= %g\n",
      centreT.x, centreT.y, centreT.z, size2T);

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
    Tree::Octree::Body::Vector octBodiesSorted;
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
  const Tree::boundary root_innerBnd = tree.root_innerBoundary();
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


  double t0 = get_wtime();
  int nfaceD_min = 1<<30;
  int nfaceD_max = 0;
  unsigned long long nfaceD= 0;
  int nfaceL_min = 1<<30;
  int nfaceL_max = 0;
  unsigned long long nfaceL= 0;
  int nfaceV_min = 1<<30;
  int nfaceV_max = 0;
  unsigned long long nfaceV= 0;

  std::vector<int> idx(n_bodies);
  for (int i = 0; i < n_bodies; i++)
    idx[i] = i;

#if 0
  t0 = get_wtime();
#endif

#if 0
  const int NFMAX = 1024;
  DirectPolyhedron<NFMAX> direct;
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
    std::random_shuffle(idx.begin(), idx.end());
    for (int j= 0; j < n_bodies; j++)
      if (i != idx[j])
        direct.push(octBodies[idx[j]].vector_pos() - ipos, idx[j], f);
    nface_min = std::min(nface_min, direct.nface());
    nface_max = std::max(nface_max, direct.nface());
    nfaceS   += direct.nface();
    //    fprintf(stderr, "i= %d  nface= %d\n", i, direct.nface());
  }
#else
  const int ngroup = groupList.size();
  double volume = 0.0;
  int np = 0;
  {
    const int NFMAX = 1024;
    DirectPolyhedron<NFMAX> direct;
    Voronoi::Site::Vector list;
    const int NF=1030;
    Voronoi::Cell<NF> *cell_ptr = new Voronoi::Cell<NF>;
    Voronoi::Cell<NF> &cell = *cell_ptr;
    std::vector<bool> used(n_bodies, false);
#if 1
    MSW lp(4*sizeT); //2.0); //;.0*size2);
#else
    SeidelLP lp(4*sizeT); //2.0); //;.0*size2);
#endif
    for (int igroup = 0; igroup < ngroup; igroup++)
    {
      const octGroup &group = groupList[igroup];
      const int ni = group.nb();
      for (int i = 0; i < ni; i++)
      {
        const int    ix = group[i].idx();
        const vec3 ipos = ptclP[ix].pos;
        if (
            ipos.x < rminD.x || ipos.x > rmaxD.x ||
            ipos.y < rminD.y || ipos.y > rmaxD.y ||
            ipos.z < rminD.z || ipos.z > rmaxD.z)
          continue;
        np++;

        assert(ipos.x > rminD.x);
        assert(ipos.y > rminD.y);
        assert(ipos.z > rminD.z);
        assert(ipos.x < rmaxD.x);
        assert(ipos.y < rmaxD.y);
        assert(ipos.z < rmaxD.z);

        direct.clear();
        std::stack<int> plist;
#if 0
        const int nj = ni;
        for (int j = 0; j < nj; j++)
        {
          const int jx = group[j].idx();
          used[jx] = true;
          const vec3 dr = ptclP[jx].pos - ipos;
          if (ix != jx)
          {
            assert (dr.norm2() > 0.0);
            direct.push(dr*0.5, jx);
            assert(ptclP[jx].id >= 0);
            ptclP[jx].id = -1-ptclP[jx].id;
            plist.push(jx);
          }
        }
        //        tree.buildDirectPolyhedron(ptcl, group[i], direct, f);
#endif 

        std::vector<Particle> ptcl1(ptclP);
#if 1
        std::sort(ptcl1.begin(), ptcl1.end(), CmpDist(ipos));
        for (int j = 1; j < n_bodies; j++)
        {
          if (ptcl1[j].id < 0) continue;
          const vec3 dr = ptcl1[j].pos - ipos;
          assert (dr.norm2() > 0.0);
          direct.push(dr*0.5, j);
        }
#else
    //    std::random_shuffle(ptcl1.begin(), ptcl1.end());
        for (int j = 0; j < n_bodies; j++)
        {
          if (ptcl1[j].id < 0) continue;
          const vec3 dr = ptcl1[j].pos - ipos;
          if (dr.norm2() > 0.0)
            direct.push(dr*0.5,j);
        }
#endif
        while(!plist.empty())
        {
          const int jx = plist.top();
          plist.pop();
          ptclP[jx].id = -1-ptclP[jx].id;
        }

        nfaceD_min = std::min(nfaceD_min, direct.nface());
        nfaceD_max = std::max(nfaceD_max, direct.nface());
        nfaceD    += direct.nface();

        lp.clear();
        const int nf = direct.nface();
#ifdef REFLECTING
        lp.push(HalfSpace(vec3( 1.0, 0.0, 0.0), vec3(1.0*(rminD.x-ipos.x), 0.0, 0.0)));
        lp.push(HalfSpace(vec3(-1.0, 0.0, 0.0), vec3(1.0*(rmaxD.x-ipos.x), 0.0, 0.0)));
        lp.push(HalfSpace(vec3( 0.0, 1.0, 0.0), vec3(0.0, 1.0*(rminD.y-ipos.y), 0.0)));
        lp.push(HalfSpace(vec3( 0.0,-1.0, 0.0), vec3(0.0, 1.0*(rmaxD.y-ipos.y), 0.0)));
        lp.push(HalfSpace(vec3( 0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0*(rminD.z-ipos.z))));
        lp.push(HalfSpace(vec3( 0.0, 0.0,-1.0), vec3(0.0, 0.0, 1.0*(rmaxD.z-ipos.z))));
#endif
        for (int j = 0; j < nf; j++)
        {
          lp.push(direct.getHalfSpace(j));
          assert(direct[j] >= 0);
          used[direct[j]] = true;
        }

        list.clear();

#if 0
        const vec3 bnd_list[6] = 
        {
          vec3(2.0*(rmin.x-ipos.x), 0.0, 0.0),
          vec3(2.0*(rmax.x-ipos.x), 0.0, 0.0),
          vec3(0.0, 2.0*(rmin.y-ipos.y), 0.0),
          vec3(0.0, 2.0*(rmax.y-ipos.y), 0.0),
          vec3(0.0, 0.0, 2.0*(rmin.z-ipos.z)),
          vec3(0.0, 0.0, 2.0*(rmax.z-ipos.z))
        };
        for (int j = 0; j < 6; j++)
        {
          const vec3 &pos = bnd_list[j];
          const HalfSpace h(pos, pos*0.5);
          const vec3 p = lp.solve(h.n);
          if (!h.outside(p))
            list.push_back(Voronoi::Site(pos, -1-j));
        }
#endif

        for (int j = 0; j < n_bodies; j++)
        {
          const vec3 pos = ptcl1[j].pos - ipos;
          if (pos.norm2() == 0.0) continue;
          const HalfSpace h(pos, pos*0.5);

#if 1 
          if (used[j]) 
          {
            list.push_back(Voronoi::Site(pos, j));
            used[j] = false;
            continue;
          }
#endif

          const vec3 p = lp.solve(h.n);
          if (!h.outside(p))
            list.push_back(Voronoi::Site(pos, j));
        }
        nfaceL_min = std::min(nfaceL_min, (int)list.size());
        nfaceL_max = std::max(nfaceL_max, (int)list.size());
        nfaceL    += (int)list.size();
        fprintf(stderr, "np= %d: nlp= %d, nc= %d nf= %d\n",np, lp.n, (int)list.size(), nf);
        //        std::random_shuffle(list.begin(), list.end());


#ifdef REFLECTING
        list.push_back(Voronoi::Site(vec3(2.0*(rminD.x-ipos.x), 0.0, 0.0), -1));
        list.push_back(Voronoi::Site(vec3(2.0*(rmaxD.x-ipos.x), 0.0, 0.0), -2));
        list.push_back(Voronoi::Site(vec3(0.0, 2.0*(rmaxD.y-ipos.y), 0.0), -3));
        list.push_back(Voronoi::Site(vec3(0.0, 2.0*(rminD.y-ipos.y), 0.0), -4));
        list.push_back(Voronoi::Site(vec3(0.0, 0.0, 2.0*(rminD.z-ipos.z)), -5));
        list.push_back(Voronoi::Site(vec3(0.0, 0.0, 2.0*(rmaxD.z-ipos.z)), -6));
#endif

        std::random_shuffle(list.begin(), list.end());
        assert(cell.build(list));
        volume += cell.volume();

        nfaceV_min = std::min(nfaceV_min, (int)cell.nb());
        nfaceV_max = std::max(nfaceV_max, (int)cell.nb());
        nfaceV    += (int)cell.nb();
      }
    }
    delete cell_ptr;
  }
  assert(np == n_bodies0);
  const real lx = rmaxD.x - rminD.x;
  const real ly = rmaxD.y - rminD.y;
  const real lz = rmaxD.z - rminD.z;
  fprintf(stderr, " volume= %g   exact= %g  diff= %g \n",
      volume, lx*ly*lz, (volume-lx*ly*lz)/(lx*ly*lz));

#endif
  fprintf(stderr, " nfaceD: min= %d  max= %d  avg= %g\n",
      nfaceD_min, nfaceD_max, 1.0*nfaceD/n_bodies0);
  fprintf(stderr, " nfaceL: min= %d  max= %d  avg= %g\n",
      nfaceL_min, nfaceL_max, 1.0*nfaceL/n_bodies0);
  fprintf(stderr, " nfaceV: min= %d  max= %d  avg= %g\n",
      nfaceV_min, nfaceV_max, 1.0*nfaceV/n_bodies0);
  const double t1 = get_wtime();
  fprintf(stderr, " -- done in %g sec -- \n", t1 - t0);

  return 0;

}
