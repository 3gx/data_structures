#ifndef __OCTREE_H__
#define __OCTREE_H__

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <list>
#include "vector3.h"

typedef float real;
typedef vector3<real> vec3;

template<class T>
static inline T __min(const T &a, const T &b) {return a < b ? a : b;}
template<class T>
static inline T __max(const T &a, const T &b) {return a > b ? a : b;}


#define SQR(x) ((x)*(x))

struct Particle
{
  typedef std::vector<Particle> Vector;
  vec3 pos;
  int  id;

  Particle(const vec3 &_pos, const int _id) :
    pos(_pos), id(_id) {}
};



struct Octree
{
  enum {NLEAF = 16};
  enum {EMPTY = -1};
  enum {BODYX = -2};

  struct Body
  {
    typedef std::vector<Body> Vector;
    vec3 pos;
    int  idx;

    Body(const vec3 &_pos, const int _idx) : pos(_pos), idx(_idx) {}
    Body(const Particle &p, const int _idx) : pos(p.pos), idx(_idx) {}
  };

  vec3 root_centre;
  real root_size;
  int depth;
  int ncell;
  int n_nodes;
  std::vector<int> node_list;

  Octree(const vec3 &_centre, const real _size, const int _n_nodes) :
    root_centre(_centre), root_size(_size), depth(0), ncell(0), n_nodes(_n_nodes)
  {
    node_list.resize(n_nodes<<3);
    for (std::vector<int>::iterator it = node_list.begin(); it != node_list.end(); it++)
      *it = -1;
  }

  void clear()
  {
    depth = ncell = 0;
    for (std::vector<int>::iterator it = node_list.begin(); it != node_list.end(); it++)
      *it = -1;
  }

  static inline int reverse_int(const int a) {return BODYX-a;}

  static inline int Octant(const vec3 &lhs, const vec3 &rhs) 
  {
    return
      (((lhs.x <= rhs.x) ? 1 : 0) +  
       ((lhs.y <= rhs.y) ? 2 : 0) + 
       ((lhs.z <= rhs.z) ? 4 : 0));
  }

  static inline vec3 compute_centre(const vec3 &c, const real s, const int oct)
  {
    return vec3(
        c.x + s * ((oct&1) ? (real)1.0 : (real)-1.0),
        c.y + s * ((oct&2) ? (real)1.0 : (real)-1.0),
        c.z + s * ((oct&4) ? (real)1.0 : (real)-1.0)
        );
  }

  void check_overflow(const int n) const
  {
#if 1
    assert(n < (int)node_list.size());
#endif
  }

  int sanity_check() const
  {
    int nbody = 0;
    int nnode = 0;
    int empty = 0;
    for (int i = 0; i < (const int)node_list.size(); i++)
      if (node_list[i] == EMPTY)
        empty++;
      else if (node_list[i] >= 0)
        nnode++;
      else
        nbody++;
    fprintf(stderr, "empty= %d nnode= %d nbody= %d\n",
        empty, nnode, nbody);
    return nbody;
  }


  void push(const Body &body, const int idx, const Body::Vector &bodies)
  {
    int child_idx = 0;   /* child idx inside a node */
    int child     = 0;   /* child */
    int locked    = 0;   /* locked cell */
    int depth     = 0;   /* depth */ 

    vec3 centre = root_centre;
    real hsize  = root_size*(real)0.5;
    while (child > EMPTY)  /* if positive, this means node is already there */
    {
      const int node = child;
      depth++;

      child_idx  = Octant(centre, body.pos);
      hsize     *= (real)0.5;
      centre     = compute_centre(centre, hsize, child_idx);

      locked = node + child_idx;
      check_overflow(locked);
      child  = node_list[locked];
    }

    if (child == EMPTY)
      node_list[locked] = reverse_int(idx);
    else
    {
      while((child = node_list[locked]) != EMPTY)
      {
        assert(child < EMPTY);
        depth++;
        ncell++;

        const int cfirst = ncell<<3;
        check_overflow(cfirst+7);
        for (int k = 0; k < 8; k++)
          node_list[cfirst + k] = EMPTY;
        node_list[locked] = cfirst;

        const int body_id = reverse_int(child);
        assert(body_id >= 0);
        child_idx = Octant(centre, bodies[body_id].pos);
        node_list[cfirst + child_idx] = child;

        child_idx = Octant(centre, body.pos);
        hsize    *= (real)0.5;
        centre    = compute_centre(centre, hsize, child_idx);
        locked    = cfirst + child_idx;
      }
      node_list[locked] = reverse_int(idx);
    }

    this->depth = __max(this->depth, depth);
  }

  void morton_dump(std::vector<int> &list) const
  {
    for (int k = 0; k < 8; k++)
      morton_dump_recursive(k, list);
  }

  void morton_dump_recursive(const int node, std::vector<int> &list) const
  {
    assert(node < (int)node_list.size());
    const int cell = node_list[node];
    if (cell == EMPTY) return;
    if (cell >  EMPTY)
    {
      for (int k = 0; k < 8; k++)
        morton_dump_recursive(cell+k, list);
    }
    else
      list.push_back(reverse_int(cell));
  }
};


#endif /* __OCTREE_H__ */
