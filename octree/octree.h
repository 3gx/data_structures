#ifndef __OCTREE_H__
#define __OCTREE_H__

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <stack>
#include <list>
#include "vector3.h"
#include "boundary.h"

typedef float real;
typedef vector3<real> vec3;
typedef Boundary<real> boundary;

template<class T>
static inline T __min(const T &a, const T &b) {return a < b ? a : b;}
template<class T>
static inline T __max(const T &a, const T &b) {return a > b ? a : b;}


#define SQR(x) ((x)*(x))


#ifdef __mySSE__
typedef int    v4si  __attribute__((vector_size(16)));
typedef float  v4sf  __attribute__((vector_size(16)));
#endif

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
    private:
#ifdef __mySSE__
    v4sf packed_pos;
#else
    vec3 _pos;
    int  _idx;
#endif

    public:
#ifdef __mySSE__
    Body(const vec3 &pos, const int idx) 
    {
      packed_pos = (v4sf){pos.x, pos.y, pos.z, (float)idx};
    }
    Body(const Particle &p, const int idx) 
    {
      packed_pos = (v4sf){p.pos.x, p.pos.y, p.pos.z, (float)idx};
    }
    int idx() const {return (int)__builtin_ia32_vec_ext_v4sf(packed_pos, 3);}
    v4sf pos() const {return packed_pos;}
#else
    Body(const vec3 &__pos, const int __idx) : _pos(__pos), _idx(__idx) {}
    Body(const Particle &p, const int __idx) : _pos(p.pos), _idx(__idx) {}
    int idx() const {return _idx;}
    const vec3& pos() const {return _pos;}
#endif

  };

  vec3 root_centre;
  real root_size;
  int depth;
  int ncell;
  int n_nodes;
  std::vector<int> node_list;
  std::vector<boundary> innerBnd;

  Octree(const vec3 &_centre, const real _size, const int _n_nodes) :
    root_centre(_centre), root_size(_size), depth(0), ncell(0), n_nodes(_n_nodes)
  {
    node_list.resize(n_nodes<<3);
    innerBnd.resize(n_nodes<<3);
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

#ifdef __mySSE__
  static inline int Octant(const v4sf lhs, const v4sf rhs) 
  {
    int mask = __builtin_ia32_movmskps(
        __builtin_ia32_cmpgeps(rhs, lhs));
    return 7 & mask;
  }
#else
  static inline int Octant(const vec3 &lhs, const vec3 &rhs) 
  {
    return
      (((lhs.x <= rhs.x) ? 1 : 0) +  
       ((lhs.y <= rhs.y) ? 2 : 0) + 
       ((lhs.z <= rhs.z) ? 4 : 0));
  }
#endif

#ifdef __mySSE__
  static inline v4sf compute_centre(const v4sf centre, const real dummy, const int oct)
  {
    static const v4sf off[8] = {
      {-0.25f, -0.25f, -0.25f, -0.5f},
      {+0.25f, -0.25f, -0.25f, -0.5f},
      {-0.25f, +0.25f, -0.25f, -0.5f},
      {+0.25f, +0.25f, -0.25f, -0.5f},
      {-0.25f, -0.25f, +0.25f, -0.5f},
      {+0.25f, -0.25f, +0.25f, -0.5f},
      {-0.25f, +0.25f, +0.25f, -0.5f},
      {+0.25f, +0.25f, +0.25f, -0.5f},
    };
    const v4sf len = __builtin_ia32_shufps(centre, centre, 0xff);
    return centre + len * off[oct];
  }
#else
  static inline vec3 compute_centre(const vec3 &c, real s, const int oct)
  {
    assert(oct >= 0);
    assert(oct  < 8);
    s *= (real)0.5;
    return vec3(
        c.x + s * ((oct&1) ? (real)1.0 : (real)-1.0),
        c.y + s * ((oct&2) ? (real)1.0 : (real)-1.0),
        c.z + s * ((oct&4) ? (real)1.0 : (real)-1.0)
        );
  }
#endif

  int sanity_check(const Body::Vector &bodies) const
  {
    int nb = 0;
    for (int k = 0; k < 8; k++)
      if (node_list[k] != EMPTY)
      {
        vec3 c = root_centre;
        real h = root_size*(real)0.5;
        c = compute_centre(c, h, k);
        h *= (real)0.5;
        sanity_check_recursive(k, bodies, c, h, nb);
      }
    return nb;
  }

  void sanity_check_recursive(
      const int node, 
      const Body::Vector &bodies,
      const vec3 c, const real h, 
      int &nb) const
  {
    assert(node < (int)node_list.size());
    const int cell = node_list[node];
    assert(cell != EMPTY);

    if (cell > EMPTY)
    {
      for (int k = 0; k < 8; k++)
        if (node_list[cell+k] != EMPTY) 
        {
          const vec3 c1 = compute_centre(c, h, k);
          const real h1 = h * (real)0.5;
          sanity_check_recursive(cell+k, bodies, c1, h1, nb);
        }
    }
    else
    {
      const vec3 jpos = bodies[reverse_int(cell)].pos();
      const real h1 = h*1.001;
      const bool overlap_x = (jpos.x < c.x+h1) && (jpos.x > c.x-h1);
      const bool overlap_y = (jpos.y < c.y+h1) && (jpos.y > c.y-h1);
      const bool overlap_z = (jpos.z < c.z+h1) && (jpos.z > c.z-h1);
      const bool overlap   = overlap_x && overlap_y && overlap_z;
      assert(overlap_x);
      assert(overlap_y);
      assert(overlap_z);
      assert(overlap);
      nb++;
    }
  }
  
  int sanity_check1(const Body::Vector &bodies) const
  {
    int nb = 0;
    for (int k = 0; k < 8; k++)
      sanity_check_recursive(k, bodies, nb);
    return nb;
  }

  void sanity_check_recursive(
      const int node, 
      const Body::Vector &bodies,
      int &nb) const
  {
    assert(node < (int)node_list.size());
    const int cell = node_list[node];
    if (cell == EMPTY) return;

    if (cell > EMPTY)
    {
      for (int k = 0; k < 8; k++)
        if (node_list[cell+k] != EMPTY) 
          sanity_check_recursive(cell+k, bodies, nb);
    }
    else
    {
      const vec3 jpos = bodies[reverse_int(cell)].pos();
      assert(overlapped(innerBnd[node], jpos));
      nb++;
    }
  }

  /********/


  void push(const Body &body, const int idx, const Body::Vector &bodies)
  {
    int child_idx = 0;   /* child idx inside a node */
    int child     = 0;   /* child */
    int locked    = 0;   /* cell that needs to be updated */
    int depth     = 0;   /* depth */ 
    const int _n_nodes = node_list.size();

#ifdef __mySSE__
    v4sf centre = {root_centre.x, root_centre.y, root_centre.z, root_size};
#else
    vec3 centre = root_centre;
#endif
    real hsize  = root_size*(real)0.5; /* not being used in SSE version */

    /* walk the tree to find the first leaf or empty cell */
    while (child > EMPTY)  /* if non-negative, means it is a tree-node */
    {
      const int node = child;
      depth++;

      child_idx  = Octant(centre, body.pos());
      centre     = compute_centre(centre, hsize, child_idx);
      hsize     *= (real)0.5;

      locked = node + child_idx;
      assert(locked < _n_nodes);
      child  = node_list[locked];
    }

    /* locked on the cell that needs to be updated */

    if (child != EMPTY)
      while((child = node_list[locked]) != EMPTY)
      {
        assert(child < EMPTY);
        depth++;
        ncell++;

        const int cfirst = ncell<<3;
        assert(cfirst+7 < _n_nodes);
        node_list[locked] = cfirst;

        const int body_id = reverse_int(child);
        assert(body_id >= 0);
        child_idx = Octant(centre, bodies[body_id].pos());
        node_list[cfirst + child_idx] = child;

        child_idx = Octant(centre, body.pos());
        centre    = compute_centre(centre, hsize, child_idx);
        hsize    *= (real)0.5;
        locked    = cfirst + child_idx;
      }
    assert(node_list[locked] == EMPTY);
    node_list[locked] = reverse_int(idx);

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

#if 0
  int range_search(const vec3 &pos, const real s, const Body::Vector &bodies) const
  {
    int nb = 0;
#ifdef __mySSE__
#error "SSE version is not implemented yet"
    const v4sf ipos   = (v4sf){pos.x, pos.y, pos.z, 0.0f};
    const v4sf centre = (v4sf){root_centre.x, root_centre.y, root_centre.z, 0.5f*root_size};
    for (int k = 0; k < 8; k++)
      range_search_recursive(k, ipos, centre, nb);
#else
    for (int k = 0; k < 8; k++)
      if (node_list[k] != EMPTY)
      {
        vec3 c = root_centre;
        real h = root_size*(real)0.5;
        c = compute_centre(c, h, k);
        h *= (real)0.5;
        range_search_recursive(k, pos, s, bodies, c, h, nb);
      }
#endif
    return nb;
  }

  void range_search_recursive(
      const int node, 
      const vec3 &ipos,   const real s,
      const Body::Vector &bodies,
      vec3 c, real h, 
      int &nb) const
  {
    assert(node < (int)node_list.size());
    const int cell = node_list[node];
    assert(cell != EMPTY);
    const bool overlap_x = (ipos.x-s < c.x+h) && (ipos.x+s > c.x - h);
    const bool overlap_y = (ipos.y-s < c.y+h) && (ipos.y+s > c.y - h);
    const bool overlap_z = (ipos.z-s < c.z+h) && (ipos.z+s > c.z - h);
    const bool overlap   = overlap_x && overlap_y && overlap_z;
    if (!overlap) return;

    if (cell > EMPTY)
    {
      for (int k = 0; k < 8; k++)
        if (node_list[cell+k] != EMPTY) 
        {
          const vec3 c1 = compute_centre(c, h, k);
          const real h1 = h * (real)0.5;
          range_search_recursive(cell+k, ipos, s, bodies, c1, h1, nb);
        }
    }
    else
    {
      const vec3 jpos = bodies[reverse_int(cell)].pos();
      if ((ipos - jpos).norm2() < s*s)
        nb++;
    }
  }
#else
  struct walkStack
  {
    int  node;
    vec3 c;
    real h;
    walkStack(const int _node, const vec3 &_c, const real _h) :
      node(_node), c(_c), h(_h) {}
  };
  int range_search(const vec3 &ipos, const real s, const Body::Vector &bodies) const
  {
    int nb = 0;
    std::stack<walkStack> stack;

    for (int k = 0; k < 8; k++)
      if (node_list[k] != EMPTY)
      {
        vec3 c = root_centre;
        real h = root_size*(real)0.5;
        c = compute_centre(c, h, k);
        h *= (real)0.5;
        stack.push(walkStack(k, c, h));
      }

    while(!stack.empty())
    {
      const walkStack top = stack.top();
      stack.pop();

      const int  node = top.node;
      const vec3   &c = top.c;
      const real   &h = top.h;

      assert(node < (int)node_list.size());
      const int cell = node_list[node];
      assert(cell != EMPTY);
      const bool overlap_x = (ipos.x-s < c.x+h) && (ipos.x+s > c.x - h);
      const bool overlap_y = (ipos.y-s < c.y+h) && (ipos.y+s > c.y - h);
      const bool overlap_z = (ipos.z-s < c.z+h) && (ipos.z+s > c.z - h);
      const bool overlap   = overlap_x && overlap_y && overlap_z;
      if (!overlap) continue;

      if (cell > EMPTY)
      {
        for (int k = 0; k < 8; k++)
          if (node_list[cell+k] != EMPTY) 
          {
            const vec3 c1 = compute_centre(c, h, k);
            const real h1 = h * (real)0.5;
            stack.push(walkStack(cell+k, c1, h1));
          }
      }
      else
      {
        const vec3 jpos = bodies[reverse_int(cell)].pos();
        if ((ipos - jpos).norm2() < s*s)
          nb++;
      }
    }

    return nb;
  }
#endif

#if 1
  int range_search1(const vec3 &pos, const real h, const Body::Vector &bodies) const
  {
    int nb = 0;
    for (int k = 0; k < 8; k++)
      range_search_recursive(k, pos, h*h, boundary(pos, h), bodies, nb);
    return nb;
  }

  void range_search_recursive(
      const int node, 
      const vec3 &ipos, const real h2, const boundary &ibnd,
      const Body::Vector &bodies,
      int &nb) const
  {
    assert(node < (int)node_list.size());
    const int cell = node_list[node];
    if (cell == EMPTY) return;
    if (not_overlapped(ibnd, innerBnd[node])) return;

    if (cell > EMPTY)
    {
      for (int k = 0; k < 8; k++)
        range_search_recursive(cell+k, ipos, h2, ibnd, bodies, nb);
    }
    else
    {
      const vec3 jpos = bodies[reverse_int(cell)].pos();
      if ((ipos - jpos).norm2() < h2)
        nb++;
    }
  }
#else
  int range_search1(const vec3 &ipos, const real h, const Body::Vector &bodies) const
  {
    int nb = 0;
    std::stack<int> stack;
    const boundary ibnd(ipos, h);
    const real h2 = h*h;

    for (int k = 0; k < 8; k++)
      if (node_list[k] != EMPTY)
        stack.push(k);

    while(!stack.empty())
    {
      const int node = stack.top();
      stack.pop();

      assert(node < (int)node_list.size());
      const int cell = node_list[node];
      assert(cell != EMPTY);
      if (not_overlapped(ibnd, innerBnd[node])) continue;

      if (cell > EMPTY)
      {
        for (int k = 0; k < 8; k++)
          if (node_list[cell+k] != EMPTY)
            stack.push(cell+k);
      }
      else
      {
        const vec3 jpos = bodies[reverse_int(cell)].pos();
        if ((ipos - jpos).norm2() < h2)
          nb++;
      }
    }

    return nb;
  }
#endif

  boundary inner_boundary(const Body::Vector &bodies)
  {
    boundary bnd;
    for (int k = 0; k < 8; k++)
    {
      innerBnd[k] = inner_boundary_recursive(k, bodies);
      bnd.merge(innerBnd[k]);
    }
    return bnd;
  }
  boundary inner_boundary_recursive(const int node, const Body::Vector &bodies) 
  {
    assert(node < (int)node_list.size());
    const int cell = node_list[node];
    if (cell == EMPTY) return boundary();

    innerBnd[node] = boundary();

    if (cell > EMPTY)
    {
      for (int k = 0; k < 8; k++)
        if (node_list[cell+k] != EMPTY) 
          innerBnd[node].merge(inner_boundary_recursive(cell+k, bodies));
    }
    else
    {
      const vec3 jpos = bodies[reverse_int(cell)].pos();
      innerBnd[node] = boundary(jpos);
    }

    return innerBnd[node];
  }

};


#endif /* __OCTREE_H__ */
