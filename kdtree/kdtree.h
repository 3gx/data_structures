#ifndef __KDTREE_H__
#define __KDTREE_H__

#include <cassert>
#include <vector>
#include <algorithm>
#include <stack>
#include <queue>
#include "vector3.h"

typedef float real;
typedef vector3<real> vec3;

template<class T>
static inline T __min(const T &a, const T &b) {return a < b ? a : b;}
template<class T>
static inline T __max(const T &a, const T &b) {return a > b ? a : b;}

  template <class T>
T prevPow2(const T v)
{
  int k;
  if (v == 0)
    return 1;
  for (k = sizeof(T) * 8 - 1; ((static_cast<T>(1U) << k) & v) == 0; k--);
  return static_cast<T>(1U) << k;
}

struct kdStack
{
  int beg, end;
  int node, depth;

  kdStack(const int _beg, const int _end, const int _node, const int _depth) :
    beg(_beg), end(_end), node(_node), depth(_depth) {}
  

#ifdef __mySSE__ 
  const kdStack operator=(const kdStack &rhs)
  {
    typedef float v4sf __attribute__ ((vector_size(16)));
    v4sf *lp =(v4sf *)this;
    v4sf *rp =(v4sf *)(&rhs);
    lp[0] = rp[0];
    return *this;
  }
#endif
};


class Particle
{
  public:
  typedef std::vector<Particle>         Vector;
  typedef Vector::iterator            Iterator;
  typedef Vector::const_iterator constIterator;

  public:

  vec3 pos;
  real mass;

  Particle(const vec3 &_pos, const real _mass) : pos(_pos), mass(_mass) {}
};


class kdBody
{
  public:
  typedef std::vector<kdBody>         Vector;
  typedef Vector::iterator            Iterator;
  typedef Vector::const_iterator constIterator;

  private:
  vec3 _pos;
  unsigned int _idx;

  public:

  kdBody(const Particle &ptcl, const int __idx) : _pos(ptcl.pos), _idx(__idx) 
  {
    assert(sizeof(kdBody) == sizeof(float)*4);
  }
  const vec3& pos() const {return _pos;}
  unsigned int idx() const {return _idx;}
   
#ifdef __mySSE__ 
  const kdBody operator = (const kdBody &rhs)
  {
    typedef float v4sf __attribute__ ((vector_size(16)));
    v4sf *lp =(v4sf *)this;
    v4sf *rp =(v4sf *)(&rhs);
    lp[0] = rp[0];
    return *this;
  }
#endif
};

class CompareBodies
{
  unsigned int _split_dim;
  public:
  CompareBodies(const unsigned int val) : _split_dim(val%3) {};
  bool operator() (const kdBody &lhs, const kdBody &rhs) const
  {
    return lhs.pos()[_split_dim] < rhs.pos()[_split_dim];
  }
  unsigned int split_dim() const {return _split_dim;}
};

#ifdef __mySSE__
namespace std
{
  template <> 
    inline void iter_swap <kdBody::Iterator, kdBody::Iterator> (kdBody::Iterator a, kdBody::Iterator b)
    {
      typedef float v4sf __attribute__ ((vector_size(16)));
      v4sf *ap =(v4sf *)&(*a);
      v4sf *bp =(v4sf *)&(*b);
      v4sf tmpa0 = ap[0];
      v4sf tmpb0 = bp[0];
      ap[0] = tmpb0;
      bp[0] = tmpa0;
    }
}
#endif

class kdNode
{
  unsigned int packed_data;
  vec3     _pos;

  public:

  typedef std::vector<kdNode>           Vector;
  typedef Vector::iterator            Iterator;
  typedef Vector::const_iterator constIterator;

  kdNode() : packed_data(0xFFFFFFFF) {}
  kdNode(const kdBody &body, const unsigned int split_dim):
    packed_data(body.idx() | (split_dim) << 29), _pos(body.pos()) 
  {
    assert(split_dim  < 3);
    assert(body.idx() < 0x20000000); /* only 536M bodies are supported */
  }

  int split_dim() const 
  {
    return packed_data >> 29;
  }
  int body_id() const
  {
    return packed_data & 0x1FFFFFFF;
  }

  const vec3 &pos() const {return _pos;}
};

class kdTree
{
  enum TREE_TYPE {BALLANCED, LEFT_BALLANCED, UNBALLANCED};
  std::vector<kdNode> nodes;
  int depth;

  public:

  kdTree(const Particle::Vector &ptcl, const TREE_TYPE type = LEFT_BALLANCED) :
    depth(0)

  {
    switch(type)
    {
      case BALLANCED:
        build_ballanced_tree(ptcl);
        break;
      case UNBALLANCED:
        build_unballanced_tree(ptcl);
        break;
      default: /* LEFT_BALLANCED */
        build_left_ballanced_tree(ptcl);
    }
  }

  private:

  void build_ballanced_tree(const Particle::Vector &ptcl)
  {
    assert(false);
  }

  void build_unballanced_tree(const Particle::Vector &ptcl)
  {
    assert(false);
  }


  void recursively_build_left_ballanced_tree(
      const int n_node,
      const int depth,
      kdBody::Iterator bodies_beg,
      kdBody::Iterator bodies_end)
  {
    const int n  = bodies_end - bodies_beg;
    if (n <= 0) return;

    const int m  = prevPow2(n);
    const int r  = n - (m - 1);
    const int lt = (m-2)/2 + (r <= m/2 ? r : m/2);

    const int median = lt;

    const int split_dim = depth%3;
    std::nth_element(bodies_beg, bodies_beg+median, bodies_end, CompareBodies(split_dim));
    this->depth = __max(this->depth, depth);

    const kdBody &body = *(bodies_beg + median);
    assert(n_node < (int)nodes.size());
    nodes[n_node] = kdNode(body, split_dim);

    recursively_build_left_ballanced_tree( n_node<<1,    depth+1, bodies_beg,          bodies_beg+median);
    recursively_build_left_ballanced_tree((n_node<<1)+1, depth+1, bodies_beg+median+1, bodies_end       );
  }
  
  void iteratively_build_left_ballanced_tree(kdBody::Vector &bodies)
  {
#if 0  /* depth first */
    std::stack<kdStack> stack;
    stack.push(kdStack(0, bodies.size(), 1, 0));

    while (!stack.empty())
    {
      const kdStack  s = stack.top();
      const int      n = s.end - s.beg;
      const int n_node = s.node; 
      const int  depth = s.depth;
      stack.pop();
    
      const int m  = prevPow2(n);
      const int r  = n - (m - 1);
      const int lt = (m-2)/2 + (r <= m/2 ? r : m/2);

      const int median = s.beg+lt;
      const int split_dim = depth%3;
      std::nth_element(bodies.begin()+s.beg, bodies.begin()+median, bodies.begin()+s.end, CompareBodies(split_dim));
      this->depth = __max(this->depth, depth);

      const kdBody &body = bodies[median];
      assert(n_node < (int)nodes.size());
      nodes[n_node] = kdNode(body, split_dim);

      const int n_left  = lt;
      const int n_right = n - lt - 1;
      if (n_left  > 0) stack.push(kdStack(s.beg,    median, (n_node<<1)+0, depth+1));
      if (n_right > 0) stack.push(kdStack(median+1, s.end,  (n_node<<1)+1, depth+1));
    }
#else  /* breadth first */
    std::vector<kdStack> *list      = new std::vector<kdStack>;
    std::vector<kdStack> *list_next = new std::vector<kdStack>;

    list_next->push_back(kdStack(0, bodies.size(), 1, 0));

    while (!list_next->empty())
    {
      std::swap(list, list_next);
      const int nlist = list->size();
      if (depth < 4)
        for (int i = 0; i < nlist; i++)
        {
          const kdStack &s = (*list)[i];
          const int      n = s.end - s.beg;
          const int n_node = s.node; 
          const int  depth = s.depth;

          const int m  = prevPow2(n);
          const int r  = n - (m - 1);
          const int lt = (m-2)/2 + (r <= m/2 ? r : m/2);

          const int median    =s.beg +  lt;
          const int split_dim = depth%3;
          std::nth_element(bodies.begin()+s.beg, bodies.begin()+median, bodies.begin()+s.end, CompareBodies(split_dim));
          this->depth = __max(this->depth, depth);

          const kdBody &body = bodies[median];
          assert(n_node < (int)nodes.size());
          nodes[n_node] = kdNode(body, split_dim);

          const int n_left  = lt;
          const int n_right = n - lt - 1;
          if (n_left  > 0) list_next->push_back(kdStack(s.beg,    median, (n_node<<1)+0, depth+1));
          if (n_right > 0) list_next->push_back(kdStack(median+1, s.end,  (n_node<<1)+1, depth+1));
        }
      else  
      {
#pragma omp parallel for
        for (int i = 0; i < nlist; i++)
        {
          const kdStack &s = (*list)[i];
          const int      n = s.end - s.beg;
          const int n_node = s.node; 
          const int  depth = s.depth;

          const int m  = prevPow2(n);
          const int r  = n - (m - 1);
          const int lt = (m-2)/2 + (r <= m/2 ? r : m/2);

#if 0
          const int median    = lt;
          const int split_dim = depth%3;
          kdBody::Vector lbodies(bodies.begin()+s.beg, bodies.begin() + s.end);
          std::nth_element(lbodies.begin(), lbodies.begin()+median, lbodies.end(), CompareBodies(split_dim));
          this->depth = __max(this->depth, depth);

          const kdBody &body = lbodies[median];
          assert(n_node < (int)nodes.size());
          nodes[n_node] = kdNode(body, split_dim);

          const int n_left  = lt;
          const int n_right = n - lt - 1;
#pragma omp critical
          {
            if (n_left  > 0) list_next->push_back(kdStack(s.beg,          s.beg+median, (n_node<<1)+0, depth+1));
            if (n_right > 0) list_next->push_back(kdStack(s.beg+median+1, s.end,        (n_node<<1)+1, depth+1));
          }
#else
          const int median    =s.beg +  lt;
          const int split_dim = depth%3;
          std::nth_element(bodies.begin()+s.beg, bodies.begin()+median, bodies.begin()+s.end, CompareBodies(split_dim));
          this->depth = __max(this->depth, depth);

          const kdBody &body = bodies[median];
          assert(n_node < (int)nodes.size());
          nodes[n_node] = kdNode(body, split_dim);

          const int n_left  = lt;
          const int n_right = n - lt - 1;
#pragma omp critical
          {
            if (n_left  > 0) list_next->push_back(kdStack(s.beg,    median, (n_node<<1)+0, depth+1));
            if (n_right > 0) list_next->push_back(kdStack(median+1, s.end,  (n_node<<1)+1, depth+1));
          }
#endif
        }
      }
      list->clear();
    }
    delete list;
    delete list_next;
#endif
  }

  void build_left_ballanced_tree(const Particle::Vector &ptcl)
  {
    kdBody::Vector bodies;
    bodies.reserve(ptcl.size());
    nodes.resize(ptcl.size()+2);
    for (unsigned int i = 0; i < ptcl.size(); i++)
      bodies.push_back(kdBody(ptcl[i], i));

#if 0
    recursively_build_left_ballanced_tree(1, 0, bodies.begin(), bodies.end());
#else
    iteratively_build_left_ballanced_tree(bodies);
#endif
  }



};

#endif /* __KDTREE_H__ */
