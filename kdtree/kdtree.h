#ifndef __KDTREE_H__
#define __KDTREE_H__

#include <cassert>
#include <vector>
#include <algorithm>
#include <stack>
#include <parallel/partition.h>
#include "vector3.h"

typedef float real;
typedef vector3<real> vec3;

template<class T>
static inline T __min(const T &a, const T &b) {return a < b ? a : b;}
template<class T>
static inline T __max(const T &a, const T &b) {return a > b ? a : b;}


#define SQR(x) ((x)*(x))

#if 1
  template <class T>
T prevPow2(const T v)
{
  int k;
  if (v == 0)
    return 1;
  for (k = sizeof(T) * 8 - 1; ((static_cast<T>(1U) << k) & v) == 0; k--);
  return static_cast<T>(1U) << k;
}
#else
template <class T>
T prevPow2(const T v) {
  if (v == 0)
    return 1;

}
#endif

struct MinDist
{
  int   id;
  float s2;

  MinDist() : id(-1), s2(HUGE) {}
  MinDist(const int _id, const float _s2) : id(_id), s2(_s2) {}

  bool operator<(const MinDist &rhs) const
  {
    return s2 < rhs.s2;
  }
  bool operator<=(const MinDist &rhs) const
  {
    return s2 <= rhs.s2;
  }
  bool operator>(const MinDist &rhs) const
  {
    return s2 > rhs.s2;
  }
  bool operator()(const MinDist &lhs, const MinDist &rhs) const
  {
    return lhs < rhs;
  }

};

template<class T, const int K>
class SortedList
{
  T list[K+1];

  public:

  SortedList()
  {
    for (int i = 0; i < K; i++)
      list[i] = T();
  }

  const T& operator[](const int i) const {return list[i];}

  bool push(const T &x)
  {
    if (x > list[K-1]) return false;
    int lo = 0;
    int hi = K;
    while (lo < hi) 
    {
      const int mid = (lo + hi) >> 1;
      if (list[mid] <= x) lo = mid + 1;
      else                hi = mid;
    }
    assert(lo <= K);

    for (int j = K-1; j >= lo; j--)
      list[j] = list[j-1];
    list[lo] = x;
    return true;
  }

  const T& back() const {return list[K-1];}
};

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
  int body_idx() const
  {
    return packed_data & 0x1FFFFFFF;
  }

  const vec3& pos() const {return _pos;}
 
#ifdef __mySSE__ 
  const kdNode operator=(const kdNode &rhs)
  {
    typedef float v4sf __attribute__ ((vector_size(16)));
    v4sf *lp =(v4sf *)this;
    v4sf *rp =(v4sf *)(&rhs);
    lp[0] = rp[0];
    return *this;
  }
#endif
};

class kdTree
{
  enum TREE_TYPE {BALLANCED, LEFT_BALLANCED, UNBALLANCED};
  std::vector<kdNode> nodes;
  int depth;

  public:

  const kdNode& operator[](const int i) const { return nodes[i];}


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

  void nth_element_omp(
      kdBody::Iterator beg, 
      kdBody::Iterator med, 
      kdBody::Iterator end,
      const int split_dim)
  {
#ifdef _OPENMP
    __gnu_parallel::__parallel_nth_element(beg, med, end, CompareBodies(split_dim));
#else
    std::nth_element(beg, med, end, CompareBodies(split_dim));
#endif
  }

  class CompareBodies
  {
    unsigned int _split_dim;
    public:
    CompareBodies(const unsigned int val) : _split_dim(val%3) {};
    bool operator() (const kdBody &lhs, const kdBody &rhs) const
    {
      return lhs.pos()[_split_dim] < rhs.pos()[_split_dim];
    }
  };


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

#define OMP_BODIES_MIN 65536
 
#if 0 
  void sanity_check()
  {
    std::stack<int> stack;
    stack.push(1);

    vec3 pos(HUGE);

    while(!stack.empty())
    {
      const int node_id = stack.top();
      stack.pop();
      
      const kdNode &node = nodes[node_id];
      const int left  = node << 1
      const int right = left  + 1;



    }
  }

#endif


  void iteratively_build_left_ballanced_tree(kdBody::Vector &bodies)
  {
#if 1  /* depth first */
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
#if 1
      if (n > OMP_BODIES_MIN)
        nth_element_omp(bodies.begin()+s.beg, bodies.begin()+median, bodies.begin()+s.end, split_dim);
      else
        std::nth_element(bodies.begin()+s.beg, bodies.begin()+median, bodies.begin()+s.end, CompareBodies(split_dim));
#else
        std::sort(bodies.begin()+s.beg, bodies.begin()+s.end, CompareBodies(split_dim));
#endif
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
          nth_element_omp(bodies.begin()+s.beg, bodies.begin()+median, bodies.begin()+s.end, split_dim);
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

#if 1
    recursively_build_left_ballanced_tree(1, 0, bodies.begin(), bodies.end());
#else
    iteratively_build_left_ballanced_tree(bodies);
#endif
  }


  public:
  int find_nnb(const vec3 &pos) const
  {
    real s2min     = HUGE;
    int  inode_min = -1;
    find_recursively_nnb(pos, 1, s2min, inode_min);

    return inode_min;
  }
  template<const int K>
  void find_knb(const vec3 &pos, int knb_list[K]) const
  {
    SortedList<MinDist, K> list;
    find_recursively_knb(pos, 1, list);
    for (int k = 0; k < K; k++)
      knb_list[k] = list[k].id;
  }
#if 1
  int find_range_nb(const vec3 &pos, const real s) const
  {
    int nb = 0;
    vec3 rmin(-HUGE);
    vec3 rmax(+HUGE);
    find_recursively_range_nb(pos, 1, s, nb, rmin, rmax);
    return nb;
  }
#else
  int find_range_nb(const vec3 &pos, const real s) const
  {
    const int n_nodes = nodes.size();
    std::vector<int> kd_stack(n_nodes);

    int n_stack = 0;
    kd_stack[n_stack] = 1;

    int nb = 0;
    while (n_stack > -1)
    {
      int node_id = kd_stack[n_stack--];
      const kdNode &node = nodes[node_id];

      const int split = node.split_dim();
      node_id *= 2;
      if (node_id <= n_nodes)
      {
        if (pos[split] - s < node.pos()[split]) kd_stack[++n_stack] = node_id;
        if (pos[split] + s > node.pos()[split]) kd_stack[++n_stack] = node_id + 1;
      }

#if 0
      if ((pos - node.pos()).norm2() < s*s)
#endif
        nb++;
      
    }
    return nb;
  }
#endif
  private:

  void find_recursively_nnb(
      const vec3 &pos,
      const int   inode,
      real &s2min,
      int  &inode_min) const
  {
    if (inode >= (int)nodes.size()) return;

    const kdNode &node = nodes[inode];
    const int split_dim = node.split_dim();
    const bool go_left = pos[split_dim] <= node.pos()[split_dim];
    const int left  = (inode << 1);
    const int right = (inode << 1) + 1;
    const int near = go_left ? left  : right;
    const int far  = go_left ? right : left;

    find_recursively_nnb(pos, near, s2min, inode_min);
    if (SQR(pos[split_dim] - node.pos()[split_dim]) <= s2min)
      find_recursively_nnb(pos, far, s2min, inode_min);

    const real s2 = (pos - node.pos()).norm2();
    if (s2 <= s2min && s2 > 0.0)
    {
      s2min      = s2;
      inode_min = inode;
    }
  }

  template<const int K>
    void find_recursively_knb(
        const vec3 &pos,
        const int inode,
        SortedList<MinDist, K> &list) const
    {
      if (inode >= (int)nodes.size()) return;

      const kdNode &node  = nodes[inode];
      const int split_dim = node.split_dim();
      const bool go_left  = pos[split_dim] <= node.pos()[split_dim];
      const int left  = (inode << 1);
      const int right = (inode << 1) + 1;
      const int near = go_left ? left  : right;
      const int far  = go_left ? right : left;

      find_recursively_knb<K>(pos, near, list);
      if (SQR(pos[split_dim] - node.pos()[split_dim]) <= list.back().s2)
        find_recursively_knb<K>(pos, far, list);

      list.push(MinDist(inode, (pos - node.pos()).norm2()));
    }
  
    void find_recursively_range_nb(
        const vec3 &pos,
        const int inode,
        const real s,
        int &nb,
        vec3 &rmin,
        vec3 &rmax) const
    {
      if (inode >= (int)nodes.size()) return;

      const kdNode &node  = nodes[inode];
      const int split_dim = node.split_dim();
      const int left  = inode << 1;
      const int right = left + 1; 

#if 1
      if (pos[split_dim] - s < node.pos()[split_dim]) 
        find_recursively_range_nb(pos, left, s, nb, rmin, rmax);
      if (pos[split_dim] + s > node.pos()[split_dim]) 
        find_recursively_range_nb(pos, right, s, nb, rmin, rmax);
#else
      if (pos[split_dim] - s < node.pos()[split_dim]) 
      {
        rmax[split_dim] = node.pos()[split_dim];
        find_recursively_range_nb(pos, left, s, nb, rmin, rmax);
      }
#endif

#if 0
      if ((pos - node.pos()).norm2() < s*s)
#endif
        nb++;
    }
};

#endif /* __KDTREE_H__ */
