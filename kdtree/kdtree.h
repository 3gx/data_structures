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

  template <class T>
T prevPow2(const T v)
{
  int k;
  if (v == 0)
    return 1;
  for (k = sizeof(T) * 8 - 1; ((static_cast<T>(1U) << k) & v) == 0; k--);
  return static_cast<T>(1U) << k;
}

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
};



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
};

template<const int NLEAF>
class kdLeaf
{
  public:
    typedef std::vector<kdLeaf> Vector;

  private:
    kdBody data[NLEAF];
    int n, iPad[3];

  public:
    kdLeaf() : n(0) {}
    kdLeaf(const kdBody::Vector &bodies) : n(bodies.size())
    {
      assert(n <= NLEAF);
      for (int i = 0; i < n; i++)
        data[i] = bodies[i];
    }
    void insert(const kdBody &body) {assert(n<NLEAF); data[n++] = body;}
    bool isEmpty() const {return  n == 0;}
    bool isFull()  const {return n == NLEAF;}
    const kdBody& operator[](const int i) const { return data[i]; }
};

class kdTree
{
  public:
    enum {NLEAF=16};
  private:
    kdNode       ::Vector nodes;
    kdLeaf<NLEAF>::Vector leaves;

  int depth;

  public:

  const kdNode& operator[](const int i) const { return nodes[i];}


  kdTree(const Particle::Vector &ptcl) : depth(0)
  {
    build_left_ballanced_tree(ptcl);
  }

  private:
  
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

  void build_left_ballanced_tree(const Particle::Vector &ptcl)
  {
    kdBody::Vector bodies;
    bodies.reserve(ptcl.size());
    nodes.resize(ptcl.size()+2);
    for (unsigned int i = 0; i < ptcl.size(); i++)
      bodies.push_back(kdBody(ptcl[i], i));

    recursively_build_left_ballanced_tree(1, 0, bodies.begin(), bodies.end());
  }


  public:
  int find_nnb(const vec3 &pos) const
  {
    real s2min     = HUGE;
    int  inode_min = -1;
    find_recursively_nnb(pos, 1, s2min, inode_min);

    return inode_min;
  }

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
    
    const real s2 = (pos - node.pos()).norm2();
    if (s2 <= s2min && s2 > 0.0)
    {
      s2min      = s2;
      inode_min = inode;
    }

    find_recursively_nnb(pos, near, s2min, inode_min);
    if (SQR(pos[split_dim] - node.pos()[split_dim]) <= s2min)
      find_recursively_nnb(pos, far, s2min, inode_min);

  }

};

#endif /* __KDTREE_H__ */
