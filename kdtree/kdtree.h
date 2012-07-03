#ifndef __KDTREE_H__
#define __KDTREE_H__

#include <cassert>
#include <vector>
#include <algorithm>
#include <stack>
//#include <parallel/partition.h>
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
  real s2;

  MinDist() : id(-1), s2(HUGE) {}
  MinDist(const int _id, const real _s2) : id(_id), s2(_s2) {}

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

class Particle
{
  public:
    typedef std::vector<Particle>         Vector;
    typedef Vector::iterator            Iterator;
    typedef Vector::const_iterator constIterator;

  public:

    vec3 pos;
    real mass;

    Particle() {}
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

    kdBody() {}
    kdBody(const Particle &ptcl, const int __idx) : _pos(ptcl.pos), _idx(__idx)  {}
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
  kdNode(const int leaf) {packed_data = -1-leaf;}
  kdNode(const kdBody &body, const unsigned int split_dim):
    packed_data(body.idx() | (split_dim) << 28), _pos(body.pos()) 
  {
    assert(split_dim  < 3);
    assert(body.idx() < 0x10000000); /* only 536M bodies are supported */
    assert(!isLeaf());
  }

  int getLeaf() const {return -1-packed_data;}
  bool isLeaf() const {return packed_data & 0x80000000;}

  int split_dim() const 
  {
    assert(!isLeaf());
    return packed_data >> 28;
  }
  int body_idx() const
  {
    assert(!isLeaf());
    return packed_data & 0xFFFFFFF;
  }

  const vec3& pos() const {return _pos;}
};

template<const int NLEAF>
class kdLeafT
{
  public:
    typedef std::vector<kdLeafT> Vector;

  private:
    kdBody data[NLEAF];
    int n, iPad[3];

  public:
    kdLeafT() : n(0) {}
    kdLeafT(const kdBody::Vector &bodies) : n(bodies.size())
    {
      assert(n <= NLEAF);
      for (int i = 0; i < n; i++)
        data[i] = bodies[i];
    }
    kdLeafT(const kdBody::constIterator beg, const kdBody::constIterator end) : n(end-beg)
    {
      assert(n <= NLEAF);
      int i = 0;
      for (kdBody::constIterator it = beg; it != end; it++, i++)
        data[i] = *it;
    }
    void insert(const kdBody &body) {assert(n<NLEAF); data[n++] = body;}
    bool isEmpty() const {return  n == 0;}
    bool isFull()  const {return n == NLEAF;}
    const kdBody& operator[](const int i) const { return data[i]; }
    int size() const {return n;}
};

class kdTree
{
  public:
    enum {NLEAF=16};
    typedef kdLeafT<NLEAF> Leaf;
  private:
    kdNode::Vector nodes;
      Leaf::Vector leaves;

  int depth;

  public:

  const kdNode& operator[](const int i) const { return nodes[i];}

  kdTree(const Particle::Vector &ptcl) : depth(0)
  {
    build_left_ballanced_tree(ptcl);
  }
  int getDepth() const {return depth;}

  int nLeaf() const { return  leaves.size(); }
  const Leaf& getLeaf(const int i) const {return leaves[i];}

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

    if (n <= NLEAF)
    {
      leaves.push_back(Leaf(bodies_beg, bodies_end));
      nodes[n_node] = kdNode(leaves.size() - 1);
      return;
    }

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

  void dump(std::vector<int> &bodies, const int inode = 1) const
  {
    if (inode >= (int)nodes.size()-1) return;
    
    const kdNode &node = nodes[inode];
    if (node.isLeaf())
    {
      const Leaf &leaf = leaves[node.getLeaf()];
      for (int i = 0; i < leaf.size(); i++)
        bodies.push_back(leaf[i].idx());
    }
    else
    {
      bodies.push_back(node.body_idx());
      const int left  = inode << 1;
      const int right = left + 1;
      dump(bodies, left);
      dump(bodies, right);
    }
  } 

  int find_nnb(const vec3 &pos) const
  {
    real smin = HUGE;
    int  body = -1;
    find_recursively_nnb(pos, 1, smin, body);

    return body;
  }

  int find_nnb_inner(const vec3 n, const real h) const
  {
    real smin = HUGE;
    int  body = -1;

    vec3 min(-HUGE);
    vec3 max(+HUGE);
    /* revert normal, so that we can reuse outer half-space code */
    /* WARNING: the equation of plane that is passed here is n.r-h = 0 */
    /* while the recursive walks uses n.r+h = 0 for convenience */
    /* so, only 'n' changes sign */
    find_recursively_nnb_outer(n*(-1.0), h, 1, min, max, smin, body);
    return body;
  }
  private:

  void find_recursively_nnb(
      const vec3 &pos,
      const int   inode,
      real &smin,
      int  &body) const
  {
    if (inode >= (int)nodes.size()-1) return;

    const kdNode &node = nodes[inode];
    if (node.isLeaf())
    {
      const Leaf &leaf = leaves[node.getLeaf()];
      real s2min = SQR(smin);
      for (int i = 0; i < leaf.size(); i++)
      {
        const real s2 = (pos - leaf[i].pos()).norm2();
        if (s2 <= s2min && s2 > 0.0)
        {
          s2min = s2;
          body  = leaf[i].idx();
        }
      }
      smin = std::sqrt(s2min);
    }
    else
    {
      const int split_dim = node.split_dim();
      const real dist     = pos[split_dim] - node.pos()[split_dim];
      const bool go_left  = dist <= 0.0;
      const int left  = inode << 1;
      const int right = left + 1;

      const real s2 = (pos - node.pos()).norm2();
      if (s2 <= SQR(smin) && s2 > 0.0)
      {
        smin = std::sqrt(s2);
        body = node.body_idx();
      }

      const int near = go_left ? left  : right;
      const int far  = go_left ? right : left;
      find_recursively_nnb(pos, near, smin, body);
      if (std::abs(dist) <  smin)
        find_recursively_nnb(pos, far, smin, body);
    }

  }

  bool split_node_outer(
      const vec3 &min, const vec3 &max,
      const vec3 &n, const real h,
      const real smin) const
  {
    const vec3 v0(min[0], min[1], min[2]);
    const vec3 v1(min[0], min[1], max[2]);
    const vec3 v2(min[0], max[1], min[2]);
    const vec3 v3(min[0], max[1], max[2]);
    const vec3 v4(max[0], min[1], min[2]);
    const vec3 v5(max[0], min[1], max[2]);
    const vec3 v6(max[0], max[1], min[2]);
    const vec3 v7(max[0], max[1], max[2]);
    const real d0 = n*v0 + h;
    const real d1 = n*v1 + h;
    const real d2 = n*v2 + h;
    const real d3 = n*v3 + h;
    const real d4 = n*v4 + h;
    const real d5 = n*v5 + h;
    const real d6 = n*v6 + h;
    const real d7 = n*v7 + h;

    const real dmin = __min(
      __min(__min(d0, d1), __min(d2, d3)),
      __min(__min(d4, d5), __min(d6, d7)));

    const real dmax = __max(
      __max(__max(d0, d1), __max(d2, d3)),
      __max(__max(d4, d5), __max(d6, d7)));
   
    if (dmax > 0.0 && dmin < smin) return true;

    return false;
  }
 
  void find_recursively_nnb_outer(
      const vec3 &n, const real h,
      const int inode,
      const vec3 min, 
      const vec3 max,
      real &smin,
      int  &body) const
  {
    if (inode >= (int)nodes.size()-1) return;

    const kdNode &node = nodes[inode];
    if (node.isLeaf())
    {
      const Leaf &leaf = leaves[node.getLeaf()];
      for (int i = 0; i < leaf.size(); i++)
      {
        const real d = n*leaf[i].pos() + h;
        if (d > 0.0 && d < smin)
        {
          smin = d;
          body = leaf[i].idx();
        }
      }
    }
    else
    {
      const real d = n*node.pos() + h;
      if (d > 0.0 && d < smin)
      {
        smin = d;
        body = node.body_idx();
      }

      const int split_dim = node.split_dim();
      const int left  = inode << 1;
      const int right = left + 1;

      const real pnt = node.pos()[split_dim];
      vec3 lmax(max);
      lmax[split_dim] = pnt;
      if (split_node_outer(min, lmax, n, h, smin))
        find_recursively_nnb_outer(n, h, left, min, lmax, smin, body);

      vec3 rmin(min);
      rmin[split_dim] = pnt;
      if (split_node_outer(rmin, max, n, h, smin))
        find_recursively_nnb_outer(n, h, right, rmin, max, smin, body);
    }

  }

};

#endif /* __KDTREE_H__ */
