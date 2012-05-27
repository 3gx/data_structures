#ifndef __OCTREE_H__
#define __OCTREE_H__

#if 0
#define MEMALIGN
#endif

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <stack>
#ifdef MEMALIGN
#include "memalign_allocator.h"
#endif

template<class T>
static inline T __min(const T &a, const T &b) {return a < b ? a : b;}
template<class T>
static inline T __max(const T &a, const T &b) {return a > b ? a : b;}


#define SQR(x) ((x)*(x))

#ifdef __mySSE1__
#ifndef __mySSE__
#define __mySSE__
#endif
#endif


#ifdef __mySSE__
typedef int    v4si  __attribute__((vector_size(16)));
typedef float  v4sf  __attribute__((vector_size(16)));
#endif


#include "boundary.h"
#include "vector3.h"

typedef float real;
typedef vector3<real> vec3;
typedef Boundary<real> boundary;


struct Particle
{
#ifdef MEMALIGN
  typedef std::vector<Particle, __gnu_cxx::malloc_allocator<Particle, 64> > Vector;
#else
  typedef std::vector<Particle> Vector;
#endif
  vec3 pos;
  int  id;

  Particle(const vec3 &_pos, const int _id) :
    pos(_pos), id(_id) {}
};


struct Octree
{
  enum {NLEAF = 64};
  enum {EMPTY = -1};
  enum {BODYX = -2};

  struct Body
  {
#ifdef MEMALIGN
    typedef std::vector<Body, __gnu_cxx::malloc_allocator<Body, 64> > Vector;
#else
    typedef std::vector<Body> Vector;
#endif
    private:
    float4 packed_pos;

    public:
    Body() {}
    Body(const vec3 &pos, const int idx) : packed_pos(float4(pos.x, pos.y, pos.z, (float)idx)) {}
    Body(const Particle &p, const int idx) : packed_pos(float4(p.pos.x, p.pos.y, p.pos.z, (float)idx)) {}
    Body(const vec3 &pos, const float h) : packed_pos(pos.x, pos.y, pos.z, h) {};
    int idx() const {return (int)packed_pos.w();}
    float4 pos() const {return packed_pos;}
  };

  struct Cell
  {
    enum CellType {ROOT};
#ifdef MEMALIGN
    typedef std::vector<Cell, __gnu_cxx::malloc_allocator<Cell, 64> > Vector;
#else
    typedef std::vector<Cell> Vector;
#endif
    private:
    int _addr;
    int _id;

    public:
    Cell() : _addr(EMPTY), _id(EMPTY) {}
    Cell(const CellType type) : _addr(EMPTY+1), _id(EMPTY) {}
    Cell(const int addr, const int id) : _addr(addr), _id(id) {setTouched();}
    bool operator==(const Cell &v) const {return _addr == v._addr && _id == v._id;}
    bool isClean() const { return _addr == EMPTY && _id == EMPTY;}
    bool isEmpty() const { return _addr == EMPTY;}
    bool isLeaf () const { return _addr < EMPTY;}
    bool isNode () const { return _addr > EMPTY;}

    void set_addr(const int addr) {_addr = addr;}
    int addr() const {return _addr;}
    int id() const {return _id & 0x7FFFFFFF;}
    bool    isTouched() const {return _id & 0x80000000;}
    void   setTouched()       { _id |= 0x80000000;}
    void unsetTouched()       { _id = id();}

    int leafIdx() const 
    { 
      assert(_addr < EMPTY); 
      return reverse_int(_addr); 
    }
  };

  struct Leaf
  {
#ifdef MEMALIGN
    typedef std::vector<Leaf, __gnu_cxx::malloc_allocator<Leaf, 128> > Vector;
#else
    typedef std::vector<Leaf> Vector;
#endif
    private:
      int _nb;           /* number of bodies in the list */
      Body _list[NLEAF];  /* idx of the body */

    public:
      Leaf() : _nb(0) {}
      Leaf(const Body &body) : _nb(0) {insert(body);}
      void insert(const Body &body) 
      {
        assert(_nb < NLEAF);
        _list[_nb++] = body;
      }
      void remove(const int idx)
      {
        assert(idx < _nb);
        _nb--;
        for (int i = idx; i < _nb; i++)
          _list[i] = _list[i+1];
      }
      bool  isFull() const {return NLEAF == _nb;}
      bool isEmpty() const {return     0 == _nb;}
      const Body& operator[](const int i) const { return _list[i];}
      int nb() const {return _nb;}
  };



  private:
  float4 root_centre;  /* root's geometric centre & size */
  int depth;         /* tree-depth */
  int nnode;         /* number of tree nodes */
  int ncell;         /* number of tree cells: node + leaf */
  bool treeReady;    /* this is a flag if tree is ready for walks */
  Cell    ::Vector cellList;    /* list of tree cells (node+leaf) */
  Leaf    ::Vector leafList;    /* list of leaves */
  std::stack<int>  leafPool;    /* pool of available leaves */
  boundary::Vector innerBnd;    /* list of cell inner boundaries */
  boundary::Vector outerBnd;    /* list of cell outer boundaries */

  public:

  Octree(const vec3 &_centre, const real _size, const int n_nodes) :
    depth(0), nnode(0), ncell(0), treeReady(false)
  {
    root_centre = float4(_centre.x, _centre.y, _centre.z, _size);
    cellList.resize(n_nodes<<3);
  }

  int get_depth() const { return depth;}
  int get_nnode() const { return nnode;}
  int get_nleaf() const { return leafList.size();}
  int get_ncell() const { return ncell;}
  const vec3 get_rootCentre() const { return vec3(root_centre.x(), root_centre.y(), root_centre.z());}
  real get_rootSize() const { return root_centre.w();}


  void clear(const int n_nodes = 0)
  {
    depth = nnode = ncell = 0;
    treeReady = false;
    const int nsize = n_nodes <= 0 ? cellList.size() : n_nodes << 3;
    cellList.clear();
    cellList.resize(nsize);
    leafList.clear();
    leafPool = std::stack<int>();
    innerBnd.clear();
  }

  bool isTreeReady() const {return treeReady;}

  static inline int reverse_int(const int a) {return BODYX-a;}

  inline int Octant(const float4 lhs, const float4 rhs) 
  {
#ifdef __mySSE1__
    int mask = __builtin_ia32_movmskps(
        __builtin_ia32_cmpgeps(rhs, lhs));
    return 7 & mask;
#else
    return
      (((lhs.x() <= rhs.x()) ? 1 : 0) +  
       ((lhs.y() <= rhs.y()) ? 2 : 0) + 
       ((lhs.z() <= rhs.z()) ? 4 : 0));
#endif
  }

  inline float4 child_centre(const float4 centre, const float4 ppos)
  {
#ifdef __mySSE1__
    const v4si mask = {(int)0x80000000, (int)0x80000000, (int)0x80000000, 0};
    const v4sf off  = {0.25f, 0.25f, 0.25f, 0.5f};
    const v4sf len  = __builtin_ia32_shufps(centre, centre, 0xff);

    v4sf tmp = __builtin_ia32_cmpgeps(ppos, centre); // mask bits
    tmp = __builtin_ia32_andps(tmp, (v4sf)mask);    // sign mask
    tmp = __builtin_ia32_orps(tmp, off*len);        // offset;
    return (v4sf)centre - tmp;
#else
    const int oct = Octant(centre, ppos);
    const float s = centre.w()*0.25f;
    return float4(
        centre.x() + s * ((oct&1) ? (real)1.0 : (real)-1.0),
        centre.y() + s * ((oct&2) ? (real)1.0 : (real)-1.0),
        centre.z() + s * ((oct&4) ? (real)1.0 : (real)-1.0),
        centre.w()*0.5f
        );
#endif
  }



  /********/

  void deleteLeaf(const int leafIdx)
  {
    leafPool.push(leafIdx);
  }
  template<const bool WITHBODY>
  Cell newLeaf(const Body &body = Body())
  {
    if (leafPool.empty())
    {
      leafList.push_back(WITHBODY ? Leaf(body) : Leaf());
      return Cell(reverse_int(leafList.size() - 1), ncell++);
    }
    else
    {
      const int addr = leafPool.top();
      leafPool.pop();
      leafList[addr] = WITHBODY ? Leaf(body) : Leaf();
      return Cell(reverse_int(addr), ncell++);
    }
  }

  void insert(const Body &new_body)
  {
    int  child_idx = 0;                      /* child idx inside a node */
    Cell child     = Cell(Cell::ROOT);       /* child */
    int locked     = 0;                      /* cell that needs to be updated */
    int depth      = 0;                      /* depth */ 
#ifndef NDEBUG
    const int n_nodes_max = cellList.size();
#endif

    float4 centre =  root_centre;

    /* walk the tree to find the first leaf or empty cell */
    while (child.isNode())  /* if non-negative, means it is a tree-node */
    {
      const Cell node = child;
      depth++;

      child_idx  = Octant(centre, new_body.pos());
      centre     = child_centre(centre, new_body.pos());

      locked = node.addr() + child_idx;
      assert(locked < n_nodes_max);
      child  = cellList[locked];
    }

    /* locked on the cell that needs to be updated */

    if (child.isEmpty())  /* the cell is empty, make it a leaf */
    {
      assert(child.isClean());
      cellList[locked] = newLeaf<true>(new_body);
    }
    else          /* this is already a leaf */
    {
      assert(child.isLeaf());
      int leaf_idx = child.leafIdx();
      assert(leaf_idx < (int)leafList.size());
      while(leafList[leaf_idx].isFull())   /* if leaf is full split it */
      {
        depth++;
        nnode++;
        const int cfirst = nnode<<3;
        assert(cfirst+7 < n_nodes_max);
        cellList[locked].set_addr(cfirst);

        assert(leaf_idx >= 0);
        assert(leaf_idx < get_nleaf());
        const Leaf leaf = leafList[leaf_idx];
        deleteLeaf(leaf_idx);

        for (int i = 0; i < leaf.nb(); i++)
        {
          const Body &body = leaf[i];
          child_idx  = Octant(centre, body.pos());
          Cell &cell = cellList[cfirst + child_idx];
          if (cell.isEmpty())
          {
            cell = newLeaf<true>(body);
          }
          else
          {
            assert(cell.isLeaf());
            assert(!leafList[cell.leafIdx()].isFull());
            leafList[cell.leafIdx()].insert(body);
          }
        }

        child_idx = Octant(centre, new_body.pos());
        centre    = child_centre(centre, new_body.pos());

        locked = cfirst + child_idx;
        Cell &child  = cellList[locked];
        if (child.isEmpty())
          child = newLeaf<false>();
        assert(child.isLeaf());
        leaf_idx = child.leafIdx();
      }
      leafList[leaf_idx].insert(new_body);  /* push particle into a leaf */
    }

    this->depth = __max(this->depth, depth);
  }

  /**************/

  bool remove(const Body &body)
  {
    float4 centre =  root_centre;
    
#ifndef NDEBUG
    const int n_nodes_max = cellList.size();
#endif
    std::stack<int> path;

    Cell child     = Cell(Cell::ROOT);       /* child */
    int locked     = 0;                      /* cell that needs to be updated */
    int child_idx  = 0;
    while (child.isNode()) 
    {
      const Cell node = child;
      depth++;

      child_idx  = Octant(centre, body.pos());
      centre     = child_centre(centre, body.pos());

      locked = node.addr() + child_idx;
      assert(locked < n_nodes_max);
      child  = cellList[locked];
      path.push(locked);
      cellList[locked].setTouched();
    }

    assert(!child.isEmpty());
    if (child.isEmpty()) return false;

    assert(child.isLeaf());
    const int leaf_idx = child.leafIdx();
    assert(leaf_idx < (int)leafList.size());
    Leaf &leaf = leafList[leaf_idx];
    int idx = 0;
    for (idx = 0; idx < leaf.nb(); idx++)
      if (body.idx() == leaf[idx].idx())
        break;

    assert(idx < leaf.nb());
    if (idx == leaf.nb()) return false;
    leaf.remove(idx);

    if (!leaf.isEmpty()) 
    {
      while(!path.empty())
      {
        path.pop();
      }
    }

    return true;
  } 

  /**************/

  template<const bool ROOT>  /* must be ROOT = true on the root node (first call) */
    void tree_dump(std::vector<int> &list, const int node = 0) const
    {
      if (ROOT)
      {
        for (int k = 0; k < 8; k++)
          tree_dump<false>(list, k);
      }
      else
      {
        assert(node < (int)cellList.size());
        const Cell cell = cellList[node];
        if (cell.isEmpty()) return;
        if (cell.isNode())
        {
          for (int k = 0; k < 8; k++)
            tree_dump<false>(list, cell.addr()+k);
        }
        else
        {
          assert(cell.isLeaf());
          const int leaf_idx = cell.leafIdx();
          assert(leaf_idx < get_nleaf());
          const Leaf &leaf = leafList[leaf_idx];
          for (int i = 0; i < leaf.nb(); i++)
            list.push_back(leaf[i].idx());
        }
      }
    }

  /**************/

  boundary rootBoundary()
  {
    boundary bnd;
    for (int k = 0; k < 8; k++)
      if (!cellList[k].isEmpty())
      {
        assert(!cellList[k].isTouched());
        bnd.merge(innerBnd[cellList[k].id()]);
      }
    return bnd;
  }

  template<const bool ROOT>  /* must be ROOT = true on the root node (first call) */
    boundary inner_boundary(const int addr = 0)
    {
      boundary bnd;
      if (ROOT)
      {
        innerBnd.resize(ncell);
        for (int k = 0; k < 8; k++)
          if (!cellList[k].isEmpty() && cellList[k].isTouched())
            bnd.merge(innerBnd[cellList[k].id()] = inner_boundary<false>(k));
        treeReady = true;
        return bnd;
      }
      else
      {
        assert(addr < (int)cellList.size());
        Cell &cell = cellList[addr];
        assert (!cell.isEmpty());
        assert(cell.id() >= 0);

        innerBnd[cell.id()] = bnd;

        if (cell.isNode())
        {
          for (int k = cell.addr(); k < cell.addr()+8; k++)
            if (!cellList[k].isEmpty() && cellList[k].isTouched())
              innerBnd[cell.id()].merge(inner_boundary<false>(k));
        }
        else
        {
          const Leaf &leaf = leafList[cell.leafIdx()];
          for (int i= 0; i < leaf.nb(); i++)
          {
            const vec3 jpos = leaf[i].pos();
            innerBnd[cell.id()].merge(jpos);
          }
        }
        cell.unsetTouched();
        return innerBnd[cell.id()];
      }
    }

  /**************/

  template<const bool ROOT>  /* must be ROOT = true on the root node (first call) */
    int sanity_check(const int addr = 0, const boundary &parent_bnd = boundary(), int nb = 0) const
    {
      if (ROOT)
      {
        assert(isTreeReady());
        for (int k = 0; k < 8; k++)
          if (!cellList[k].isEmpty())
            nb = sanity_check<false>(k, innerBnd[cellList[k].id()], nb);
      }
      else
      {
        assert(addr < (int)cellList.size());
        const Cell &cell = cellList[addr];
        if (cell.isEmpty()) return nb;
        assert(cell.id() >= 0);

        if (cell.isNode())
        {
          for (int k = cell.addr(); k < cell.addr()+8; k++)
            nb = sanity_check<false>(k, innerBnd[cell.id()], nb);
        }
        else
        {
          assert(cell.isLeaf());
          assert(cell.leafIdx() < get_nleaf());
          const Leaf &leaf = leafList[cell.leafIdx()];
          for (int i = 0; i < leaf.nb(); i++)
          {
            const vec3 jpos = leaf[i].pos();
            assert(overlapped(innerBnd[cell.id()], jpos));
            assert(overlapped(parent_bnd, jpos));
            nb++;
          }
        }
      }
      return nb;
    }


  /**************/

  template<const bool ROOT>  /* must be ROOT = true on the root node (first call) */
    int range_search(const Body &body,
        const int addr = 0, const boundary &ibnd = boundary(), int nb = 0) const
    {
      if (ROOT)
      {
        assert(isTreeReady());
        const boundary ibnd(body.pos());
        for (int k = 0; k < 8; k++)
          if (!cellList[k].isEmpty())
            if (!not_overlapped(ibnd, innerBnd[cellList[k].id()]))
              nb = range_search<false>(body, k, ibnd, nb);
      }
      else
      {
        const Cell cell = cellList[addr];
        if (cell.isNode())
        {
          for (int k = 0; k < 8; k++)
            if (!cellList[cell.addr()+k].isEmpty())
              if (!not_overlapped(ibnd, innerBnd[cellList[cell.addr()+k].id()]))
                nb = range_search<false>(body, cell.addr()+k, ibnd, nb);
        }
        else
        {
          const Leaf &leaf  = leafList[cell.leafIdx()];
          const float4 ipos = body.pos();
          const real     h2 = ipos.w()*ipos.w();
          for (int i = 0; i < leaf.nb(); i++)
          {
            const float4 jpos = leaf[i].pos();
            if ((ipos - jpos).norm2() < h2)
              nb++;
          }
        }
      }

      return nb;
    }
};


#endif /* __OCTREE_H__ */
