#ifndef __OCTREE_H__
#define __OCTREE_H__

#if 0
#define MEMALIGN
#endif

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <stack>
#include <list>
#ifdef MEMALIGN
#include "memalign_allocator.h"
#endif

template<class T>
static inline T __min(const T &a, const T &b) {return a < b ? a : b;}
template<class T>
static inline T __max(const T &a, const T &b) {return a > b ? a : b;}


#define SQR(x) ((x)*(x))

#ifdef __mySSEX__
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
    int addr;
    int id;
   
    Cell() : addr(EMPTY), id(EMPTY) {}
    Cell(const CellType type) : addr(EMPTY+1), id(EMPTY) {}
    Cell(const int _addr, const int _id) : addr(_addr), id(_id) {}
    bool operator==(const Cell &v) const {return addr == v.addr && id == v.id;}
    bool isClean() const { return addr == EMPTY && id == EMPTY;}
    bool isEmpty() const { return addr == EMPTY;}
    bool isLeaf () const { return addr < EMPTY;}
    bool isNode () const { return addr > EMPTY;}

    const int leafIdx() const 
    { 
      assert(addr < EMPTY); 
      return reverse_int(addr); 
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
      Leaf(const Body &body) : _nb(0) {push(body);}
      void push(const Body &body) 
      {
        assert(_nb < NLEAF);
        _list[_nb++] = body;
      }
      bool isFull() const {return _nb == NLEAF;}
      const Body& operator[](const int i) const { return _list[i];}
      int nb() const {return _nb;}
  };



  private:
  float4 root_centre;  /* root's geometric centre & size */
  int depth;         /* tree-depth */
  int nnode;         /* number of tree nodes */
  int ncell;         /* number of tree cells: node + leaf */
  bool treeReady;    /* this is a flag if tree is ready for walks */
  Cell    ::Vector cellList;
  Leaf    ::Vector leafList;
  boundary::Vector cellBnd;

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
    cellBnd.clear();
  }

  bool isTreeReady() const {return treeReady;}

  static inline int reverse_int(const int a) {return BODYX-a;}

  inline int Octant(const float4 lhs, const float4 rhs) 
  {
#ifdef __mySSEX__
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
#ifdef __mySSEX__
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

  Cell newLeaf()
  {
    leafList.push_back(Leaf());
    return Cell(reverse_int(leafList.size() - 1), ncell++);
  }
  Cell newLeaf(const Body &body)
  {
    leafList.push_back(Leaf(body));
    return Cell(reverse_int(leafList.size() - 1), ncell++);
  }

  void push(const Body &new_body)
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

      locked = node.addr + child_idx;
      assert(locked < n_nodes_max);
      child  = cellList[locked];
    }

    /* locked on the cell that needs to be updated */

    if (child.isEmpty())  /* the cell is empty, make it a leaf */
    {
      assert(child.isClean());
      cellList[locked] = newLeaf(new_body);
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
        cellList[locked].addr = cfirst;

        assert(leaf_idx >= 0);
        assert(leaf_idx < get_nleaf());
        const Leaf leaf = leafList[leaf_idx];

        for (int i = 0; i < leaf.nb(); i++)
        {
          const Body &body = leaf[i];
          child_idx = Octant(centre, body.pos());
          Cell &cell = cellList[cfirst + child_idx];
          if (cell.isEmpty())
          {
            cell = newLeaf(body);
          }
          else
          {
            assert(cell.isLeaf());
            assert(!leafList[cell.leafIdx()].isFull());
            leafList[cell.leafIdx()].push(body);
          }
        }

        child_idx = Octant(centre, new_body.pos());
        centre    = child_centre(centre, new_body.pos());

        locked = cfirst + child_idx;
        Cell &child  = cellList[locked];
        if (child.isEmpty())
          child = newLeaf();
        assert(child.isLeaf());
        leaf_idx = child.leafIdx();
      }
      leafList[leaf_idx].push(new_body);  /* push particle into a leaf */
    }

    this->depth = __max(this->depth, depth);
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
            tree_dump<false>(list, cell.addr+k);
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

  template<const bool ROOT>  /* must be ROOT = true on the root node (first call) */
    boundary inner_boundary(const int addr = 0)
    {
      boundary bnd;
      if (ROOT)
      {
        cellBnd.resize(ncell);
        for (int k = 0; k < 8; k++)
          if (!cellList[k].isEmpty())
            bnd.merge(cellBnd[cellList[k].id] = inner_boundary<false>(k));
        treeReady = true;
        return bnd;
      }
      else
      {
        assert(addr < (int)cellList.size());
        const Cell &cell = cellList[addr];
        assert (!cell.isEmpty());
        assert(cell.id >= 0);

        cellBnd[cell.id] = bnd;

        if (cell.isNode())
        {
          for (int k = cell.addr; k < cell.addr+8; k++)
            if (!cellList[k].isEmpty())
              cellBnd[cell.id].merge(inner_boundary<false>(k));
        }
        else
        {
          const Leaf &leaf = leafList[cell.leafIdx()];
          for (int i= 0; i < leaf.nb(); i++)
          {
            const vec3 jpos = leaf[i].pos();
            cellBnd[cell.id].merge(jpos);
          }
        }
        return cellBnd[cell.id];
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
            nb = sanity_check<false>(k, cellBnd[cellList[k].id], nb);
      }
      else
      {
        assert(addr < (int)cellList.size());
        const Cell &cell = cellList[addr];
        if (cell.isEmpty()) return nb;
        assert(cell.id >= 0);

        if (cell.isNode())
        {
          for (int k = cell.addr; k < cell.addr+8; k++)
            nb = sanity_check<false>(k, cellBnd[cell.id], nb);
        }
        else
        {
          assert(cell.isLeaf());
          assert(cell.leafIdx() < get_nleaf());
          const Leaf &leaf = leafList[cell.leafIdx()];
          for (int i = 0; i < leaf.nb(); i++)
          {
            const vec3 jpos = leaf[i].pos();
            assert(overlapped(cellBnd[cell.id], jpos));
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
            if (!not_overlapped(ibnd, cellBnd[cellList[k].id]))
              nb = range_search<false>(body, k, ibnd, nb);
      }
      else
      {
        const Cell cell = cellList[addr];
        if (cell.isNode())
        {
          for (int k = 0; k < 8; k++)
            if (!cellList[cell.addr+k].isEmpty())
              if (!not_overlapped(ibnd, cellBnd[cellList[cell.addr+k].id]))
                nb = range_search<false>(body, cell.addr+k, ibnd, nb);
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
