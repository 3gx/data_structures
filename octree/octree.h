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


#include "vector3.h"

typedef float real;
typedef vector3<real> vec3;
#if 0
#include "boundary.h"
typedef Boundary<real> boundary;
#else
#include "boundary4.h"
#endif


struct Particle
{
  typedef std::vector<Particle> Vector;
  vec3 pos;
  int  id;
  real h;
  int nb;
  int iPad1, iPad2;

  Particle() {}
  Particle(const vec3 &_pos, const int _id, const real _h = 0.0) :
    pos(_pos), id(_id), h(_h), nb(0) {}
};


struct Octree
{
  enum {NLEAF =  64};
  enum {EMPTY =  -1};
  enum {BODYX =  -2};

  struct Body
  {
#ifdef MEMALIGN
    typedef std::vector<Body, __gnu_cxx::malloc_allocator<Body, 64> > Vector;
#else
    typedef std::vector<Body> Vector;
#endif
    private:
    float4 Packed_pos;
    int _idx;
    int iPad1, iPad2, iPad3; /* padding */

    public:
    Body() {}
    Body(const vec3 &pos, const int idx, const float h = 0.0f) : Packed_pos(float4(pos.x, pos.y, pos.z, h)), _idx(idx) {}
    Body(const Particle &p, const int idx) : Packed_pos(float4(p.pos.x, p.pos.y, p.pos.z, p.h)), _idx(idx) {}
    Body(const vec3 &pos, const float h, const int idx = -1) : Packed_pos(pos.x, pos.y, pos.z, h), _idx(idx) {};
    int idx() const {return _idx;}
    float4 packed_pos() const {return Packed_pos;}
    float          h()  const {return Packed_pos.w();}
    vec3   vector_pos() const {return vec3(Packed_pos.x(), Packed_pos.y(), Packed_pos.z());}
#ifdef __mySSE__
    Body& operator=(const Body &rhs)
    {
      v4sf *dst = (v4sf*)this;
      v4sf *src = (v4sf*)&rhs;
#if 1
      dst[0] = src[0];
      dst[1] = src[1];
#else
      const v4sf a = src[0];
      const v4sf b = src[1];
      __builtin_ia32_movntps((float*)(dst+0), a);
      __builtin_ia32_movntps((float*)(dst+1), b);
#endif
      return *this;
    }
#endif
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
    int _np;
    int iPad;

    public:
    Cell() : _addr(EMPTY), _id(EMPTY), _np(0) {}
    Cell(const CellType type) : _addr(EMPTY+1), _id(EMPTY), _np(0) {}
    Cell(const int addr, const int id, const int np = 0) : _addr(addr), _id(id), _np(np) {setTouched();}
    bool operator==(const Cell &v) const {return _addr == v._addr && _id == v._id;}
    bool isClean() const { return _addr == EMPTY && _id == EMPTY;}
    bool isEmpty() const { return _addr == EMPTY;}
    bool isLeaf () const { return _addr < EMPTY;}
    bool isNode () const { return _addr > EMPTY;}

    int operator++(const int) {_np++; return _np;}
    int operator--(const int) {_np--; return _np;}
    int np() const {return _np;}

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

  struct Boundaries
  {
    typedef std::vector<Boundaries> Vector;
    private:
    boundary _inner;
    boundary _outer;

    public:
    Boundaries(const boundary &inner = boundary(), const boundary &outer = boundary()) :
      _inner(inner), _outer(outer) {}
    Boundaries& merge(const Boundaries &rhs) 
    {
      _inner.merge(rhs._inner);
      _outer.merge(rhs._outer);
      return *this;
    }

    const boundary& inner() const { return _inner; }
    const boundary& outer() const { return _outer; }

  };

  template<const int N>
  struct GroupT
  {
    typedef std::vector<GroupT> Vector;
    private:
      int _nb;           /* number of bodies in the list */
      Body _list[N];  /* idx of the body */

    public:
      GroupT() : _nb(0) {}
      GroupT(const Body &body) : _nb(0) {insert(body);}
      void insert(const Body &body) 
      {
        assert(_nb < N);
        _list[_nb++] = body;
      }
      void remove(const int idx)
      {
        assert(idx < _nb);
        _nb--;
        for (int i = idx; i < _nb; i++)
          _list[i] = _list[i+1];
      }
      bool  isFull() const {return N == _nb;}
      bool isEmpty() const {return 0 == _nb;}
      const Body& operator[](const int i) const { return _list[i];}
      int nb() const {return _nb;}

      boundary innerBoundary() const
      {
        boundary bnd;
        for (int i = 0; i < _nb; i++)
          bnd.merge(_list[i].vector_pos());
        return bnd;
      }
      boundary outerBoundary() const
      {
        boundary bnd;
        for (int i = 0; i < _nb; i++)
          bnd.merge(_list[i].packed_pos());
        return bnd;
      }
      Boundaries computeBoundaries() const
      {
        return Boundaries(innerBoundary(), outerBoundary());
      }
  };
  typedef GroupT<NLEAF>  Leaf;

  private:
  float4 root_centre;  /* root's geometric centre & size */
  int depth;         /* tree-depth */
  int nnode;         /* number of tree nodes */
  int ncell;         /* number of tree cells: node + leaf */
  bool treeReady;    /* this is a flag if tree is ready for walks */
  Cell    ::Vector   cellList;    /* list of tree cells (node+leaf) */
  Leaf    ::Vector   leafList;    /* list of leaves */
  std::stack<int>    leafPool;    /* pool of available leaves */
  std::stack<int>    cellPool;    /* pool of available cells */
  Boundaries::Vector bndsList;    /* list of cell  boundaries */

  std::vector<int> leafList_addr;

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
    bndsList.clear();
    leafList_addr.clear();
    leafPool = std::stack<int>();
  }

  bool isTreeReady() const {return treeReady;}
  
  static inline int reverse_int(const int a) {return BODYX-a;}

  inline int Octant(const float4 lhs, const float4 rhs) 
  {
#ifdef __mySSE__
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
#ifdef __mySSE__
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

  void deleteCell(const int cellId)
  {
    cellPool.push(cellId);
  }
  void deleteLeaf(const int leafIdx)
  {
    leafPool.push(leafIdx);
  }
  template<const bool WITHBODY>
  Cell newLeaf(const Body &body = Body())
  {
    const int cell = cellPool.empty() ? ncell++ : cellPool.top();
    if (!cellPool.empty()) cellPool.pop();

    if (leafPool.empty())
    {
      const int addr = leafList.size();
      leafList.push_back(WITHBODY ? Leaf(body) : Leaf());
      return Cell(reverse_int(addr), cell, WITHBODY ? 1 : 0);
    }
    else
    {
      const int addr = leafPool.top();
      leafPool.pop();
      leafList[addr] = WITHBODY ? Leaf(body) : Leaf();
      return Cell(reverse_int(addr), cell, WITHBODY ? 1 : 0);
    }
  }

  int get_np() const 
  {
    int np = 0;
    for (int k = 0; k < 8; k++)
      np += cellList[k].np();
    return np;
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

      child_idx  =       Octant(centre, new_body.packed_pos());
      centre     = child_centre(centre, new_body.packed_pos());

      locked = node.addr() + child_idx;
      assert(locked < n_nodes_max);
      child  = cellList[locked];
      cellList[locked]++;
    }

    /* locked on the cell that needs to be updated */

    if (child.isEmpty())  /* the cell is empty, make it a leaf */
    {
      assert(child.isClean());
      cellList[locked] = newLeaf<true>(new_body);
      assert(cellList[locked].np() == leafList[cellList[locked].leafIdx()].nb());
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
          child_idx  = Octant(centre, body.packed_pos());
          Cell &cell = cellList[cfirst + child_idx];
          cell++;
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

        child_idx =       Octant(centre, new_body.packed_pos());
        centre    = child_centre(centre, new_body.packed_pos());

        locked = cfirst + child_idx;
        Cell &child  = cellList[locked];
        if (child.isEmpty())
          child = newLeaf<false>();
        child++;
        assert(child.isLeaf());
        leaf_idx = child.leafIdx();
      }
      leafList[leaf_idx].insert(new_body);  /* push particle into a leaf */
      assert(leafList[leaf_idx].nb() == cellList[locked].np());
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

      child_idx  =       Octant(centre, body.packed_pos());
      centre     = child_centre(centre, body.packed_pos());

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

    if (leaf.isEmpty()) 
    {
      while(!path.empty())
      {
        path.pop();
      }
      deleteLeaf(leaf_idx);
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

  boundary root_innerBoundary()
  {
    boundary bnd;
    for (int k = 0; k < 8; k++)
      if (!cellList[k].isEmpty())
      {
        assert(!cellList[k].isTouched());
        bnd.merge(bndsList[cellList[k].id()].inner());
      }
    return bnd;
  }
  boundary root_outerBoundary()
  {
    boundary bnd;
    for (int k = 0; k < 8; k++)
      if (!cellList[k].isEmpty())
      {
        assert(!cellList[k].isTouched());
        bnd.merge(bndsList[cellList[k].id()].outer());
      }
    return bnd;
  }

  template<const bool ROOT>  /* must be ROOT = true on the root node (first call) */
    Boundaries computeBoundaries(const int addr = 0)
    {
      Boundaries bnds;
      if (ROOT)
      {
        bndsList.resize(ncell);
        for (int k = 0; k < 8; k++)
          if (!cellList[k].isEmpty() && cellList[k].isTouched())
            bnds.merge(bndsList[cellList[k].id()] = computeBoundaries<false>(k));
        treeReady = true;
      }
      else
      {
        assert(addr < (int)cellList.size());
        Cell &cell = cellList[addr];
        assert (!cell.isEmpty());
        assert(cell.id() >= 0);


        if (cell.isNode())
        {
          for (int k = cell.addr(); k < cell.addr()+8; k++)
            if (!cellList[k].isEmpty() && cellList[k].isTouched())
              bnds.merge(computeBoundaries<false>(k));
        }
        else
        {
          const Leaf &leaf = leafList[cell.leafIdx()];
          bnds.merge(leaf.computeBoundaries());
        }
        bndsList[cell.id()] = bnds;
        cell.unsetTouched();
      }
      return bnds;
    }

  /**************/

  template<const bool ROOT>  /* must be ROOT = true on the root node (first call) */
    int sanity_check(const int addr = 0) const
    {
      int nb = 0;
      if (ROOT)
      {
        assert(isTreeReady());
        for (int k = 0; k < 8; k++)
          if (!cellList[k].isEmpty())
            nb += sanity_check<false>(k);
      }
      else
      {
        assert(addr < (int)cellList.size());
        const Cell &cell = cellList[addr];
        if (cell.isEmpty()) return nb;
        assert(cell.id() >= 0);

        if (cell.isNode())
        {
          assert(cell.np() > NLEAF);
          for (int k = cell.addr(); k < cell.addr()+8; k++)
            nb += sanity_check<false>(k);
        }
        else
        {
          assert(cell.np() <= NLEAF);
          assert(cell.isLeaf());
          assert(cell.leafIdx() < get_nleaf());
          const Leaf &leaf = leafList[cell.leafIdx()];
          for (int i = 0; i < leaf.nb(); i++)
          {
            const vec3 jpos = leaf[i].vector_pos();
            assert(overlapped(bndsList[cell.id()].inner(), jpos));
            assert(overlapped(bndsList[cell.id()].outer(), leaf[i].packed_pos()));
            nb++;
          }
        }
      }
      return nb;
    }

  /**************/

  int nLeaf() const {return leafList_addr.size();}
  const Leaf& getLeaf(const int i) const {return leafList[i];}
  template<const bool ROOT>  /* must be ROOT = true on the root node (first call) */
    void buildLeafList(const int addr = 0) 
    {
      if (ROOT)
      {
        leafList_addr.clear();
        leafList_addr.reserve(128);
        assert(isTreeReady());
        for (int k = 0; k < 8; k++)
          if (!cellList[k].isEmpty())
            buildLeafList<false>(k);
      }
      else
      {
        assert(addr < (int)cellList.size());
        const Cell &cell = cellList[addr];
        if (cell.isEmpty()) return;
        assert(cell.id() >= 0);

        if (cell.isNode())
        {
          for (int k = cell.addr(); k < cell.addr()+8; k++)
            buildLeafList<false>(k);
        }
        else
        {
          assert(cell.isLeaf());
          assert(cell.leafIdx() < get_nleaf());
          leafList_addr.push_back(cell.leafIdx());
        }
      }
    }

  /**************/
  
  template<const bool ROOT, const int N>  /* must be ROOT = true on the root node (first call) */
    void buildGroupList(
        std::vector< GroupT<N> > &groupList,
        const int addr = 0) const
    {
      assert(NLEAF <= N);
      if (ROOT)
      {
        assert(isTreeReady());
        for (int k = 0; k < 8; k++)
          if (!cellList[k].isEmpty())
            buildGroupList<false>(groupList, k);
      }
      else
      {
        assert(addr < (int)cellList.size());
        const Cell &cell = cellList[addr];
        if (cell.isEmpty()) return;
        assert(cell.id() >= 0);

        if (cell.np() <= N)
        {
          groupList.push_back(GroupT<N>());
          extractBodies(groupList.back(), addr);
        }
        else
        {
          assert(!cell.isLeaf());
          for (int k = cell.addr(); k < cell.addr()+8; k++)
            buildGroupList<false>(groupList, k);
        }
      }
    }

  template<const int N>
    void extractBodies(GroupT<N> &group, const int addr) const
    {
      assert(addr < (int)cellList.size());
      const Cell &cell = cellList[addr];
      if (cell.isEmpty()) return;

      if (cell.isNode())
      {
        for (int k = 0; k < 8; k++)
          extractBodies(group, cell.addr() + k);
      }
      else
      {
        assert(cell.isLeaf());
        const Leaf &leaf = leafList[cell.leafIdx()];
        for (int i = 0; i < leaf.nb(); i++)
          group.insert(leaf[i]);
      }
    }

  /**************/

  template<const bool ROOT>  /* must be ROOT = true on the root node (first call) */
    int range_search(const Body &ibody, const boundary &ibnd = boundary(), const int addr = 0, int nb = 0) const
    {
      if (ROOT)
      {
        assert(isTreeReady());
        const boundary ibnd(ibody.packed_pos());
        for (int k = 0; k < 8; k++)
          if (!cellList[k].isEmpty())
            if (!not_overlapped(ibnd, bndsList[cellList[k].id()].inner()))
              nb = range_search<false>(ibody, ibnd, k, nb);
      }
      else
      {
        const Cell cell = cellList[addr];
        if (cell.isNode())
        {
          for (int k = 0; k < 8; k++)
            if (!cellList[cell.addr()+k].isEmpty())
              if (!not_overlapped(ibnd, bndsList[cellList[cell.addr()+k].id()].inner()))
                nb = range_search<false>(ibody, ibnd, cell.addr()+k, nb);
        }
        else
        {
          const Leaf &leaf  = leafList[cell.leafIdx()];
          const float4 ipos = ibody.packed_pos();
          const real     h2 = ibody.h()*ibody.h();
          for (int i = 0; i < leaf.nb(); i++)
          {
            const float4 jpos = leaf[i].packed_pos();
            const float  r2   = (ipos - jpos).norm2();
            if (r2 < h2)
              nb++;
          }
        }
      }

      return nb;
    }

  template<const bool ROOT, const int N>  /* must be ROOT = true on the root node (first call) */
    void range_search(
        int nb[N],
        const GroupT<N> &igroup, const boundary &ibnd = boundary(), const int addr = 0) const
    {
      if (ROOT)
      {
        assert(isTreeReady());
        const boundary ibnd(igroup.outerBoundary());
        for (int i = 0; i < N; i++)
          nb[i] = 0;
        for (int k = 0; k < 8; k++)
          if (!cellList[k].isEmpty())
            if (!not_overlapped(ibnd, bndsList[cellList[k].id()].inner()))
              range_search<false>(nb, igroup, ibnd, k);
      }
      else
      {
        const Cell cell = cellList[addr];
        if (cell.isNode())
        {
          for (int k = 0; k < 8; k++)
            if (!cellList[cell.addr()+k].isEmpty())
              if (!not_overlapped(ibnd, bndsList[cellList[cell.addr()+k].id()].inner()))
                range_search<false>(nb, igroup, ibnd, cell.addr()+k);
        }
        else
        {
#ifdef __mySSE__
          const Leaf &leaf        = leafList[cell.leafIdx()];
          const boundary &leafBnd = bndsList[cell.     id()].inner();
          const v4sf  jmin = leafBnd.min;
          const v4sf  jmax = leafBnd.max;
          const int   ni = igroup.nb();
          const int   nj =   leaf.nb();
          const v4sf *ib = (const v4sf*)&igroup[0];
          const v4sf *jb = (const v4sf*)&  leaf[0];
          const int nj2 = nj<<1;
          for (int i = 0; i < ni; i += 4)
          {
            asm("#eg01");
            const int i2 = i<<1;
            const v4sf ip0 = *(ib + i2 + 0);
            const v4sf ip1 = *(ib + i2 + 2);
            const v4sf ip2 = *(ib + i2 + 4);
            const v4sf ip3 = *(ib + i2 + 6);

            /* check if these i-particles overlap with the leaf */

            const v4sf h0  = __builtin_ia32_shufps(ip0, ip0, 0xFF);
            const v4sf h1  = __builtin_ia32_shufps(ip1, ip1, 0xFF);
            const v4sf h2  = __builtin_ia32_shufps(ip2, ip2, 0xFF);
            const v4sf h3  = __builtin_ia32_shufps(ip3, ip3, 0xFF);

            /*   0     1     2     3   */
            /*  0x00 ,0x55, 0xaa, 0xff */

            const v4sf imin    = __builtin_ia32_minps(__builtin_ia32_minps(ip0-h0, ip1-h1), __builtin_ia32_minps(ip2-h2,ip3-h3));
            const v4sf imax    = __builtin_ia32_maxps(__builtin_ia32_maxps(ip0+h0, ip1+h1), __builtin_ia32_maxps(ip2+h2,ip3+h3));

            const bool skip    = __builtin_ia32_movmskps(__builtin_ia32_orps(
                  __builtin_ia32_cmpltps(jmax, imin),
                  __builtin_ia32_cmpltps(imax, jmin))) & 7;
            if (skip && i+4 < ni) continue;

            /* they do overlap, now proceed to the interaction part */

            const v4sf t0 = __builtin_ia32_unpcklps(ip0, ip2);
            const v4sf t1 = __builtin_ia32_unpckhps(ip0, ip2);
            const v4sf t2 = __builtin_ia32_unpcklps(ip1, ip3);
            const v4sf t3 = __builtin_ia32_unpckhps(ip1, ip3);

            const v4sf ipx = __builtin_ia32_unpcklps(t0, t2);
            const v4sf ipy = __builtin_ia32_unpckhps(t0, t2);
            const v4sf ipz = __builtin_ia32_unpcklps(t1, t3);
            const v4sf iph = __builtin_ia32_unpckhps(t1, t3);
            const v4sf iph2 = iph*iph;

            v4sf inb = {0.0f,0.0f,0.0f,0.0f};
            for (int j = 0; j < nj2; j += 2)
            {
              const v4sf jp = *(jb + j);

#if 0 /*makes it slow*/
              const bool skip    = __builtin_ia32_movmskps(__builtin_ia32_orps(
                    __builtin_ia32_cmpltps(jp,   imin),
                    __builtin_ia32_cmpltps(imax, jp))) & 7;
              if (skip) continue;
#endif

              const v4sf jpx = __builtin_ia32_shufps(jp, jp, 0x00);
              const v4sf jpy = __builtin_ia32_shufps(jp, jp, 0x55);
              const v4sf jpz = __builtin_ia32_shufps(jp, jp, 0xAA);


              const v4sf dx = jpx - ipx;
              const v4sf dy = jpy - ipy;
              const v4sf dz = jpz - ipz;
              const v4sf r2 = dx*dx + dy*dy + dz*dz;

              const v4sf mask = __builtin_ia32_cmpltps(r2, iph2);
#if 0
              const int imask = __builtin_ia32_movmskps(mask);
              if (imask == 0) continue;
#endif
              inb += __builtin_ia32_andps((v4sf){1,1,1,1}, mask);

            }
            const float4 fnb(inb);
            nb[i+0] += (int)fnb.x();
            nb[i+1] += (int)fnb.y();
            nb[i+2] += (int)fnb.z();
            nb[i+3] += (int)fnb.w();
            asm("#eg02");
          }
#else
          const Leaf     &leaf    = leafList[cell.leafIdx()];
          const boundary &leafBnd = bndsList[cell.     id()].inner();
          for (int i = 0; i < igroup.nb(); i++)
          {
            const Body &ibody = igroup[i];
            const float4 ipos = ibody.packed_pos();
            const real     h  = ibody.h();
            const real     h2 = h*h;
            if (overlapped(boundary(ipos), leafBnd))
              for (int j = 0; j < leaf.nb(); j++)
              {
                const float4 jpos = leaf[j].packed_pos();
                const float  r2   = (ipos - jpos).norm2();
                if (r2 < h2)
                  nb[i]++;
              }
          }
#endif
        }
      }
    }
};


#endif /* __OCTREE_H__ */
