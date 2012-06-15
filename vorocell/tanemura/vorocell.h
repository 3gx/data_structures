#ifndef __VOROCELL_H__
#define __VOROCELL_H__

/* Tanemura's algorithm */

#include <cassert>
#include <cvector>
#include <stack>
#include <map>

template<class T, const int N>
class Array
{
  private:
    int n;
    T data[N];

  public:
    Array(const int _n = 0) : n(_n) {}
    Array(const std::vector<T> &_data)
    {
      n = _data.size();
      assert(n <= N);
      for (int i = 0; i < n; i++)
        data[i] = _data[i];
    }
    void clear()
    {
      n = 0;
    }
    const T& operator[](const int i) const {return data[i];}
          T& operator[](const int i)       {return data[i];}

    void push_back(const T &t) {assert(n<N); data[n++] = t;}
    int size() const { return n; }
    void resize(const int size) 
    {
      assert (size <= N);
      n = size;
    }

    int capacity() const {return N;}

    bool erase(const T &t) 
    {
      for (int i = 0; i < n; i++)
        if (data[i] == t)
        {
          n--;
          std::swap(data[i], data[n]);
          return true;
        }
      return false;
    }

    T erase(const int i)
    {
      assert(i >= 0);
      assert(i < n);
      n--;
      std::swap(data[i], data[n]);
      return data[n];
    }
};

template<const int N, const int NCRIT>
struct ConnectivityMatrix2D
{
  private:
    int matrix[N*(N+1)/2];
    int _Flags[N];

  public:

    ConnectivityMatrix2D() 
    {
      clear();
    }

    clear()
    {
      for (int i = 0; i < N*(N+1)/2; i++)
        matrix[i] = 0;

      for (int i = 0; i < N; i++)
        connected[i] = 0;
    }
    
    /* i: x
     * j: y */
    unsigned int  operator()(const int i, const int j) const { return matrix[map(i,j)]; }
    unsigned int  inc       (const int i, const int j) 
    {
      unsigned int &val = matrix[map(i,j)];
      val++;
      const int inc = (int)(val == NCRIT) - (int)(val == NCRIT+1);
      _flags[i] += inc;
      _flags[j] += inc;
    }

    /************/

    bool isFullyConnected(const int i) const {return _flags[i] == 0;}

  private:

    /* i: x
     * j: y */
    int map(const int i, const int j) const
    {
#if 1 /* LOWER_PACKED */
        return i + ((2*N-j)*(j-1)>>1);
#else /* UPPER_PACKED */
        return i + (j*(j-1)>>1);
    }
};

template<const int N>
struct Vorocell
{
  typedef ConnectivityMatrix2D<N, 1>      Edges;
  typedef ConnectivityMatrix2D<N, 2>  Triangles;
  private:
    Edges     edges;
    Triangles triangles;

    std::deque<    int     > vertexQueue;
    std::vector<Tetrahedron>  tetrahedra;
    std::vector<   int     >      nbList;

    std::vector<bool> vertexQueued;
    std::vector<int >        edges;
    std::vector<int >    triangles;

    Array<N>    faceVtx[N];
    int vertexCompleted[N];

  public:
    Vorocell(const int Particle::Vector &ptcl)
    {
      clear();
      assert((int)ptcl.size() <= N);
    }
    
    clear()
    {
      edges      .clear();
      triangles  .clear();
      vertexQueue.clear();
      tetrahedra .clear();
      for (int i = 0; i < N; i++)
      {
        vertexCompleted[i] = 0;
        faceVtx[i].clear();
      }
    }


  private:

    void buildCell()
    {
    }

    void completeCell()
    {
      while (!vertexQueue.empty())
      {
        /* step 4.2:
         *  extract vertex from the list
         */
        const int iVertex = vertexQueue.front();
        vertexQueue.pop_front();

        nbList.push_back(iVertex);
        
        if (edges.isFullyConnected(iVertex))
          vertexCompleted[iVertex] = 1;

        /* step 4.3:
         *  the vertex is completed, proceed to the next one
         */
        if (vertexCompleted[iVertex]) continue;


        /* step 4.4:
         *  find a tetrahedron (i, iVertex, jVertex, kVertex) with at least one
         *  incomplete face
         */

        /*    a): get a single counted edge from the iVertex , and obtain 2nd vertex */
        Edge::Iterator &edgeIt = edgeList.find(Edge(iVertex, Edge::ANY_VERTEX));
        assert(edgeIt != edgeList.end());
        const Edge & edge = edgeIt->first;
        int jVertex = edge.jVertex();
        edgeIt->second++;
        assert(2 == edgeIt->second);
        edgeList.erase(edge);  /* since edge becomes double connected, remove it */

        /*    b) find tetrahedron that shares this edge, and obtain the 3d vertex */
        Tetrahedron::Iterator &tetraIt = tatraList.find(Tetrahedron(iVertex, jVertex, Tetrahedron::ANY_VERTEX));
        assert(tetraIt != tetraList.end());
        const Tetrahedron &tetrahedron = tetraIt->first;
        int kVertex = tetrahedron.kVertex();

        /* step 4.5-4.7 */

        while(triangles[map(iVertex, jVertex)] < 2)
        {

          /* step 4.5 - 4.6: 
           *  search a vertex on the opposite side of the kVertex 
           *  (in the half-space bounded by ijFace that does not contain kVertex)
           */
          const vec3 &ipos = siteList[iVertex];
          const vec3 &jpos = siteList[jVertex];
          const vec3 &kpos = siteList[kVertex];
          const Plane plane(ipos, jpos);
          const int  sideK = plane(kpos) > 0.0;
          real largeNum = +1e10;
          vec3  cpos(0.0);
          int lVertex = -1;
          for (int i = 0; i < nVertex; i++)
          {
            const vec3 &pos = siteList[i];
            const int  side = plane(pos) > 0.0;
            const real dist = pos*(pos + cpos) + largeNum;
            const bool chck = !vertexCompleted[i] &&
              (i != iVertex) && (i != jVertex) && (i != kVertex);
            
            if (c < 0.0 && side^sideK = 1 && chck)
            {
              real radius = 0.0;
              cpos = sphere(ipos, jpos, pos, radius);
              largeNum = 0.0;
              lVertex = i;
            }
          }
          assert(lVertex >= 0);

          /* step 4.7:
           *  register new tetrahedron (iVertex, jVertex, lVertex) 
           */
          tetrahedronList.push_back(Tetrahedron(iVertex, jVertex, lVertex));
          if (!queuedVertex[lVertex])
          {
            vertexQueue.push_back(lVertex);
            queuedVertex[lVertex] = 1;
          }
          faceVtx[iVertex].push_back(tetrahedronList.size()-1);

          edges.inc(iVertex, jVertex);
          edges.inc(iVertex, lVertex);
          edges.inc(jVertex, lVertex);
 
          triangles.inc(iVertex, jVertex);
          triangles.inc(iVertex, lVertex);
          triangles.inc(jVertex, lVertex);
          kVertex = jVertex;
          jVertex = lVertex;
        }

        /* step 4.8: */
        vertexCompleted[iVertex] = 1;
      }
    }



};

#endif /* __VOROCELL_H__ */
