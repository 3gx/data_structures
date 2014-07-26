#include "qhull.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

    
#define N 3
using QHull = QHull_t<N>;

template<>
template<>
void QHull::findExtremeSimplex<2>(
    const Vertex *pos, 
    Simplex &simplex,
    const int np)
{
#if 0
  vec_t rMax(-HUGE);
  real_t xMin = +HUGE, xMax = -HUGE;
  // foreach
  for (int i = 0; i < np; i++)
  {
    const auto &p = pos[i];
    if (p[0] < xMin)
    {
      xMin       = p[0];
      simplex[0] = p;
    }
    if (p[0] > xMax)
    {
      xMax       = p[0];
      simplex[1] = p;
    }
  }
#else
  vec_t rMin(+HUGE), rMax(-HUGE);
  std::array<Vertex,NDIM> pMin,pMax;

  for (int i = 0; i < np; i++)
  {
    const auto &p = pos[i];
    for (int l = 0; l < NDIM; l++)
    {
      if (p[l] < rMin[l])
      {
        rMin[l] = p[l];
        pMin[l] = p;
      }
      if (p[l] > rMax[l])
      {
        rMax[l] = p[l];
        pMax[l] = p;
      }
    }
  }
  real_t distMax = 0;
  for (int l = 0; l < NDIM; l++)
    for (int ll = 0; ll < NDIM; ll++)
    {
      const real_t dist = norm2(pMax[l].pos - pMin[ll].pos);
      if (dist > distMax)
      {
        distMax = dist;
        simplex[0] = pMin[l];
        simplex[1] = pMax[ll];
      }

    }
#endif
}

QHull::Vertex::vector readData(std::istream &in)
{
  std::string line;
  if (!std::getline(in, line)) 
  {
    std::cerr << "io: dimension line" << std::endl;
    exit(-1);
  }
  std::istringstream iss;
  size_t dimension;
  {
    iss.str(line);
    if (!(iss >> dimension)) {
      std::cerr << "io: dimension" << std::endl;
      exit(-1);
    }
    {
      using char_type = typename std::string::value_type;
      std::cerr << "#command line:";
      std::istreambuf_iterator<char_type> const ibeg(iss), iend;
      std::copy(ibeg, iend, std::ostreambuf_iterator<char_type>(std::cerr));
      std::cerr << '\n';
    }
    iss.clear();
  }
  if (!std::getline(in, line)) {
    std::cerr << "io: count line" << std::endl;
    exit(-1);
  }

  typedef QHull::Vertex vtx_t;
  size_t count;
  assert(dimension == QHull::NDIM);
  {
    iss.str(line);
    if (!(iss >> count)) {
      std::cerr << "io: count" << std::endl;
      exit(-1);
    }
    iss.clear();
  }
  vtx_t::vector pos(count);
  for (size_t i = 0; i < count; ++i) 
  {
    if (!std::getline(in, line)) {
      std::cerr << "io: line count error" << std::endl;
      exit(-1);
    }
    vtx_t &p = pos[i];
    p.idx = i;
    {
      iss.str(line);
      for (size_t j = 0; j < dimension; ++j) {
        if (!(iss >> p[j])) {
          std::cerr << "io: bad value at line " << j << " of data" << std::endl;
          exit(-1);
        }
      }
      iss.clear();
    }
  }
  std::cerr << "#D = " << dimension << '\n';
  std::cerr << "#N = " << count << '\n';

  return pos;
}

template<int NN>
void dumpSimplex(const QHull::Simplex &simplex, std::ostream &out);

template<> void dumpSimplex<1>(const QHull::Simplex &simplex, std::ostream &out)
{
  const int NDIM = 1;
  for (int i = 0; i < NDIM+1; i++)
  {
    const QHull::Vertex &p = simplex[i%(NDIM+1)];
    for (int l = 0; l < NDIM; l++)
      out << p[l] << " ";
    out << "\n";
  }
}
template<> void dumpSimplex<2>(const QHull::Simplex &simplex, std::ostream &out)
{
  const int NDIM = 2;
  for (int i = 0; i < NDIM+2; i++)
  {
    const QHull::Vertex &p = simplex[i%(NDIM+1)];
    for (int l = 0; l < NDIM; l++)
      out << p[l] << " ";
    out << "\n";
  }
}
template<> void dumpSimplex<3>(const QHull::Simplex &simplex, std::ostream &out)
{
  const int vertexList[4][4] =  { {0,1,2,0}, {0,1,3,0}, {1,2,3,1}, {2,0,3,2} };
  for (auto &vtx : vertexList)
    for (auto i : vtx)
    {
      for (int l = 0; l < 3 ; l++)
        out << simplex[i][l] << " ";
      out << "\n";
    }
}


void dumpGnuPlot(const QHull::Vertex::vector &pos, const QHull &q, std::ostream &out)
{
  out << "clear\n";
  out << "set autoscale\n";
  const int NDIM = QHull::NDIM;
  switch (NDIM)
  {
    case 2:
      out << "plot";
      break;
    case 3:
      out << "splot";
      break;
    default:
      assert(0);
  };
  out << " '-' with points notitle, '-' with lines lc 3 notitle, '-' with lines lc 2 notitle";
  out << ";\n";

  {
    const int np = pos.size();
    for (int i = 0; i < np; i++)
    {
      const QHull::Vertex &p = pos[i];
      for (int l = 0; l < NDIM; l++)
        out << p[l] << " ";
      out << "\n";
    }
  }
  out << "e\n";
#if 0
  dumpSimplex<NDIM>(q.extremeSimplex, out);
#else
  {
    const QHull::Simplex &simplex = q.extremeSimplex;
    const int vertexList[3][4] =  { {0,1,3,0}, {1,2,3,1}, {2,0,3,2} };
    for (auto &vtx : vertexList)
      for (auto i : vtx)
      {
        for (int l = 0; l < 3 ; l++)
          out << simplex[i][l] << " ";
        out << "\n";
      }
    out << "e\n";
    const int vtx[4] =  {0,1,2,0};
    for (auto i : vtx)
    {
      for (int l = 0; l < 3 ; l++)
        out << simplex[i][l] << " ";
      out << "\n";
    }
  }
#endif
  out << "e\n";
    
}

#if 0
#define TESTPLANEQ
void testPlaneEquations()
{
  using Vertex = QHull::Vertex;
  using Basis  = QHull::Basis;
  Basis v;
  //
  v[0][0] = 0.0;
  v[0][1] = 0.0;
  v[0][2] = 0.0;
  //
  v[1][0] = 1.0;
  v[1][1] = 0.0;
  v[1][2] = 0.0;
  //
  v[2][0] = 0.5;
  v[2][1] = 1.0;
  v[2][2] = 0.5;

  Vertex p;
  p[0] = 0.5;
  p[1] = 0.5;
  p[2] = 0.5;

  std::cout << "splot";
  std::cout << " '-' with points, '-' with lines lc 3 notitle, '-' with vector lc 2 notitle";
  std::cout << ";\n";
  std::cout << p[0] << " ";
  std::cout << p[1] << " ";
  std::cout << p[2] << " ";
  std::cout << "\ne\n";
  const int vtx[4] =  {0,1,2,0};
  for (auto i : vtx)
  {
    for (int l = 0; l < 3 ; l++)
      std::cout << v[i][l] << " ";
    std::cout << "\n";
  }
  std::cout << "e\n";
  const auto &plane = QHull::planeEquation(v, p, +1.0);
#if 0
  std::cout << " 0.25 0.25 0.0 ";
  std::cout << 0.25+plane.first[0] << " ";
  std::cout << 0.25+plane.first[1] << " ";
  std::cout << plane.first[2] << " \n";
#else
  std::cout << " 0 0 0.0 ";
  std::cout << plane.first[0] << " ";
  std::cout << plane.first[1] << " ";
  std::cout << plane.first[2] << " \n";
#endif
}
#endif


int main(int argc, char*argv[])
{
#ifdef TESTPLANEQ
  testPlaneEquations();
  return 0;
#endif
  std::ifstream ifs;
  if (argc > 1)
  {
    ifs.open(argv[1]);
    if (!ifs.is_open())
    {
      std::cerr << "cannot open the file" << std::endl;
      exit(-1);
    }
  }

  std::ofstream ofs;
  if (argc > 2)
  {
    ofs.open(argv[2]);
    if (!ofs.is_open())
    {
      std::cerr << "cannot open the file" << std::endl;
      exit(-1);
    }
  }
  std::istream &in  = ifs.is_open() ? ifs : std::cin;
  std::ostream &out = ofs.is_open() ? ofs : std::cout;
  std::cerr << "#read  file: " << (!ifs.is_open() ? "stdin"  : argv[1]) << std::endl;
  std::cerr << "#write file: " << (!ofs.is_open() ? "stdout" : argv[2]) << std::endl;

  const QHull::Vertex::vector pos = readData(in);

  QHull q;

  q.convexHull(pos);

  const int nFaces = q.getNumFacets();

  dumpGnuPlot(pos, q, out);
  fprintf(stderr, "nFaces= %d\n", nFaces);

  ifs.close();
  ofs.close();
  return 0;
};
