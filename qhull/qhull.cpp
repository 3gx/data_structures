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
    const typename Vertex::vector &pos, 
    Simplex &simplex)
{
  real_t xMin = +HUGE, xMax = -HUGE;
  const int np = pos.size();
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

  typedef QHull::real_t real_t;
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


int main(int argc, char*argv[])
{
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
