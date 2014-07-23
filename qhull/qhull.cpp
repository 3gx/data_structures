#include "qhull.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>


#define N 2
using QHull = QHull_t<N>;


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

void dumpGnuPlot(const QHull::Vertex::vector &pos, const QHull &q, std::ostream &ofs)
{
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
