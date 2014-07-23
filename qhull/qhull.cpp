#include "qhull.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#if 0
template<> 
void 
QHull::extremeSimplexR<2>::eval(const QHull::pos_t::vector &pos, QHull::Simplex_t &simplex)
{
  real_t xMin = +HUGE, xMax = -HUGE;
  const int np = pos.size();
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
#endif

QHull::pos_t::vector readData(std::istream &in)
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
  typedef QHull::pos_t  pos_t;
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
  pos_t::vector pos(count);
  for (size_t i = 0; i < count; ++i) 
  {
    if (!std::getline(in, line)) {
      std::cerr << "io: line count error" << std::endl;
      exit(-1);
    }
    pos_t &p = pos[i];
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
  std::istream &in = ifs.is_open() ? ifs : std::cin;
  std::cerr << "#read file: " << (ifs.is_open() ? "stdin" : argv[1]) << std::endl;

  const QHull::pos_t::vector pos = readData(in);

  QHull q;
  return 0;
};
