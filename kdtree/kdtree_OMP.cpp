#include "kdtree.h"

class CompareBodiesOMP
{
  unsigned int _split_dim;
  public:
  CompareBodiesOMP(const unsigned int val) : _split_dim(val%3) {};
  bool operator() (const kdBody &lhs, const kdBody &rhs) const
  {
    return lhs.pos()[_split_dim] < rhs.pos()[_split_dim];
  }
};

void kdTree::nth_element_omp(
    kdBody::Iterator beg, 
    kdBody::Iterator med,  
    kdBody::Iterator end,
    const int split_dim)
{
  std::nth_element(beg, med, end, CompareBodiesOMP(split_dim));
}
