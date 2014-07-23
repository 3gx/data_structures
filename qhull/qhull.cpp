#include "qhull.h"

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


int main(int argc, char*argv[])
{
  QHull q;
  return 0;
};
