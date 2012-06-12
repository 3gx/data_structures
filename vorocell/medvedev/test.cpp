#include <algorithm>
#include <vector>
#include <stack>
#include <cstdio>
#include <cstdlib>
#include <cassert>
// Our template functor
template <typename T1, typename T2>
struct t_unpair
{
  T1& a1;
  T2& a2;
  explicit t_unpair( T1& a1, T2& a2 ): a1(a1), a2(a2) { }
  t_unpair<T1,T2>& operator = (const std::pair<T1,T2>& p)
  {
    a1 = p.first;
    a2 = p.second;
    return *this;
  }
};

// Our functor helper (creates it)
  template <typename T1, typename T2>
inline t_unpair<T1,T2> unpair( T1& a1, T2& a2 )
{
  return t_unpair<T1,T2>( a1, a2 );
}


std::pair<int, int> test(const int n)
{
  return std::make_pair(n, n*n);
}


int main(int argc, char * argv[])
{
  assert(argc > 1);
  const int n = atoi(argv[1]);
  fprintf(stderr, " n = %d\n", n);

  std::vector<double> vec;
  std::stack <double> stk;

  for (int i = 0; i < n; i++)
  {
    vec.push_back(drand48());
    stk.push(vec.back());
  }

  std::sort(vec.begin(), vec.end());

  int n1, n2;
  unpair(n1, n2) = test(n);
  fprintf(stderr, "n1= %d n2= %d \n", n1, n2);

  return 0;
}
