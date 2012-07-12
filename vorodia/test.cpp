template<const int N>
struct Tuple
{
  int tuple[N];

#if 0
  inline Tuple(const int&);
  inline Tuple(const int&, const int&);
  inline Tuple(const int&, const int&, const int&);
#else
  template<>
  ::Tuple<1>::Tuple(const int);
#endif

};

template<>
inline Tuple<1>::Tuple(const int& v0) {tuple[0] = v0;}
template<>
inline Tuple<2>::Tuple(const int& v0, const int &v1) {tuple[0] = v0;tuple[1] = v1;}


#include <iostream>
#include <cstdio>
#include <cassert>

typedef float real;

int main(int argc, char * argv[])
{
  assert(argc > 1);
  const int n = atoi(argv[1]);

  fprintf(stderr, " n = %d\n", n);

  int i = -1;
  while (++i < n)
  {
    fprintf(stderr, "i= %d\n", i);
  }
  return 0;
}
