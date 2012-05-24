#include <iostream>
#include <cassert>
#include <cstdio>
#include <cstdlib>
  
  template <class T>
T prevPow2(const T v)
{
  int k;
  if (v == 0)
    return 1;
  for (k = sizeof(T) * 8 - 1; ((static_cast<T>(1U) << k) & v) == 0; k--);
  return static_cast<T>(1U) << k;
}

  template <class T>
T nearPow2(const T v)
{
  int k;
  if (v == 0)
    return 1;
  for (k = sizeof(T) * 8 - 1; ((static_cast<T>(1U) << k) & v) == 0; k--);
  if (((static_cast<T>(1U) << (k - 1)) & v) == 0)
    return static_cast<T>(1U) << k;
  return static_cast<T>(1U) << (k + 1);
}
  
  template <class T>
T nextPow2(const T v)
{
  int k;
  if (v == 0)
    return 1;
  for (k = sizeof(T) * 8 - 1; ((static_cast<T>(1U) << k) & v) == 0; k--);
  return static_cast<T>(1U) << (k + 1);
}

int main(int  argc, char * argv[])
{
  assert(argc > 1);
  const int n = atoi(argv[1]);

  fprintf(stderr, "n= %d:  low2= %d  near2= %d high2= %d\n",
      n, prevPow2(n), nearPow2(n), nextPow2(n));
}


