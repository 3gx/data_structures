#include <iostream>
#include <cstdio>
#include <cassert>

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
