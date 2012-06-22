#include <cstdio>
#include <cstdlib>

int main()
{
  int n = 2000;
  for (int i = 0; i < n; i++)
  {
    fprintf(stdout, " %d %g %g %g \n",
        i+1,
        drand48(),
        drand48(),
        drand48());
  }
  return 0;
}
