#pragma once

#include <array>
  
template<typename real_t, int N>
static 
std::array<std::array<real_t,N>,N>
choleskyInverse(const std::array<std::array<real_t,N>,N> &in)
{
  std::array<std::array<real_t,N>,N> out;
  /* extract upper triangle from the [symmetric] input matrix */

  for (int j = 0; j < N; j++)
    for (int i = j; i < N; i++)
      out[j][i] = out[i][j] = in[j][i];

  /* Cholesky factorization, stored in the lower triangle */

  for (int k = 0; k < N; k++)
  {
    if (out[k][k] <= 0) assert(0);  /* matrix is not positive definite */
    out[k][k] = std::sqrt(out[k][k]);
    const real_t ainv = real_t(1.0)/out[k][k];
    for (int i = k+1; i < N; i++)
      out[i][k] *= ainv;
    for (int j = k+1; j < N; j++)
      for (int i = j; i < N; i++)
        out[i][j] -= out[i][k]*out[j][k];
  }

  for (int j = 0; j < N; j++)
    for (int i = j+1; i < N; i++)
      out[j][i] = 0.0;
  
  /* determinant */
  
  real_t det = 1.0;
  for (int i = 0; i < N; i++)
    det *= out[i][i];
  det *= det;
  assert(det > 0.0);

  /* invert lower triangular matrix */

  for (int k = 0; k < N; k++)
    out[k][k] = real_t(1.0)/out[k][k];

  for (int i = 1; i < N; i++)
    for (int j = 0; j < i; j++)
    {
      real_t sum = 0.0;
      for (int k = j; k < i; k++)
        sum += out[i][k] * out[k][j];
      out[i][j] = -out[i][i]*sum;
    }

  /* compute inverse by multiplying inverse of L and its transpose */

  for (int j = 0; j < N; j++)
    for (int i = j; i < N; i++)
    {
      real_t sum = 0;
      for (int k = i; k < N; k++)
        sum += out[k][j]*out[k][i];
      out[j][i] = sum;
    }

  for (int j = 0; j < N; j++)
    for (int i = j+1; i < N; i++)
      out[i][j] = out[j][i];


  return out;
}

  template<typename real_t, int D>
static std::array<real_t,D> linSolve(const std::array<std::array<real_t,D>,D> &m, const std::array<real_t,D> &b)
{
  const auto mInv = choleskyInverse<real_t,D>(m);
  std::array<real_t,D> x;
  for (int l = 0; l < D; l++)
  {
    x[l] = 0;
    for (int ll = 0; ll < D; ll++)
      x[l] += mInv[l][ll]*b[ll];
  }
  return x;
}
