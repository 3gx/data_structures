#pragma once

  template<typename real_t, int D>
static inline std::array<real_t,D> linSolve(const std::array<std::array<real_t,D>,D> &m, const std::array<real_t,D> &b)
{
  std::array<real_t,D> x = b;
  return x;
}
