#pragma once

template<typename real_t, int N>
class VectorNaive_t
{
  private:
    std::array<real_t,N> x;
  public:
    VectorNaive_t() {}
    VectorNaive_t(const real_t y) 
    {
      for (int l = 0; l < N; l++)
        x[l] = y;
    }
    VectorNaive_t(const real_t *y) 
    {
      for (int l = 0; l < N; l++)
        x[l] = y[l];
    }
    real_t& operator[](const int i)       { return x[i]; }
    real_t  operator[](const int i) const { return x[i]; }

    /******************/
    friend real_t dot(const VectorNaive_t &a, const VectorNaive_t &b)
    {
      real_t sum = 0;
      for (int l = 0; l < N; l++)
        sum += a[l]*b[l];
      return sum;
    }
    friend real_t norm2(const VectorNaive_t &a)
    {
      return dot(a,a);
    }
    friend real_t norm(const VectorNaive_t &a)
    {
      return std::sqrt(norm2(a));
    }
    /******************/
    VectorNaive_t& operator*=(const VectorNaive_t &a)
    {
      for (int l = 0; l < N; l++)
        x[l] *= a[l];
      return *this;
    }
    VectorNaive_t& operator*=(const real_t a)
    {
      for (int l = 0; l < N; l++)
        x[l] *= a;
      return *this;
    }
    friend VectorNaive_t operator*(const VectorNaive_t &a, const VectorNaive_t &b) 
    {
      VectorNaive_t res;
      for (int l = 0; l < N; l++)
        res[l] = a[l]*b[l];
      return res;
    }
    friend VectorNaive_t operator*(const VectorNaive_t &a, const real_t b) 
    {
      VectorNaive_t res;
      for (int l = 0; l < N; l++)
        res[l] = a[l]*b;
      return res;
    }
    friend VectorNaive_t operator*(const real_t a, const VectorNaive_t &b) 
    {
      return b*a;
    }
    /******************/
    VectorNaive_t& operator+=(const VectorNaive_t &a)
    {
      for (int l = 0; l < N; l++)
        x[l] += a[l];
      return *this;
    }
    VectorNaive_t& operator+=(const real_t a)
    {
      for (int l = 0; l < N; l++)
        x[l] += a;
      return *this;
    }
    friend VectorNaive_t operator+(const VectorNaive_t &a, const VectorNaive_t &b) 
    {
      VectorNaive_t res;
      for (int l = 0; l < N; l++)
        res[l] = a[l]+b[l];
      return res;
    }
    friend VectorNaive_t operator+(const VectorNaive_t &a, const real_t b) 
    {
      VectorNaive_t res;
      for (int l = 0; l < N; l++)
        res[l] = a[l]+b;
      return res;
    }
    friend VectorNaive_t operator+(const real_t a, const VectorNaive_t &b) 
    {
      return b+a;
    }
    /******************/
    VectorNaive_t& operator-=(const VectorNaive_t &a)
    {
      for (int l = 0; l < N; l++)
        x[l] -= a[l];
      return *this;
    }
    VectorNaive_t& operator-=(const real_t a)
    {
      for (int l = 0; l < N; l++)
        x[l] -= a;
      return *this;
    }
    friend VectorNaive_t operator-(const VectorNaive_t &a, const VectorNaive_t &b) 
    {
      VectorNaive_t res;
      for (int l = 0; l < N; l++)
        res[l] = a[l]-b[l];
      return res;
    }
    friend VectorNaive_t operator-(const VectorNaive_t &a, const real_t b) 
    {
      VectorNaive_t res;
      for (int l = 0; l < N; l++)
        res[l] = a[l]-b;
      return res;
    }
    friend VectorNaive_t operator+(const real_t a, const VectorNaive_t &b) 
    {
      VectorNaive_t res;
      for (int l = 0; l < N; l++)
        res[l] = a - b[l];
      return res;
    }
};

template<typename real_t, int N> using Vector_t = VectorNaive_t<real_t,N>;
