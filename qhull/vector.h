#pragma once

template<typename real_t, int N>
class Vector_t
{
  private:
    std::array<real_t,N> x;
  public:
    Vector_t() {}
    Vector_t(const real_t y) 
    {
      for (int l = 0; l < N; l++)
        x[l] = y;
    }
    Vector_t(const real_t *y) 
    {
      for (int l = 0; l < N; l++)
        x[l] = y[l];
    }
    real_t& operator[](const int i)       { return x[i]; }
    real_t  operator[](const int i) const { return x[i]; }

    /******************/
    friend real_t dot(const Vector_t &a, const Vector_t &b)
    {
      real_t sum = 0;
      for (int l = 0; l < N; l++)
        sum += a[l]*b[l];
      return sum;
    }
    friend real_t norm2(const Vector_t &a)
    {
      return dot(a,a);
    }
    friend real_t norm(const Vector_t &a)
    {
      return std::sqrt(norm2(a));
    }
    /******************/
    Vector_t& operator*=(const Vector_t &a)
    {
      for (int l = 0; l < N; l++)
        x[l] *= a[l];
      return *this;
    }
    Vector_t& operator*=(const real_t a)
    {
      for (int l = 0; l < N; l++)
        x[l] *= a;
      return *this;
    }
    friend Vector_t operator*(const Vector_t &a, const Vector_t &b) 
    {
      Vector_t res;
      for (int l = 0; l < N; l++)
        res[l] = a[l]*b[l];
      return res;
    }
    friend Vector_t operator*(const Vector_t &a, const real_t b) 
    {
      Vector_t res;
      for (int l = 0; l < N; l++)
        res[l] = a[l]*b;
      return res;
    }
    friend Vector_t operator*(const real_t a, const Vector_t &b) 
    {
      return b*a;
    }
    /******************/
    Vector_t& operator+=(const Vector_t &a)
    {
      for (int l = 0; l < N; l++)
        x[l] += a[l];
      return *this;
    }
    Vector_t& operator+=(const real_t a)
    {
      for (int l = 0; l < N; l++)
        x[l] += a;
      return *this;
    }
    friend Vector_t operator+(const Vector_t &a, const Vector_t &b) 
    {
      Vector_t res;
      for (int l = 0; l < N; l++)
        res[l] = a[l]+b[l];
      return res;
    }
    friend Vector_t operator+(const Vector_t &a, const real_t b) 
    {
      Vector_t res;
      for (int l = 0; l < N; l++)
        res[l] = a[l]+b;
      return res;
    }
    friend Vector_t operator+(const real_t a, const Vector_t &b) 
    {
      return b+a;
    }
    /******************/
    Vector_t& operator-=(const Vector_t &a)
    {
      for (int l = 0; l < N; l++)
        x[l] -= a[l];
      return *this;
    }
    Vector_t& operator-=(const real_t a)
    {
      for (int l = 0; l < N; l++)
        x[l] -= a;
      return *this;
    }
    friend Vector_t operator-(const Vector_t &a, const Vector_t &b) 
    {
      Vector_t res;
      for (int l = 0; l < N; l++)
        res[l] = a[l]-b[l];
      return res;
    }
    friend Vector_t operator-(const Vector_t &a, const real_t b) 
    {
      Vector_t res;
      for (int l = 0; l < N; l++)
        res[l] = a[l]-b;
      return res;
    }
    friend Vector_t operator+(const real_t a, const Vector_t &b) 
    {
      Vector_t res;
      for (int l = 0; l < N; l++)
        res[l] = a - b[l];
      return res;
    }
};
