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
    friend VectorNaive_t operator-(const real_t a, const VectorNaive_t &b) 
    {
      VectorNaive_t res;
      for (int l = 0; l < N; l++)
        res[l] = a - b[l];
      return res;
    }
};


/*===================================================================*/
template <typename real_t,int N> class VectorET_t;
template<typename LHS,typename OP,typename RHS>
class mathVecExpr
{
  private:
    mathVecExpr _l, _r;
  public:
    mathVecExpr() = delete;
    mathVecExpr(LHS l,RHS r) : 
      _l(std::forward<LHS>(l)),
      _r(std::forward<RHS>(r)) { }

    /* prohibit copying */
    mathVecExpr(mathVecExpr const&) = delete;
    mathVecExpr& operator=(mathVecExpr const&) = delete;

    /* allow moves */
    mathVecExpr(mathVecExpr&&) = default;
    mathVecExpr& operator=(mathVecExpr&&) = default;

    template<typename RE>
      auto operator+(RE&& re) const ->
      mathVecExpr<mathVecExpr<LHS,OP,RHS> const&,OP,decltype(std::forward<RE>(re))>
      {
        return
          mathVecExpr<
            mathVecExpr<LHS,OP,RHS> const&,
            OP,
            decltype(std::forward<RE>(re))>(*this, std::forward<RE>(re));
      }

    auto le() -> typename std::add_lvalue_reference<LHS>::type { return _l; }
    auto le() const -> typename std::add_lvalue_reference<
                  typename std::add_const<LHS>::type>::type
                  {
                    return _l;
                  }
    auto re() -> typename std::add_lvalue_reference<RHS>::type { return _r; }
    auto re() const -> typename std::add_lvalue_reference<
                  typename std::add_const<RHS>::type>::type
                  {
                    return _r;
                  }

    auto operator[](const int index) const ->
      decltype(OP::apply(this->le()[index], this->re()[index]))
      {
        return OP::apply(le()[index], re()[index]);
      }

};
template <typename T>
struct PlusOp
{
  static T apply(T const& a, T const& b)
  {
    return a + b;
  }

  static T apply(T&& a, T const& b)
  {
    a += b;
    return std::move(a);
  }

  static T apply(T const& a, T&& b)
  {
    b += a;
    return std::move(b);
  }

  static T apply(T&& a, T&& b)
  {
    a += b;
    return std::move(a);
  }
};
/*===================================================================*/

template<typename real_t, int N> using Vector_t = VectorNaive_t<real_t,N>;
