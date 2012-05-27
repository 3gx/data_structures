#ifndef __VECTOR3_H
#define __VECTOR3_H

#include <iostream>
#include <fstream>
#include <cmath>

struct float4
{
#ifdef __mySSEX__
  v4sf  v;
  float4(const v4sf _v) : v(_v) {}
  float4(const float _x, const float _y, const float _z, const float _w) 
  {
    v = (v4sf){_x, _y, _z, _w};
  }
  operator v4sf() const {return v;}
  float4 operator-(const float4 rhs) const
  {
    return (v4sf)rhs - v;
  }
  float norm2() const
  {
    const v4si mask = {0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x0};
    const v4sf r    = __builtin_ia32_andps(v, (v4sf)mask);
    const v4sf r2   = r*r;
    const v4sf tmp  = __builtin_ia32_haddps(r2,  r2);
    const v4sf res  = __builtin_ia32_haddps(tmp, tmp);
    return __builtin_ia32_vec_ext_v4sf(res, 0);
  }
  float x() const {return __builtin_ia32_vec_ext_v4sf(v, 0);}
  float y() const {return __builtin_ia32_vec_ext_v4sf(v, 1);}
  float z() const {return __builtin_ia32_vec_ext_v4sf(v, 2);}
  float w() const {return __builtin_ia32_vec_ext_v4sf(v, 3);}
#else
  float _x, _y, _z, _w;
  float4(const float x, const float y, const float z, const float w) 
  {
    _x = x;
    _y = y;
    _z = z;
    _w = w;
  }
  float4 operator-(const float4 v) const
  {
    return float4(_x-v._x, _y-v._y, _z-v._z, _w-v._w);
  }
  float norm2() const
  {
    return _x*_x+_y*_y+_z*_z;
  }
  float x() const {return _x;}
  float y() const {return _y;}
  float z() const {return _z;}
  float w() const {return _w;}
#endif
  float4() {}
};


template <class REAL> struct vector3{
public:
#ifdef __mySSE__
  union
  {
    struct{ REAL x, y, z, w; };
    v4sf v;
  };
  vector3 (const float4 r) 
  {
    assert(sizeof(float)  == sizeof(REAL));
    x = r.x();
    y = r.y();
    z = r.z();
    w = r.w();
  }
  vector3 (const v4sf _v) : v(_v) {}
  operator v4sf() const {return v;}
#else
  REAL x, y, z, w;
  vector3 (const float4 r) : x(r.x()), y(r.y()), z(r.z()), w(r.w())  {}
#endif

	vector3(){
		x = y = z = REAL(0);
    w = REAL(0);
	}
	vector3(const REAL &r){
		x = y = z = r;
    w = REAL(0);
	}
	vector3(const REAL &_x, const REAL &_y, const REAL &_z){
		x = _x;  y = _y;  z = _z;
    w = REAL(0);
	}
	vector3(const REAL *p){
		x = p[0]; y = p[1]; z = p[2];
    w = REAL(0);
	}
	~vector3(){}

	REAL &operator [](int i){
	  return (&x)[i];
	}
	const REAL &operator [](int i) const{
	  return (&x)[i];
	}
	template <class real> 
		operator vector3<real> () const {
			return vector3<real> (real(x), real(y), real(z));
		}
	operator REAL *(){
		return &x;
	}
	REAL (*toPointer())[3]{
		return (REAL (*)[3])&x;
	}
	typedef REAL (*pArrayOfReal3)[3];
	operator pArrayOfReal3(){
		return toPointer();
	}

	void outv(std::ostream &ofs = std::cout) const{
		ofs << "(" << x << ", " << y << ", " << z << ")" << std::endl;
	}
	bool are_numbers () const{
		// returns false if *this has (a) NaN member(s)
		return (norm2() >= REAL(0));
	}

	REAL norm2() const{
		return (*this)*(*this);
	}
	REAL abs() const{
		return std::sqrt(norm2());
	}

	friend std::ostream &operator << (std::ostream &ofs, const vector3<REAL> &v){
		ofs << v.x << " " << v.y << " " << v.z;
		return ofs;
	}
	friend std::istream &operator >> (std::istream &ifs, vector3<REAL> &v){
		ifs >> v.x >> v.y >> v.z;
		return ifs;
	}
	const vector3<REAL> operator + (const vector3<REAL> &v) const{
		return vector3<REAL> (x+v.x, y+v.y, z+v.z);
	}
	const inline vector3<REAL> operator - (const vector3<REAL> &v) const{
		return vector3<REAL> (x-v.x, y-v.y, z-v.z);
	}
	const vector3<REAL> operator * (const REAL &s) const{
		return vector3<REAL> (x*s, y*s, z*s);
	}
	friend const vector3<REAL> operator * (const REAL &s, const vector3<REAL> &v){
		return v*s;
	}
	// dot product
	const inline REAL operator * (const vector3<REAL> &v) const{
		return (x*v.x + y*v.y + z*v.z);
	}
	// vector product
	const vector3<REAL> operator % (const vector3<REAL> &v) const{
		return vector3<REAL> (
        y*v.z - z*v.y, 
        z*v.x - x*v.z, 
        x*v.y - y*v.x);
  }
  const vector3<REAL> operator / (const REAL &s) const{
    REAL r = REAL(1)/s;
    return (*this)*r;
  }
  const vector3<REAL> operator = (const vector3<REAL> &v){
    x = v.x; y=v.y; z=v.z;
    return *this;
  }

  const vector3<REAL> operator - (){
    return vector3<REAL> (-x, -y, -z);
  }
  const vector3<REAL> &operator += (const vector3<REAL> &v){
    *this = *this + v;
    return *this;
  }
  const vector3<REAL> &operator -= (const vector3<REAL> &v){
    *this = *this - v;
    return *this;
  }
  const vector3<REAL> &operator *= (const REAL &s){
    *this = *this * s;
    return *this;
  }
  const vector3<REAL> &operator /= (const REAL &s){
    *this = *this / s;
    return *this;
  }

  friend const vector3<REAL> maxeach (const vector3<REAL> &a, const vector3<REAL> &b){
    return vector3<REAL> (std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z));
  }
  friend const vector3<REAL> mineach (const vector3<REAL> &a, const vector3<REAL> &b){
    return vector3<REAL> (std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z));
  }
  const vector3<REAL> abseach(){
    return vector3<REAL> (std::fabs(x), std::fabs(y), std::fabs(z));
  }
};

typedef vector3<double> dvec3;
typedef vector3<float>  fvec3;



#endif //  __VECTOR3_H
