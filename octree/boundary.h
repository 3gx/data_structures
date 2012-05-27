#ifndef __BOUNDARY_H__
#define __BOUNDARY_H__

#include <vector>
#include "vector3.h"
#include <cmath>

template <typename REAL>
struct Boundary{
  typedef std::vector<Boundary> Vector;
	typedef vector3<REAL> vec;

	vec min, max;
	Boundary() : min(HUGE), max(-HUGE) {}
  Boundary(const float4 pos) 
  {
    min = vec(pos.x() - pos.w(), pos.y() - pos.w(), pos.z() - pos.w());
    max = vec(pos.x() + pos.w(), pos.y() + pos.w(), pos.z() + pos.w());
  }
	Boundary(const vec &_min, const vec &_max) : min(_min), max(_max) {}
	Boundary(const vec &pos, const REAL &h = 0.0) : min(pos - vec(h)), max(pos + vec(h)) {}
	
  const vec center() const {
		return REAL(0.5) * (max + min);
	}
	const vec hlen() const {
		return REAL(0.5) * (max - min);
	}
	const REAL separation2_from(const vec &pos) const{
		vec dr = center() - pos;
		dr = dr.abseach() - hlen();
		dr = vec::maxeach(dr, vec(0.0));
		return dr.norm2();
	}

#ifdef __mySSE__
	Boundary(v4sf _min, v4sf _max)
  {
		min = _min;
		max = _max;
	}

	static const Boundary merge(const Boundary &a, const Boundary &b){
		return Boundary(
				__builtin_ia32_minps(a.min, b.min),
				__builtin_ia32_maxps(a.max, b.max));
	}
	void merge(const Boundary &b){
		*this = merge(*this, b);
	}
	friend bool not_overlapped(const Boundary &a, const Boundary &b)
  {
		return __builtin_ia32_movmskps(__builtin_ia32_orps(
         __builtin_ia32_cmpltps(a.max, b.min),
         __builtin_ia32_cmpltps(b.max, a.min)
         ));
	}
	friend bool overlapped(const Boundary &a, const Boundary &b){
		return !not_overlapped(a, b);
	}
#else
	static const Boundary merge(const Boundary &a, const Boundary &b){
		return Boundary(mineach(a.min, b.min), maxeach(a.max, b.max));
	}
	void merge(const Boundary &b){
		*this = merge(*this, b);
	}
	friend bool not_overlapped(const Boundary &a, const Boundary &b){
		return (a.max.x < b.min.x) || (b.max.x < a.min.x)
		    || (a.max.y < b.min.y) || (b.max.y < a.min.y)
		    || (a.max.z < b.min.z) || (b.max.z < a.min.z);
	}
	friend bool overlapped(const Boundary &a, const Boundary &b){
		return !not_overlapped(a, b);
	}
#endif
};
#if 0
template <>
struct Boundary<float>{
	typedef vector3<float> vec;
	typedef float v4sf __attribute__ ((vector_size(16)));

	vec min; 
	float p0;
	vec max;
	float p1;

	Boundary() : min(HUGE), p0(0.f), max(-HUGE), p1(0.f) {}
	Boundary(const vec &_min, const vec &_max) : min(_min), p0(0.f), max(_max), p1(0.f) {}
	Boundary(const vec &pos, float h = 0.f) : 
		min(pos - vec(h)), p0(0.f), max(pos + vec(h)), p1(0.f) {}
	Boundary(v4sf _min, v4sf _max){
		*(v4sf *)&min = _min;
		*(v4sf *)&max = _max;
	}

	static const Boundary merge(const Boundary &a, const Boundary &b){
		return Boundary(
				__builtin_ia32_minps(*(v4sf *)&a.min, *(v4sf *)&b.min),
				__builtin_ia32_maxps(*(v4sf *)&a.max, *(v4sf *)&b.max));
	}
	void merge(const Boundary &b){
		*this = merge(*this, b);
	}
	friend bool not_overlapped(const Boundary &a, const Boundary &b){
		return __builtin_ia32_movmskps(
				(v4sf)(__builtin_ia32_cmpltps(
						*(v4sf *)&a.max, *(v4sf *)&b.min)))
		   ||  __builtin_ia32_movmskps(
				(v4sf)(__builtin_ia32_cmpltps(
						*(v4sf *)&b.max, *(v4sf *)&a.min)));
	}
	friend bool overlapped(const Boundary &a, const Boundary &b){
		return !not_overlapped(a, b);
	}
};
#endif

#endif /* __BOUNDARY_H__ */
