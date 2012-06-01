#ifndef __BODY_H__
#define __BODY_H__
  
struct Body
{
  typedef std::vector<Body> Vector;
  private:
  float4 Packed_pos;
  int _idx;
  int iPad1, iPad2, iPad3; /* padding */

  public:
  Body() 
  {
    assert(sizeof(Body) == sizeof(float)*8);
  }
  Body(const vec3 &pos, const int idx, const float h = 0.0f) : Packed_pos(float4(pos.x, pos.y, pos.z, h)), _idx(idx) {}
  Body(const Particle &p, const int idx) : Packed_pos(float4(p.pos.x, p.pos.y, p.pos.z, p.h)), _idx(idx) {}
  Body(const vec3 &pos, const float h, const int idx = -1) : Packed_pos(pos.x, pos.y, pos.z, h), _idx(idx) {};
  int idx() const {return _idx;}
  float4 packed_pos() const {return Packed_pos;}
  float          h()  const {return Packed_pos.w();}
  vec3   vector_pos() const {return vec3(Packed_pos.x(), Packed_pos.y(), Packed_pos.z());}
#ifdef __SSE_H__
  Body& operator=(const Body &rhs)
  {
    _v4sf *dst = (_v4sf*)this;
    _v4sf *src = (_v4sf*)&rhs;
#if 1
    dst[0] = src[0];
    dst[1] = src[1];
#else
    const v4sf a = src[0];
    const v4sf b = src[1];
    __builtin_ia32_movntps((float*)(dst+0), a);
    __builtin_ia32_movntps((float*)(dst+1), b);
#endif
    return *this;
  }
#endif
};

#endif /* __BODY_H__ */
