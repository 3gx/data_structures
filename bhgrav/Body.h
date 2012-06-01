#ifndef __BODY_H__
#define __BODY_H__
  
struct Body
{
  typedef std::vector<Body> Vector;
  private:
  float4 packed_pos;
  int _idx;
  float _mass;
  int   iPad2, iPad3; /* padding */

  public:
  Body() 
  {
    assert(sizeof(Body) == sizeof(float)*8);
  }
  Body(const vec3 &pos, const int idx, const float h = 0.0f, const float m = 0.0f) : packed_pos(float4(pos.x, pos.y, pos.z, h)), _idx(idx), _mass(m) {}
  Body(const Particle &p, const int idx) : packed_pos(float4(p.pos.x, p.pos.y, p.pos.z, p.h)), _idx(idx), _mass(p.mass) {}
  Body(const vec3 &pos, const float h, const int idx = -1, const float mass = 0.0f) : packed_pos(pos.x, pos.y, pos.z, h), _idx(idx), _mass(mass) {};
  int           idx() const {return _idx;}
  float        mass() const {return _mass;}
  float4      pos_h() const {return packed_pos;}
  float           h() const {return packed_pos.w();}
  float4   pos_mass() const {return float4(packed_pos.x(), packed_pos.y(), packed_pos.z(), _mass); }
  vec3   vector_pos() const {return vec3(packed_pos.x(), packed_pos.y(), packed_pos.z());}
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
