#ifndef __SSE_H__
#define __SSE_H__

typedef int    _v4si  __attribute__((vector_size(16)));
typedef float  _v4sf  __attribute__((vector_size(16)));

inline _v4sf v4sf(const float x) {return (_v4sf){x,x,x,x};}
inline _v4si v4si(const  int  x) {return (_v4si){x,x,x,x};}

#endif /* __SSE_H__ */
