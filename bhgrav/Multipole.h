#ifndef __MULTIPOLE_H__
#define __MULTIPOLE_H__

struct Monopole
{
  typedef std::vector<Monopole> Vector;
  private:
  dvec3  _mpos;
  double _mass;

  public:

  Monopole() : _mpos(vec3(0.0)), _mass(0.0) {}
  Monopole(const dvec3 &pos, const double mass) :
    _mpos(pos*mass), _mass(mass) {}

  const Monopole& operator+=(const Monopole &m)
  {
    _mpos += m._mpos;
    _mass += m._mass;
    return *this;
  }

  const dvec3& mpos() const {return _mpos;}
  dvec3   pos() const {return _mpos*(1.0/_mass);}
  double mass() const {return _mass;}
};

struct Quadrupole
{
  typedef std::vector<Quadrupole> Vector;
  enum TYPE {UNIT};

  private:
  double _xx, _yy, _zz;
  double _xy, _xz, _yz;

  public:
  Quadrupole() : _xx(0.0), _yy(0.0), _zz(0.0), _xy(0.0), _xz(0.0), _yz(0.0) {}
  Quadrupole(const double x, const TYPE type) : _xx(x), _yy(x), _zz(x), _xy(0.0), _xz(0.0), _yz(0.0) {}
  Quadrupole(const dvec3 &pos, const double mass) :
    _xx(mass*pos.x*pos.x),
    _yy(mass*pos.y*pos.y),
    _zz(mass*pos.z*pos.z),
    _xy(mass*pos.x*pos.y),
    _xz(mass*pos.x*pos.z),
    _yz(mass*pos.y*pos.z) {}

  const Quadrupole& operator+=(const Quadrupole &q)
  {
    _xx += q._xx; _yy += q._yy; _zz += q._zz;
    _xy += q._xy; _xz += q._xz; _yz += q._yz;
    return *this;
  }
  const Quadrupole& operator*=(const double x)
  {
    _xx *= x; _yy *= x; _zz *= x;
    _xy *= x; _xz *= x; _yz *= x;
    return *this;
  }

  double xx() const {return _xx;}
  double yy() const {return _yy;}
  double zz() const {return _zz;}
  double xy() const {return _xy;}
  double xz() const {return _xz;}
  double yz() const {return _yz;}
  double trace() const {return _xx + _yy + _zz;}

};

struct Multipole
{
  typedef std::vector<Multipole> Vector;

  private:

  Monopole     _monopole;
  Quadrupole _quadrupole;

  public:

  Multipole(
      const   Monopole &m =   Monopole(),
      const Quadrupole &q = Quadrupole()) :
    _monopole(m), _quadrupole(q) {}

  Multipole(const dvec3 &pos, const double mass) :
    _monopole(pos, mass), _quadrupole(pos, mass) {}

  const Multipole& operator+=(const Multipole &m)
  {
    _monopole   += m.  _monopole;
    _quadrupole += m._quadrupole;
    return *this;
  }

  const   Monopole   monopole() const {return   _monopole;}
  const Quadrupole quadrupole() const {return _quadrupole;}

  Multipole complete() const 
  {
    Multipole m(*this);
    const double _mass =   _monopole. mass();
    const dvec3   _pos =   _monopole.  pos();
    const double trace = _quadrupole.trace();
    const double P     = trace - _mass*_pos.norm2();
    m._quadrupole += Quadrupole(_pos, -_mass);
    m._quadrupole *= 3.0;
    m._quadrupole += Quadrupole(-P, Quadrupole::UNIT);

    return m;
  }
};

  template<const bool ROOT>
Multipole computeMultipole(const Particle::Vector &ptcl, const int addr = 0)
{
  Multipole multipole;
  if (ROOT)
  {
    assert(isTreeReady());
    multipoleList.resize(ncell);
    cellCoM      .resize(ncell);
    for (int k = 0; k < 8; k++)
      if (!cellList[k].isEmpty() && cellList[k].isTouched())
        multipole += computeMultipole<false>(ptcl, k);
    multipole = multipole.complete();
  }
  else
  {
    Cell &cell = cellList[addr];
    assert(!cell.isEmpty());
    if (cell.isNode())
    {
      const int i = cell.addr();
      for (int k = 0; k < 8; k++)
        if (!cellList[i+k].isEmpty() && cellList[i+k].isTouched())
          multipole += computeMultipole<false>(ptcl, i+k);
    }
    else
    {
      const Leaf &leaf = leafList[cell.leafIdx()];
      for (int i = 0; i < leaf.nb(); i++)
      {
        const int idx = leaf[i].idx();
        assert(idx < (int)ptcl.size());
        multipole += Multipole(ptcl[idx].pos, ptcl[idx].mass);
        assert(overlapped(bndsList[cell.id()].inner(), leaf[i ].vector_pos()));
        assert(overlapped(bndsList[cell.id()].inner(), ptcl[idx].pos));
      }
      const vec3 com = multipole.monopole().pos();
    }
    multipoleList[cell.id()] = multipole.complete();

    /* compute opening criterion for this cell */

    const int id = cell.id();
    const Multipole &m = multipoleList[id];
    const vec3 com = m.monopole().pos();
    const float  s = (bndsList[id].inner().center() - com).abs();
    const vec3 len =  bndsList[id].inner().  hlen();
    const float  l = __max(len.x, __max(len.y, len.z)) * 2.0f;
    cellCoM[id] = float4(com.x, com.y, com.z, SQR(l*inv_theta + s));

#if 0
    fprintf(stderr, "dcom= %f %f %f  c= %f %f %f  l=%G\n",
        com.x, com.y, com.z,
        bndsList[id].inner().center().x,
        bndsList[id].inner().center().y,
        bndsList[id].inner().center().z, l);

    static int cnt = 0;
    if (cnt++ > 1000) exit(0);
#endif

    /* now the physical properties of the cell have been updated, unTouch it */

    cell.unsetTouched();
  }
  return multipole;
}

#endif /* __MULTIPOLE_H__ */
