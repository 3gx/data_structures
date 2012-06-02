#ifndef __GFORCE_GROUP_H__
#define __GFORCE_GROUP_H__

inline bool split_node(
    const float4 &cellCoM,     /* contains (x,y,z,w), where w is an opening criteria distance */
    const float4 &groupCentre,
    const float4 &groupSize) const

{
  const float4  dr = fabs(groupCentre - cellCoM) - groupSize;
  const float4 adr = dr + fabs(dr);  /* 0 if dr < 0, otherwise 2*dr */
  const float  ds2 = adr.norm2()*0.25f;

  return ds2 <= cellCoM.w();
}

  template<const int Ng>
std::pair<unsigned int, unsigned int> gForce(const GroupT<Ng> &group, float4 force[Ng]) const
{
  const int Ncell = NLEAF;
  const int Nptcl = NLEAF*16;
  float4      ptcl_list[Nptcl*2] __attribute__ ((aligned(32)));
  fMultipole  cell_list[Ncell*2] __attribute__ ((aligned(32)));
  int nc = 0, np = 0;
  int np_tot = 0, nc_tot = 0;

  gForce<true, Ncell, Nptcl>(group, force, cell_list, nc, ptcl_list, np, np_tot, nc_tot);

  return std::make_pair(np_tot*group.nb(), nc_tot*group.nb());
}

template<const bool ROOT, const int Nc, const int Np, const int Ng>
void gForce(
    const GroupT<Ng> &group,
    float4     force[Ng  ],         /* maximal number of particles in a group  */
    fMultipole cell_list[Nc*2], int &nc, /* list for particle-cell interactions     */
    float4     ptcl_list[Np*2], int &np, /* list for particle-particle interactions */
    int &np_tot, int &nc_tot,
    const  float4 &groupCentre = 0.0, const float4 &groupSize = 0.0, const int addr = 0) const
{
  if (ROOT)
  {
    assert(isTreeReady());
    assert(Np >= 8*NLEAF);
    assert(Nc >  8);
    const boundary bnd = group.outerBoundary();
    const float4 groupCentre =  bnd.center().           to_float4();
    const float4 groupSize   = (bnd.hlen  ()*(real)2.0).to_float4();
    nc = np = 0;
    nc_tot = np_tot = 0;
    for (int i = 0; i < Ng; i++)
      force[i] = 0.0f;
    for (int k = 0; k < 8; k++)
      if (!cellList[k].isEmpty())
        if (split_node(cellCoM[cellList[k].id()], groupCentre, groupSize))
          gForce<false, Nc, Np>(group, force, cell_list, nc, ptcl_list, np, np_tot, nc_tot,
              groupCentre, groupSize, k);

    np = particle_particle<Np>(group, force, ptcl_list, np);
    nc = particle_cell    <Nc>(group, force, cell_list, nc);
  }
  else
  {
    const Cell cell = cellList[addr];
    if (split_node(cellCoM[cell.id()], groupCentre, groupSize))
    {
      if (cell.isNode())
      {
        for (int k = 0; k < 8; k++)
          if (!cellList[cell.addr()+k].isEmpty())
            gForce<false, Nc, Np>(group, force, cell_list, nc, ptcl_list, np, np_tot, nc_tot,
                groupCentre, groupSize, cell.addr()+k);
      }
      else
      {
        const Leaf &leaf = leafList[cell.leafIdx()];
        for (int i = 0; i < leaf.nb(); i++)
          ptcl_list[np++] = leaf[i].pos_mass();
        assert(np <= 2*Np);
        np_tot += leaf.nb();
      }
    }
    else
    {
      cell_list[nc++] = multipoleList[cell.id()];
      nc_tot++;
      assert(nc <= 2*Nc);
    }

    if (np >= Np) np = particle_particle<Np>(group, force, ptcl_list, np);
    if (nc >= Nc) nc = particle_cell    <Nc>(group, force, cell_list, nc);
  }
}

template<const int Np, const int Ng>
int particle_particle(
    const  GroupT<Ng> &igroup,
    float4     force[Ng  ],       
    float4 ptcl_list[Np*2], int np) const
{
#ifdef __AVX_H__
  const GroupT<Ng> group = igroup;
  const _v8sf *ib = (const _v8sf*)&group[0];
  const _v4sf *jb = (const _v4sf*)ptcl_list;
  
  const int ni = group.nb();
  const int nj = np;

  const _v8sf veps2 = v8sf(eps2);
  for (int i = 0; i < ni; i += 8)
  {
    const _v8sf ip04 = __mergelo(*(ib+i+0), *(ib+i+4));
    const _v8sf ip15 = __mergelo(*(ib+i+1), *(ib+i+5));
    const _v8sf ip26 = __mergelo(*(ib+i+2), *(ib+i+6));
    const _v8sf ip37 = __mergelo(*(ib+i+3), *(ib+i+7));
        
    const _v8sf xy02 = __builtin_ia32_unpcklps256(ip04, ip26);
    const _v8sf xy13 = __builtin_ia32_unpcklps256(ip15, ip37);
    const _v8sf zw02 = __builtin_ia32_unpckhps256(ip04, ip26);
    const _v8sf zw13 = __builtin_ia32_unpckhps256(ip15, ip37);
    const _v8sf  ipx = __builtin_ia32_unpcklps256(xy02, xy13);
    const _v8sf  ipy = __builtin_ia32_unpckhps256(xy02, xy13);
    const _v8sf  ipz = __builtin_ia32_unpcklps256(zw02, zw13);

    _v8sf fx   = v8sf(0.0f);
    _v8sf fy   = v8sf(0.0f);
    _v8sf fz   = v8sf(0.0f);
    _v8sf gpot = v8sf(0.0f);
    for (int j = 0; j < nj; j++)
    {
      const _v8sf jp  = pack2ymm(*(jb+j), *(jb+j));
      const _v8sf jpx = __bcast<0>(jp);
      const _v8sf jpy = __bcast<1>(jp);
      const _v8sf jpz = __bcast<2>(jp);
      const _v8sf jpm = __bcast<3>(jp);

      const _v8sf dx = jpx - ipx;
      const _v8sf dy = jpy - ipy;
      const _v8sf dz = jpz - ipz;
      const _v8sf r2 = dx*dx + dy*dy + dz*dz + veps2;

      const _v8sf  rinv  = __builtin_ia32_rsqrtps256(r2);
      const _v8sf mrinv  = rinv  *  jpm;
      const _v8sf  rinv2 = rinv  *  rinv;
      const _v8sf mrinv3 = rinv2 * mrinv;

      fx   += mrinv3 * dx;
      fy   += mrinv3 * dy;
      fz   += mrinv3 * dz;
      gpot += mrinv;
    }
    gpot = -gpot;

    const _v8sf t0  = __builtin_ia32_unpcklps256(fx,   fz);
    const _v8sf t1  = __builtin_ia32_unpcklps256(fy, gpot);
    const _v8sf t2  = __builtin_ia32_unpckhps256(fx,   fz);
    const _v8sf t3  = __builtin_ia32_unpckhps256(fy, gpot);
    const _v8sf f04 = __builtin_ia32_unpcklps256(t0, t1);
    const _v8sf f15 = __builtin_ia32_unpckhps256(t0, t1);
    const _v8sf f26 = __builtin_ia32_unpcklps256(t2, t3);
    const _v8sf f37 = __builtin_ia32_unpckhps256(t2, t3);

    _v8sf* vforce = (_v8sf*)&force[i];
    *(vforce + 0) += __merge<0,0>(f04, f15);
    *(vforce + 1) += __merge<0,0>(f26, f37);
    *(vforce + 2) += __merge<1,1>(f04, f15);
    *(vforce + 3) += __merge<1,1>(f26, f37);
  }
#else
  const  GroupT<Ng> &group = igroup;
  const int ni = group.nb();
  const int nj = np;
  for (int i = 0; i < ni; i++)
  {
    const float4 ip = group[i].pos_mass();
    for (int j = 0; j < nj; j++)
    {
      const float4 jp = ptcl_list[j];
      const float4 dr = jp - ip;
      const float  r2 = dr.norm2() + eps2;
      const float  mj = jp.w();

      const float  rinv  = 1.0f/std::sqrt(r2);
      const float mrinv  = rinv*mj;
      const float mrinv3 = rinv*rinv*mrinv;

      float4 acc = dr * float4(mrinv3);
      acc.w() = -mrinv;

      force[i] = force[i] + acc;
    }
  }
#endif

  return 0;
}

template<const int Nc, const int Ng>
int particle_cell(
    const GroupT<Ng> &igroup,
    float4     force[Ng  ],       
    fMultipole cell_list[Nc*2], int nc) const
{
#ifdef __AVX_H__
  const GroupT<Ng> group = igroup;
  const _v8sf *ib = (const _v8sf*)&group[0];
  const _v4sf *jb = (const _v4sf*)cell_list;
  
  const int ni  = group.nb();
  const int nj  = nc;
  const int nj3 = nj * 3;
  
  const _v8sf veps2 = v8sf(eps2);
  for (int i = 0; i < ni; i += 8)
  {
    const _v8sf ip04 = __mergelo(*(ib+i+0), *(ib+i+4));
    const _v8sf ip15 = __mergelo(*(ib+i+1), *(ib+i+5));
    const _v8sf ip26 = __mergelo(*(ib+i+2), *(ib+i+6));
    const _v8sf ip37 = __mergelo(*(ib+i+3), *(ib+i+7));
        
    const _v8sf xy02 = __builtin_ia32_unpcklps256(ip04, ip26);
    const _v8sf xy13 = __builtin_ia32_unpcklps256(ip15, ip37);
    const _v8sf zw02 = __builtin_ia32_unpckhps256(ip04, ip26);
    const _v8sf zw13 = __builtin_ia32_unpckhps256(ip15, ip37);
    const _v8sf  ipx = __builtin_ia32_unpcklps256(xy02, xy13);
    const _v8sf  ipy = __builtin_ia32_unpckhps256(xy02, xy13);
    const _v8sf  ipz = __builtin_ia32_unpcklps256(zw02, zw13);

    _v8sf fx   = v8sf(0.0f);
    _v8sf fy   = v8sf(0.0f);
    _v8sf fz   = v8sf(0.0f);
    _v8sf gpot = v8sf(0.0f);
    for (int j = 0; j < nj3; j += 3)
    {
      const _v8sf jp  = pack2ymm(*(jb+j), *(jb+j));
      const _v8sf jpx = __bcast<0>(jp);
      const _v8sf jpy = __bcast<1>(jp);
      const _v8sf jpz = __bcast<2>(jp);
      const _v8sf jpm = __bcast<3>(jp);

      const _v8sf dx = jpx - ipx;
      const _v8sf dy = jpy - ipy;
      const _v8sf dz = jpz - ipz;
      const _v8sf r2 = dx*dx + dy*dy + dz*dz + veps2;

      const _v8sf  rinv  = __builtin_ia32_rsqrtps256(r2);
      const _v8sf mrinv  = rinv  *  jpm;
      const _v8sf  rinv2 = rinv  *  rinv;
      const _v8sf mrinv3 = rinv2 * mrinv;

      fx   += mrinv3 * dx;
      fy   += mrinv3 * dy;
      fz   += mrinv3 * dz;
      gpot += mrinv;
    }
    gpot = -gpot;

    const _v8sf t0  = __builtin_ia32_unpcklps256(fx,   fz);
    const _v8sf t1  = __builtin_ia32_unpcklps256(fy, gpot);
    const _v8sf t2  = __builtin_ia32_unpckhps256(fx,   fz);
    const _v8sf t3  = __builtin_ia32_unpckhps256(fy, gpot);
    const _v8sf f04 = __builtin_ia32_unpcklps256(t0, t1);
    const _v8sf f15 = __builtin_ia32_unpckhps256(t0, t1);
    const _v8sf f26 = __builtin_ia32_unpcklps256(t2, t3);
    const _v8sf f37 = __builtin_ia32_unpckhps256(t2, t3);

    _v8sf* vforce = (_v8sf*)&force[i];
    *(vforce + 0) += __merge<0,0>(f04, f15);
    *(vforce + 1) += __merge<0,0>(f26, f37);
    *(vforce + 2) += __merge<1,1>(f04, f15);
    *(vforce + 3) += __merge<1,1>(f26, f37);
  }
#else
  const GroupT<Ng> &group = igroup;
  const int ni = group.nb();
  const int nj = nc;

  for (int i = 0; i < ni; i++)
  {
    const float4 ip = group[i].pos_mass();
    for (int j = 0; j < nj; j++)
    {
      const fMultipole &multipole = cell_list[j];
      const fMonopole   &m =   multipole.monopole();

      const float4 jp = *(float4*)&m;
      const float4 dr = jp - ip;
      const float  r2 = dr.norm2() + eps2;
      assert(r2 > 0.0f);
      const float  mj = jp.w();

      const float  rinv  = 1.0f/std::sqrt(r2);
      const float mrinv  = rinv*mj;
      const float mrinv3 = rinv*rinv*mrinv;

      float4 acc = dr * float4(mrinv3);
      acc.w() = -mrinv;

#if 0
      const Quadrupole &q = multipole.quadrupole();
      const int rinv2 = rinv*rinv;
      const int rinv4 = rinv2*rinv2;
      const int rinv5 = rinv4*rinv;
      const int rinv7 = rinv5*rinv2;

      const double dx = dr.x();
      const double dy = dr.y();
      const double dz = dr.z();

      const double Qrx = q.xx()*dx + q.xy()*dy + q.xz()*dz;
      const double Qry = q.xy()*dx + q.yy()*dy + q.yz()*dz;
      const double Qrz = q.xz()*dx + q.yz()*dy + q.zz()*dz;

      const double rQr = Qrx*dx + Qry*dy + Qrz*dz;
      float4 qacc1 = float4(2.5f*rinv7*rQr)*dr;
      qacc1.w() = (-0.5f)*rinv5*rQr;

      qacc1 = qacc1 - v4sf(rinv5)*(_v4sf){(float)Qrx, (float)Qry, (float)Qrz, 0.0f};
      acc = acc +  qacc1;
#endif

      force[i] = force[i] + acc;
    }
  }
#endif

  return 0;
}



#endif /* __GFORCE_GROUP_H__ */
