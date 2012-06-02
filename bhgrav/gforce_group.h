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
  float4 ptcl_list[Nptcl*2];
  int    cell_list[Ncell*2];
  int nc = 0, np = 0;
  int np_tot = 0, nc_tot = 0;

  gForce<true, Ncell, Nptcl>(group, force, cell_list, nc, ptcl_list, np, np_tot, nc_tot);

  return std::make_pair(np_tot*group.nb(), nc_tot*group.nb());
}

template<const bool ROOT, const int Nc, const int Np, const int Ng>
void gForce(
    const GroupT<Ng> &group,
    float4     force[Ng  ],         /* maximal number of particles in a group  */
    int    cell_list[Nc*2], int &nc, /* list for particle-cell interactions     */
    float4 ptcl_list[Np*2], int &np, /* list for particle-particle interactions */
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
      cell_list[nc++] = cell.id();
      nc_tot++;
      assert(nc <= 2*Nc);
    }

    if (np >= Np) np = particle_particle<Np>(group, force, ptcl_list, np);
    if (nc >= Nc) nc = particle_cell    <Nc>(group, force, cell_list, nc);
  }
}

template<const int Np, const int Ng>
int particle_particle(
    const  GroupT<Ng> &group,
    float4     force[Ng  ],       
    float4 ptcl_list[Np*2], int np) const
{
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

  return 0;
}

template<const int Nc, const int Ng>
int particle_cell(
    const GroupT<Ng> &group,
    float4    force[Ng  ],       
    int   cell_list[Nc*2], int nc) const
{
  const int ni = group.nb();
  const int nj = nc;

  for (int i = 0; i < ni; i++)
  {
    const float4 ip = group[i].pos_mass();
    for (int j = 0; j < nj; j++)
    {
      const fMultipole &multipole = multipoleList[cell_list[j]];
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

  return 0;
}



#endif /* __GFORCE_GROUP_H__ */
