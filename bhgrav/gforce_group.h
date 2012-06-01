#ifndef __GFORCE_GROUP_H__
#define __GFORCE_GROUP_H__

inline bool split_node(
    const float4 cellCoM,
    const float4 groupCentre,
    const float4 groupSize)

{
  const float4  dr = fabs(groupCentre - cellCoM) - groupSize;
  const float4 adr = dr + fabs(dr);  /* 0 if dr < 0, otherwise 2*dr */
  const float  ds2 = adr.norm2()*0.25f;

  return ds2 <= cellCoM.w();
}

template<const bool ROOT, const int Ng, const int Nc, const int Np>
void gforce(
    float4     force[Ng  ],         /* maximal number of particles in a group  */
    int    cell_list[Nc*2], int &nc, /* list for particle-cell interactions     */
    float4 ptcl_list[Np*2], int &np, /* list for particle-particle interactions */
    const  GroupT<Ng> &group, 
    const  float4 &groupCentre = 0.0, const float4 &groupSize = 0.0, const int addr = 0) const
{
  if (ROOT)
  {
    assert(isTreeReady());
    assert(Np >= 8*NLEAF);
    assert(Nc >  8);
    const boundary bnd = group.outerBoundary();
    const float4 groupCentre = bnd.center();
    const float4 groupSize   = bnd.hlen()*(real)2.0;
    nc = np = 0;
    for (int i = 0; i < Ng; i++)
      force[i] = 0.0f;
    for (int k = 0; k < 8; k++)
      if (!cellList[k].isEmpty())
        if (split_node(cellCoM[cellList[k].id()], groupCentre, groupSize))
          gforce<false>(force, cell_list, nc, ptcl_list, np,  group, groupCentre, groupSize, k);
  }
  else
  {
    const Cell cell = cellList[addr];
    if (cell.isNode())
    {
      for (int k = 0; k < 8; k++)
        if (!cellList[cell.addr()+k].isEmpty())
        {
          if (split_node(cellCoM[cellList[cell.addr()+k].id()], groupCentre, groupSize))
            gforce<false>(force, cell_list, nc, ptcl_list, np, group, groupCentre, groupSize, cell.addr()+k);
        }
        else
          cell_list[nc++] = cell.addr() + k;
    }
    else
    {
      const Leaf  leaf = leafList[cell.leafIdx()];
      for (int i = 0; i < leaf.nb(); i++)
        ptcl_list[np++] = leaf[i].pos_mass();
    }

    if (nc >= Nc)
    {
      nc -= Nc;
    }
    if (np >= Np)
    {
      np -= Np;
    }
  }
}

#endif /* __GFORCE_GROUP_H__ */
