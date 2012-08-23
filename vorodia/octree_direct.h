public: template<const int N>
void buildDirectPolyhedron(Particle::Vector &ptcl, const Body &iBody, DirectPolyhedron<N> &d) const
{
  assert(isTreeReady());
  for (int k = 0; k < 8; k++)
    buildDirPol(ptcl, iBody, d, k);
}

private: template<const int N>
void buildDirPol(Particle::Vector &ptcl, const Body &iBody, DirectPolyhedron<N> &d, const int addr) const
{
  const Cell cell = cellList[addr];
  if (cell.isEmpty()) return;

  const boundary bnd = bndsList[cell.id()].inner();
  const vec3  c = bnd.center();
  const vec3  h = bnd.hlen();

  const vec3 dr = c - iBody.vector_pos();
  vec3 ds =  abseach(dr) - h;
  ds = (abseach(ds) + ds)*0.5*0.499;
  if (dr.x < 0.0) ds.x = -ds.x;
  if (dr.y < 0.0) ds.y = -ds.y;
  if (dr.z < 0.0) ds.z = -ds.z;
  assert(ds*dr >= 0.0);
  bool flag = ds.norm2() == 0.0;

  if (!flag)
  {

#define INTERSECT(x1) {\
  const dvec3 x = x1; \
  assert(x*dds> 0.0); \
  flag |= d.intersect(x * (x*dds/x.norm2()) ); }

    const dvec3 dds = ds;
    INTERSECT(dvec3(dr.x - h.x, dr.y - h.y, dr.z - h.z));
    INTERSECT(dvec3(dr.x + h.x, dr.y - h.y, dr.z - h.z));
    INTERSECT(dvec3(dr.x - h.x, dr.y + h.y, dr.z - h.z));
    INTERSECT(dvec3(dr.x + h.x, dr.y + h.y, dr.z - h.z));
    INTERSECT(dvec3(dr.x - h.x, dr.y - h.y, dr.z + h.z));
    INTERSECT(dvec3(dr.x + h.x, dr.y - h.y, dr.z + h.z));
    INTERSECT(dvec3(dr.x - h.x, dr.y + h.y, dr.z + h.z));
    INTERSECT(dvec3(dr.x + h.x, dr.y + h.y, dr.z + h.z));
    INTERSECT(dds);
    //  if (!flag) return;
  }

  /* check overlap */

  if (cell.isNode())
  {
    for (int k = 0; k < 8; k++)
      buildDirPol(ptcl, iBody, d, cell.addr() + k);
  }
  else
  {
    const Leaf &leaf = leafList[cell.leafIdx()];
    const int  ix    = iBody.idx();
    const dvec3 ipos  = ptcl[ix].pos;
    for (int j = 0; j < leaf.nb(); j++)
    {
      const int jx  = leaf[j ].idx();
      assert(overlapped(leaf[j].vector_pos(), bnd));
      if (ix != jx)
      {
        const dvec3 dr1 = ptcl[jx].pos - ipos;
        const bool flag1 = d.push(dr1*0.5, jx);
        const bool myflag = flag1 ? flag == flag1 : true;
        if (!myflag)
        {
          fprintf(stderr, " flag= %d  flag1= %d\n", flag, flag1);
          fprintf(stderr, " dr1= %g %g %g \n", dr1.x, dr1.y, dr1.z);
          fprintf(stderr, " dr= %g %g %g \n", dr.x, dr.y, dr1.z);
          fprintf(stderr, "  h= %g %g %g \n", h.x, h.y, h.z);
          fprintf(stderr, " ds= %g %g %g \n", ds.x, ds.y, ds.z);
          //        fprintf(stderr, "jpos= %g %g %g \n", ptcl[jx].pos.x, ptcl[jx].pos.y, ptcl[jx].pos.z);
          //         fprintf(stderr, "jpos= %g %g %g \n", leaf[j].vector_pos().x, leaf[j].vector_pos().y, leaf[j].vector_pos().z);
          //         fprintf(stderr, "   c= %g %g %g \n", c.x, c.y, c.z);
   
         const dvec3 dds = ds;
#define INTERSECT1(x1) {\
  const dvec3 x = x1; \
  assert(x*dds> 0.0); \
  const dvec3 y = x * (x*dds/x.norm2()); \
  const bool flag = d.intersect(y); \
  fprintf(stderr, " flag= %d  : x= %g %g %g  y= %g %g %g \n", flag, x.x,x.y,x.z, y.x,y.y,y.z); }
         INTERSECT1(dvec3(dr.x - h.x, dr.y - h.y, dr.z - h.z));
         INTERSECT1(dvec3(dr.x + h.x, dr.y - h.y, dr.z - h.z));
         INTERSECT1(dvec3(dr.x - h.x, dr.y + h.y, dr.z - h.z));
         INTERSECT1(dvec3(dr.x + h.x, dr.y + h.y, dr.z - h.z));
         INTERSECT1(dvec3(dr.x - h.x, dr.y - h.y, dr.z + h.z));
         INTERSECT1(dvec3(dr.x + h.x, dr.y - h.y, dr.z + h.z));
         INTERSECT1(dvec3(dr.x - h.x, dr.y + h.y, dr.z + h.z));
         INTERSECT1(dvec3(dr.x + h.x, dr.y + h.y, dr.z + h.z));
         INTERSECT1(dds);
         assert(0);
        }
      }
    }
  }
}

