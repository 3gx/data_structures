
namespace Voronoi
{

  struct Diagram
  {

    void cell(Polyhedron &c)
    {
      if (!poly.intersect(boundary)) return;

      if (is_leaf())
        for (int i = 0; i < np; i++)
          c.push(ptcl[i].pos, ptcl[i].id);
      else
        for (int ch = 0 ; ch < 8; ch++)
          leaf[ch].cell(c);
    }
  }
}
