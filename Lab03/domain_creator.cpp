#include <set>
#include <gmsh.h>
#include <iostream>

int main(int argc, char **argv)
{
  gmsh::initialize();

  gmsh::model::add("star");

  double lc = 0.02;

  gmsh::model::geo::addPoint(0, 4, 0, lc, 1);
  gmsh::model::geo::addPoint(0.9, 1.3, 0, lc, 2);
  gmsh::model::geo::addPoint(4, 1.3, 0, lc, 3);
  gmsh::model::geo::addPoint(1.5, -0.5, 0, lc, 4);
  gmsh::model::geo::addPoint(2.5, -3.4, 0, lc, 5);
  gmsh::model::geo::addPoint(0, -1.6, 0, lc, 6);
  gmsh::model::geo::addPoint(-2.5, -3.4, 0, lc, 7);
  gmsh::model::geo::addPoint(-1.5, -0.5, 0, lc, 8);
  gmsh::model::geo::addPoint(-4, 1.3, 0, lc, 9);
  gmsh::model::geo::addPoint(-0.9, 1.3, 0, lc, 10);
  for (unsigned i = 1; i < 10; ++i){
    gmsh::model::geo::addLine(i, i + 1, i);
  }
  gmsh::model::geo::addLine(10, 1, 10);

  gmsh::model::geo::addCurveLoop({1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, 1);
  gmsh::model::geo::addPlaneSurface({1}, 1);
  gmsh::model::geo::synchronize();

  gmsh::model::mesh::generate(2);
  gmsh::write("meshes/star.msh");

  gmsh::finalize();

  return 0;
}

