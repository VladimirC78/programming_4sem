#include <set>
#include <gmsh.h>

int main(int argc, char **argv)
{
    gmsh::initialize();
    gmsh::model::add("torus");
    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", 1.0);
    gmsh::option::setNumber("Mesh.CharacteristicLengthMin", 0.5);

    double R = 20, r = 10, h = 1;

    int pos_tor = gmsh::model::occ::addTorus(0, 0, 0, R, r);
    int neg_tor = gmsh::model::occ::addTorus(0, 0, 0, R, r - h);

    gmsh::vectorpair outDimTags;
    std::vector<gmsh::vectorpair> outDimTagsMap;

    gmsh::model::occ::cut({{3, pos_tor}}, {{3, neg_tor}}, outDimTags, outDimTagsMap, 3, true, false);

    gmsh::model::occ::synchronize();
    gmsh::model::mesh::generate(3);

    gmsh::write("torus.msh");
    gmsh::finalize();

    return 0;
}