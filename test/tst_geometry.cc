#include <memory>
#include <vector>

#include "Circle.hh"
#include "Line.hh"
#include "Region.hh"
#include "Solid_Geometry.hh"
#include "Surface.hh"

using namespace std;

void test_square()
{
    int dimension = 2;
    int num_surfaces = 4;
    int num_regions = 1;
    vector<double> origins = {-1, -1,
                              -1, 1,
                              1, 1,
                              1, -1};
    vector<double> directions = {0, 1,
                                 1, 0,
                                 0, -1,
                                 -1, 0};
    vector<double> normals = {-1, 0,
                              0, 1,
                              1, 0,
                              0, -1};
    
    vector<shared_ptr<Surface> > surfaces(num_surfaces);
    vector<shared_ptr<Region> > regions(num_regions);

    for (int i = 0; i < num_surfaces; ++i)
    {
        Surface::Surface_Type surface_type = Surface::Surface_Type::REFLECTIVE_BOUNDARY;
        vector<double> origin(dimension);
        vector<double> direction(dimension);
        vector<double> normal(dimension);

        for (int j = 0; j < dimension; ++j)
        {
            int k = j + dimension * i;

            origin[j] = origins[k];
            direction[j] = directions[k];
            normal[j] = normals[k];
        }
        
        surfaces[i] = make_shared<Line>(surface_type,
                                        origin,
                                        direction,
                                        normal);
    }
    
    vector<Surface::Relation> surface_relations(num_surfaces, Surface::Relation::NEGATIVE);

    regions[0] = make_shared<Region>();
    regions[0]->initialize(0, // material,
                           surface_relations,
                           surfaces);

    shared_ptr<Solid_Geometry> solid_geometry = make_shared<Solid_Geometry>(dimension,
                                                                            surfaces,
                                                                            regions);
}

void test_pincell()
{

}

int main()
{
    test_square();
}
