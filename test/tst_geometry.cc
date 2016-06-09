#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

#include "Circle.hh"
#include "Line.hh"
#include "Random_Number_Generator.hh"
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

    int num_tests = 10;

    Random_Number_Generator rng(-1.5, 1.5);

    int w = 12;
    
    cout << setw(w) << "x";
    cout << setw(w) << "y";
    cout << setw(w) << "Direction_x";
    cout << setw(w) << "Direction_y";
    cout << setw(w) << "Region";
    cout << setw(w) << "Surface";
    cout << setw(w) << "Distance";
    cout << setw(w) << "x_new";
    cout << setw(w) << "y_new";
    cout << endl;
    
    for (int i = 0; i < num_tests; ++i)
    {
        vector<double> position = rng.random_double_vector(dimension);
        
        vector<double> direction = rng.random_double_vector(dimension);
        double sum = 0;
        for (int i = 0; i < dimension; ++i)
        {
            sum += direction[i] * direction[i];
        }
        sum = sqrt(sum);
        for (int i = 0; i < dimension; ++i)
        {
            direction[i] /= sum;
        }

        double distance;
        vector<double> new_position;
        
        int region = solid_geometry->find_region(position);
        int surface = solid_geometry->next_intersection(position,
                                                        direction,
                                                        distance,
                                                        new_position);
        
        for (int i = 0; i < dimension; ++i)
        {
            cout << setw(w) << position[i];
        }
        for (int i = 0; i < dimension; ++i)
        {
            cout << setw(w) << direction[i];
        }
        cout << setw(w) << region;
        cout << setw(w) << surface;
        cout << setw(w) << distance;
        for (int i = 0; i < dimension; ++i)
        {
            cout << setw(w) << new_position[i];
        }
        cout << endl;
    }
}

void test_pincell()
{

}

int main()
{
    test_square();
}
