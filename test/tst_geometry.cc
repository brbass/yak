#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

#include "Circle.hh"
#include "Cylinder.hh"
#include "Line.hh"
#include "Plane.hh"
#include "Random_Number_Generator.hh"
#include "Region.hh"
#include "Solid_Geometry.hh"
#include "Sphere.hh"
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

void cube_with_sphere_minus_cylinder()
{
    int dimension = 3;
    int num_surfaces = 8;
    int num_regions = 2;
    
    vector<shared_ptr<Surface> > surfaces(num_surfaces);
    vector<shared_ptr<Region> > regions(num_regions);

    vector<double> x0 = {-1, 0, 0};
    vector<double> x1 = {1, 0, 0};
    vector<double> x2 = {0, -1, 0};
    vector<double> x3 = {0, 1, 0};
    vector<double> x4 = {0, 0, -1};
    vector<double> x5 = {0, 0, 1};
    vector<double> x6 = {0, 0, 0};
    vector<double> x7 = {0, 0, 0};

    vector<double> n0 = {-1, 0, 0};
    vector<double> n1 = {1, 0, 0};
    vector<double> n2 = {0, -1, 0};
    vector<double> n3 = {0, 1, 0};
    vector<double> n4 = {0, 0, -1};
    vector<double> n5 = {0, 0, 1};
    vector<double> n6 = {0, 1, 0};

    double r6 = 0.25;
    double r7 = 0.5;
    
    surfaces[0] = make_shared<Plane>(Surface::Surface_Type::REFLECTIVE_BOUNDARY,
                                     x0,
                                     n0);
    surfaces[1] = make_shared<Plane>(Surface::Surface_Type::REFLECTIVE_BOUNDARY,
                                     x1,
                                     n1);
    surfaces[2] = make_shared<Plane>(Surface::Surface_Type::REFLECTIVE_BOUNDARY,
                                     x2,
                                     n2);
    surfaces[3] = make_shared<Plane>(Surface::Surface_Type::REFLECTIVE_BOUNDARY,
                                     x3,
                                     n3);
    surfaces[4] = make_shared<Plane>(Surface::Surface_Type::REFLECTIVE_BOUNDARY,
                                     x4,
                                     n4);
    surfaces[5] = make_shared<Plane>(Surface::Surface_Type::REFLECTIVE_BOUNDARY,
                                     x5,
                                     n5);
    surfaces[6] = make_shared<Cylinder>(Surface::Surface_Type::INTERNAL,
                                        r6,
                                        x6,
                                        n6);
    surfaces[7] = make_shared<Sphere>(Surface::Surface_Type::INTERNAL,
                                        r7,
                                        x7);

    int num_surfaces0 = 2;
    int num_surfaces1 = 6;
    int num_regions1 = 1;

    regions[0] = make_shared<Region>();
    regions[1] = make_shared<Region>();

    vector<shared_ptr<Surface> > surfaces0 = {surfaces[6], surfaces[7]};
    vector<shared_ptr<Surface> > surfaces1 = {surfaces[0], surfaces[1], surfaces[2], surfaces[3], surfaces[4], surfaces[5]};
    vector<shared_ptr<Region> > regions1 = {regions[0]};
    
    vector<Surface::Relation> surface_relations0 = {Surface::Relation::OUTSIDE,
                                                    Surface::Relation::INSIDE};
    vector<Surface::Relation> surface_relations1(num_surfaces1, Surface::Relation::NEGATIVE);
    vector<Region::Relation> region_relations1 = {Region::Relation::OUTSIDE};
    
    regions[0]->initialize(0, // material,
                           surface_relations0,
                           surfaces0);
    regions[1]->initialize(1,
                           surface_relations1,
                           surfaces1,
                           region_relations1,
                           regions1);
    
    shared_ptr<Solid_Geometry> solid_geometry = make_shared<Solid_Geometry>(dimension,
                                                                            surfaces,
                                                                            regions);
    
    int num_tests = 30;

    Random_Number_Generator rng(-1.5, 1.5);

    int w = 12;
    
    cout << setw(w) << "x";
    cout << setw(w) << "y";
    cout << setw(w) << "z";
    cout << setw(w) << "Direction_x";
    cout << setw(w) << "Direction_y";
    cout << setw(w) << "Direction_z";
    cout << setw(w) << "Region";
    cout << setw(w) << "Surface";
    cout << setw(w) << "Distance";
    cout << setw(w) << "x_new";
    cout << setw(w) << "y_new";
    cout << setw(w) << "z_new";
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

int main()
{
    cube_with_sphere_minus_cylinder();
}
