#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

#include "Discrete_To_Moment.hh"
#include "Energy_Discretization.hh"
#include "Gauss_Legendre_Quadrature.hh"
#include "LDFE_Quadrature.hh"
#include "Math_Functions.hh"
#include "Moment_To_Discrete.hh"
#include "Quadrature_Rule.hh"
#include "Random_Number_Generator.hh"
#include "Simple_Spatial_Discretization.hh"

using namespace std;

void print(vector<double> &x0, vector<double> &x)
{
    for (unsigned i = 0; i < x.size(); ++i)
    {
        cout << setw(16) << x0[i];
        cout << setw(16) << x[i];
        cout << endl;
    }
    
    cout << endl;
}

void test_sph_harm()
{
    using namespace Math_Functions;
    
    Random_Number_Generator rng(-1, 1);

    double x = rng.random_double();
    double y = rng.random_double();
    double z = rng.random_double();

    x = 1;
    y = 1;
    z = 0.0;
    
    double sum = sqrt(x * x + y * y + z * z);

    x /= sum;
    y /= sum;
    z /= sum;

    for (int l = 0; l < 4; ++l)
    {
        for (int m = -l; m <= l; ++m)
        {
            cout << setw(16) << l;
            cout << setw(16) << m;
            cout << setw(16) << spherical_harmonic(l, m, x, y, z);
            cout << endl;
        }
    }
}


void test_gauss_legendre(int number_of_moments,
                         int number_of_ordinates)
{
    int dimension = 1;
    int number_of_cells = 2;
    int number_of_nodes = 2;
    int number_of_groups = 2;
    int number_of_boundary_cells = 2;

    Random_Number_Generator rng(0, 1);
    vector<double> x0 = rng.random_double_vector(number_of_cells * number_of_nodes * number_of_groups * number_of_moments);
    vector<double> y0 = rng.random_double_vector(number_of_cells * number_of_nodes * number_of_groups * number_of_ordinates);

    shared_ptr<Spatial_Discretization> spatial_discretization 
        = make_shared<Simple_Spatial_Discretization>(dimension,
                                                     number_of_cells,
                                                     number_of_nodes,
                                                     number_of_boundary_cells,
                                                     Spatial_Discretization::SLAB);
    shared_ptr<Angular_Discretization> angular_discretization 
        = make_shared<Gauss_Legendre_Quadrature>(dimension,
                                                 number_of_moments,
                                                 number_of_ordinates);
    shared_ptr<Energy_Discretization> energy_discretization
        = make_shared<Energy_Discretization>(number_of_groups);
    
    shared_ptr<Vector_Operator> upM 
        = make_shared<Moment_To_Discrete>(spatial_discretization,
                                          angular_discretization,
                                          energy_discretization);
    
    shared_ptr<Vector_Operator> upD 
        = make_shared<Discrete_To_Moment>(spatial_discretization,
                                          angular_discretization,
                                          energy_discretization);
    
    vector<double> x(x0);

    (*upM)(x);
    (*upD)(x);
    print(x0, x);
}

void test_ldfe(int number_of_scattering_moments,
               int rule,
               int dimension)
{
    int number_of_cells = 1;
    int number_of_nodes = 1;
    int number_of_groups = 1;
    int number_of_boundary_cells = 1;
    
    Spatial_Discretization::Geometry geometry;
    
    switch(dimension)
    {
    case 1:
        geometry = Spatial_Discretization::SLAB;
        break;
    case 2:
        geometry = Spatial_Discretization::RECTANGLE;
        break;
    case 3:
        geometry = Spatial_Discretization::CUBOID;
        break;
    }
    
    shared_ptr<Spatial_Discretization> spatial_discretization
        = make_shared<Simple_Spatial_Discretization>(dimension,
                                                     number_of_cells,
                                                     number_of_nodes,
                                                     number_of_boundary_cells,
                                                     geometry);

    shared_ptr<LDFE_Quadrature> angular_discretization 
        = make_shared<LDFE_Quadrature>(dimension,
                                       number_of_scattering_moments,
                                       rule);
    shared_ptr<Energy_Discretization> energy_discretization
        = make_shared<Energy_Discretization>(number_of_groups);
    
    shared_ptr<Vector_Operator> upM 
        = make_shared<Moment_To_Discrete>(spatial_discretization,
                                          angular_discretization,
                                          energy_discretization);
    
    shared_ptr<Vector_Operator> upD 
        = make_shared<Discrete_To_Moment>(spatial_discretization,
                                          angular_discretization,
                                          energy_discretization);

    int number_of_ordinates = angular_discretization->number_of_ordinates();
    int number_of_moments = angular_discretization->number_of_moments();
    
    vector<double> const ordinates = angular_discretization->ordinates();
    
    // for (int i = 0; i < number_of_ordinates; ++i)
    // {
    //     cout << setw(16) << i;
    //     for (int d = 0; d < dimension; ++d)
    //     {
    //         cout << setw(16) << ordinates[d + dimension * i];
    //     }
    //     cout << endl;
    // }

    Random_Number_Generator rng(0, 1);
    vector<double> x0 = rng.random_double_vector(number_of_cells * number_of_nodes * number_of_groups * number_of_moments);
    vector<double> y0 = rng.random_double_vector(number_of_cells * number_of_nodes * number_of_groups * number_of_ordinates);

    vector<double> x(x0);
    vector<double> y(y0);
    
    (*upM)(x);
    (*upD)(x);
    print(x0, x);
}

int main(int argc, char **argv)
{
    if (argc != 4)
    {
        cerr << "usage: tst_moment_discrete [num_moments num_ordinates dimension]" << endl;
        return 1;
    }
    
    int number_of_moments = atoi(argv[1]);
    int number_of_ordinates = atoi(argv[2]);
    int dimension = atoi(argv[3]);

    switch(dimension)
    {
    case 1:
        test_gauss_legendre(number_of_moments,
                            number_of_ordinates);
        break;
    case 2:
        test_ldfe(number_of_moments,
                  number_of_ordinates,
                  dimension);
        break;
    case 3:
        test_ldfe(number_of_moments,
                  number_of_ordinates,
                  dimension);
        break;
    default:
        cerr << "dimension must be 1, 2, or 3" << endl;
        return 1;
        break;
    }

    

    
}
