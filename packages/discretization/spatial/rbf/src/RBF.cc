#include "RBF.hh"

#include <iostream>
#include <string>

#include "Check.hh"

using namespace std;

RBF::
RBF(int number_of_dimensions,
    vector<double> const &position,
    vector<double> const &shape_parameter):
    number_of_dimensions_(number_of_dimensions),
    position_(position),
    shape_(shape_parameter)
{
    check_class_invariants();
}

void RBF::
check_class_invariants()
{
    Assert(shape_.size() == number_of_dimensions());
}

// double RBF::
// basis(vector<double> const &r) const
// {
//     Check(r.size() == number_of_dimensions());

//     double const dist2 = get_distance_squared(r);
//     double const dist = sqrt(dist2);

//     switch(basis_type())
//     {
//     case GAUSSIAN:
//         return exp(-dist2);
//     case MULTIQUADRIC:
//         return sqrt(1 + dist2);
//     case INVERSE_MULTIQUADRIC:
//         return 1/sqrt(1 + dist2);
//     case WENDLAND30:
//         if(dist<1)
//         {
//             return pow(1-dist, 2);
//         }
//         return 0;
//     case WENDLAND31:
//         if (dist<1)
//         {
//             return pow(1-dist, 4)*(4*dist+1);
//         }
//         return 0;
//     case WENDLAND32:
//         if (dist<1)
//         {
//             return pow(1-dist, 6)*(35*dist*dist+18*dist+3);
//         }
//         return 0;
//     case WENDLAND33:
//         if (dist<1)
//         {
//             return pow(1-dist, 8)*(32*dist*dist*dist+18*dist*dist+8*dist+1);
//         }
//         return 0;
//     default:
//         AssertMsg(false, "no such type of basis function");
//         return 0;
//     }
// }

// double RBF::
// dbasis(int const dim,
//        vector<double> const &r) const
// {
//     Check(r.size() == number_of_dimensions());

//     double const dist2 = get_distance_squared(r);
//     double const dist = sqrt(dist2);
//     double const kr2 = shape2(dim)*r[dim];
    
//     switch(basis_type())
//     {
//     case GAUSSIAN:
//         return -2*kr2*exp(-dist2);
//     case MULTIQUADRIC:
//         return kr2/sqrt(1 + dist2);
//     case INVERSE_MULTIQUADRIC:
//         return -kr2*pow(1 + dist2, -1.5);
//     case WENDLAND30:
//         if (dist==0.)
//         {
//             return -2 * shape(dim);
//         }
//         else if(dist<1)
//         {
//             return kr2*2*(1-1/dist);
//         }
//         return 0;
//     case WENDLAND31:
//         if (dist<1)
//         {
//             return kr2*20*pow(dist-1, 3);
//         }
//         return 0;
//     case WENDLAND32:
//         if (dist<1)
//         {
//             return kr2*56*pow(dist-1, 5)*(5*dist+1);
//         }
//         return 0;
//     case WENDLAND33:
//         if (dist<1)
//         {
//             return kr2*22*pow(dist-1, 7)*(1+dist*(16*dist+7));
//         }
//         return 0;
//     default:
//         AssertMsg(false, "no such type of basis function");
//         return 1;
//     }
// }

// double RBF::
// ddbasis(int const dim,
//         vector<double> const &r) const
// {
//     Check(r.size() == number_of_dimensions());
    
//     double const dist2 = get_distance_squared(r);
//     double const dist = sqrt(dist2);
//     double const k2r2 = shape2(dim)*r[dim]*r[dim];

//     switch(basis_type())
//     {
//     case GAUSSIAN:
//         return 2*(2*k2r2-1)*shape2(dim)*exp(-dist2);
//     case MULTIQUADRIC:
//         return shape2(dim)*(1+dist2-k2r2)*pow(1+dist2, -1.5);
//     case INVERSE_MULTIQUADRIC:
//         return -shape2(dim)*(1+dist2-3*k2r2)*pow(1 + dist2, -2.5);
//     case WENDLAND30:
//         if (dist==0.)
//         {
//             return -2 * shape2(dim);
//         }
//         else if(dist<1)
//         {
//             return 2 * shape2(dim) * (dist2 * (dist - 1) + k2r2) / (dist2 * dist);
//         }
//         return 0;
//     case WENDLAND31:
//         if (dist==0.)
//         {
//             return -20 * shape2(dim);
//         }
//         else if (dist<1)
//         {
//             return 20 * shape2(dim) * pow(dist - 1, 2)* (dist*(dist-1) + 3 * k2r2) / dist;
//         }
//         return 0;
//     case WENDLAND32:
//         if (dist<1)
//         {
//             return 56./3.*shape2(dim) * pow(dist-1, 4)*(-1 + dist*(-4+5*dist)+30*k2r2);
//         }
//         return 0;
//     case WENDLAND33:
//         if (dist<1)
//         {
//             return 22 * shape2(dim)* pow(dist - 1, 6)*((dist-1)*(1+dist*(7+16*dist)) + 24*k2r2*(1+6*dist));
//         }
//         return 0;
//     default:
//         AssertMsg(false, "no such type of basis function");
//         return 1;
//     }
// }

double RBF::
get_distance_squared(vector<double> const &r) const
{
    Check(r.size() == number_of_dimensions_);

    double sum = 0;
    
    for (int i = 0; i < number_of_dimensions(); ++i)
    {
        double dist = shape_[i] * (r[i] - position_[i]);
        
        sum += dist * dist;
    }

    return sum;
}

void RBF::
multiply_shape_parameter(double t)
{
    for (int d = 0; d < number_of_dimensions_; ++d)
    {
        shape_[d] = shape_[d] * t;
    }
}
