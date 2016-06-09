#include "Cylinder.hh"
#include "Vector_Functions_3D.hh"

#include <cmath>

using namespace std;

namespace vf3 = Vector_Functions_3D;

Cylinder::
Cylinder(Surface_Type surface_type,
         double radius,
         vector<double> const &origin,
         vector<double> const &direction):
    Surface(3, // dimension
            surface_type),
    radius_(radius),
    origin_(origin),
    direction_(direction)
{
}

Cylinder::Relation Cylinder::
relation(vector<double> const &particle_position) const
{
    return Relation::EQUAL;
}

bool Cylinder::
intersection(vector<double> const &particle_position,
             vector<double> const &particle_direction,
             double &distance,
             vector<double> &position) const
{
    return true;
}

bool Cylinder::
normal_direction(vector<double> const &position,
                 vector<double> &normal) const
{
    return true;
}
