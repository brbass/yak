#include "Plane.hh"
#include "Vector_Functions_3D.hh"

#include <cmath>

using namespace std;

namespace vf3 = Vector_Functions_3D;

Plane::
Plane(Surface_Type surface_type,
      vector<double> const &origin,
      vector<double> const &normal):
    Surface(3, // dimension
            surface_type),
    origin_(origin),
    normal_(normal)
{
}

Plane::Relation Plane::
relation(vector<double> const &particle_position) const
{
    double const k = vf3::dot(normal_,
                              vf3::subtract(particle_position,
                                            origin_));
    
    if (k > 0)
    {
        return Relation::POSITIVE;
    }
    else if (k < 0)
    {
        return Relation::NEGATIVE;
    }
    else
    {
        return Relation::EQUAL;
    }
}

bool Plane::
intersection(vector<double> const &particle_position,
             vector<double> const &particle_direction,
             double &distance,
             vector<double> &position) const
{
    double direction_dot = vf3::dot(particle_direction,
                                    normal_);
    
    if (direction_dot == 0)
    {
        return false;
    }
    
    distance = vf3::dot(vf3::subtract(origin_,
                                      particle_position), normal_) / direction_dot;

    if (distance <= 0)
    {
        return false;
    }
    
    position = vf3::add(particle_position,
                        vf3::multiply(particle_direction,
                                      distance));
    
    return true;
}

bool Plane::
normal_direction(vector<double> const &position,
                 vector<double> &normal) const
{
    // Check whether point lies on plane
    if (vf3::dot(normal_,
                 vf3::subtract(position, origin_))
        > tolerance_)
    {
        return false;
    }
    
    normal = normal_;
    
    return true;
}
