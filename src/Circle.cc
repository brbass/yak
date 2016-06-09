#include "Circle.hh"

#include <cmath>

#include "Vector_Functions_2D.hh"

using namespace std;

namespace vf2 = Vector_Functions_2D;

/*
  Describes a circle with radius "radius" centered at the point "origin."
  The formula used is
  
  ||x - x0||^2 = r^2.
*/

Circle::
Circle(Surface_Type surface_type,
       double radius,
       vector<double> const &origin):
    Surface(2, // dimension
            surface_type),
    radius_(radius),
    origin_(origin)
{
}

Circle::Relation Circle::
relation(vector<double> const &particle_position) const
{
    vector<double> const x = vf2::subtract(particle_position, origin_);
    
    double const r = vf2::magnitude(x);
    
    if (r < radius_)
    {
        return Relation::INSIDE;
    }
    else if (r > radius_)
    {
        return Relation::OUTSIDE;
    }
    else // r == radius_
    {
        return Relation::EQUAL;
    }
}

bool Circle::
intersection(vector<double> const &particle_position,
             vector<double> const &particle_direction,
             double &distance,
             vector<double> &position) const
{
    vector<double> const k0 = vf2::subtract(particle_position,
                                            origin_);

    double const l0 = vf2::magnitude_squared(k0) - radius_ * radius_;
    double const l1 = vf2::dot(k0,
                               particle_direction);
    double const l2 = l1 * l1 - l0;

    if (l2 <= 0)
    {
        return false;
    }

    double const l3 = sqrt(l2);
    
    double const s1 = l1 + l3;
    double const s2 = l1 - l3;
    
    if (s2 > 0)
    {
        distance = s2;
    }
    else if (s1 > 0)
    {
        distance = s1;
    }
    else
    {
        return false;
    }
    
    position = vf2::add(particle_position,
                        vf2::multiply(particle_direction,
                                      distance));
    
    return true;
}

bool Circle::
normal_direction(vector<double> const &position,
                 vector<double> &normal) const
{
    // Check if point lies on circle

    vector<double> const k0 = vf2::subtract(position,
                                            origin_);
    
    if (abs(vf2::magnitude_squared(k0) - radius_ * radius_) > tolerance_)
    {
        return false;
    }
    
    normal = vf2::normalize(k0);
    
    return true;
}
