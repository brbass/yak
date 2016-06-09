#include "Circle.hh"

#include <cmath>

#include "Vector_Functions_2D.hh"

using namespace std;

namespace vf2 = Vector_Functions_2D;

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
    vector<double> const center_distance = vf2::subtract(origin_, particle_position);
    
    double const cross_val = vf2::cross(particle_direction, center_distance);
    double const k1 = radius_ * radius_ - cross_val * cross_val;
    
    if (k1 <= 0)
    {
        return false;
    }

    double const k2 = sqrt(k1);
    double const k3 = vf2::dot(center_distance, particle_direction);
    
    double const t1 = k3 + k2;
    double const t2 = k3 - k2;
    
    if (t2 > 0)
    {
        distance = t2;
    }
    else if (t1 > 0)
    {
        distance = t1;
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
    if (abs(pow(position[0] - origin_[0], 2) + pow(position[1] - origin_[1], 2)
            - radius_ * radius_) > tolerance_)
    {
        return false;
    }
    
    normal = vf2::subtract(position,
                           origin_);
    
    return true;
}
