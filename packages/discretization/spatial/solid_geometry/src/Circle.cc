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
    
    if (abs(r - radius_) <= relation_tolerance_)
    {
        return Relation::EQUAL;
    }
    else if (r < radius_)
    {
        return Relation::INSIDE;
    }
    else // if (r > radius_)
    {
        return Relation::OUTSIDE;
    }
}

double Circle::
distance(vector<double> const &position) const
{
    vector<double> const x = vf2::subtract(position, origin_);
    
    double const r = vf2::magnitude(x);
    
    return abs(r - radius_);
}

Circle::Intersection Circle::
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
    double const l2 = vf2::magnitude_squared(particle_direction);
    double const l3 = l1 * l1 - l0 * l2;

    if (l3 < 0)
    {
        return Intersection::NONE;
    }
    else if (l2 <= intersection_tolerance_)
    {
        return Intersection::PARALLEL;
    }
    
    double const l4 = sqrt(l3);
    
    double const s1 = (-l1 + l4) / l2;
    double const s2 = (-l1 - l4) / l2;
    
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
        return Intersection::NEGATIVE;
    }
    
    position = vf2::add(particle_position,
                        vf2::multiply(particle_direction,
                                      distance));
    
    if (l3 <= intersection_tolerance_)
    {
        return Intersection::TANGEANT;
    }
    else
    {
        return Intersection::INTERSECTS;
    }
}

bool Circle::
normal_direction(vector<double> const &position,
                 vector<double> &normal) const
{
    // Check if point lies on circle

    vector<double> const k0 = vf2::subtract(position,
                                            origin_);

    if (check_normal_)
    {
        if (abs(vf2::magnitude_squared(k0) - radius_ * radius_) > normal_tolerance_ * radius_ * radius_)
        {
            return false;
        }
    }
    
    normal = vf2::normalize(k0);
    
    return true;
}
