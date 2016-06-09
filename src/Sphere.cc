#include "Sphere.hh"

#include <cmath>

#include "Vector_Functions_3D.hh"

using namespace std;

namespace vf3 = Vector_Functions_3D;

Sphere::
Sphere(Surface_Type surface_type,
       double radius,
       vector<double> const &origin):
    Surface(3, // dimension
            surface_type),
    radius_(radius),
    origin_(origin)
{
}

Sphere::Relation Sphere::
relation(vector<double> const &particle_position) const
{
    vector<double> const k0 = vf3::subtract(particle_position,
                                            origin_);
    
    double const r = vf3::magnitude(k0);
    
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

Sphere::Intersection Sphere::
intersection(vector<double> const &particle_position,
             vector<double> const &particle_direction,
             double &distance,
             vector<double> &position) const
{
    vector<double> const k0 = vf3::subtract(particle_position,
                                            origin_);

    double const l0 = vf3::magnitude_squared(k0) - radius_ * radius_;
    double const l1 = vf3::dot(k0,
                               particle_direction);
    double const l2 = l1 * l1 - l0;

    if (l2 < 0)
    {
        return Intersection::NONE; // no intersection for line
    }
    
    double const l3 = sqrt(l2);
    
    double const s1 = -l1 + l3;
    double const s2 = -l1 - l3;

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
        return Intersection::NEGATIVE; // intersection behind current point
    }
    
    position = vf3::add(particle_position,
                        vf3::multiply(particle_direction,
                                      distance));
    
    if (l2 == 0)
    {
        return Intersection::TANGEANT;
    }
    else
    {
        return Intersection::INTERSECTS;
    }
}

bool Sphere::
normal_direction(vector<double> const &position,
                 vector<double> &normal) const
{
    // Check if point lies on sphere
    
    vector<double> const k0 = vf3::subtract(position,
                                            origin_);
    
    if (abs(vf3::magnitude_squared(k0) - radius_ * radius_) > tolerance_)
    {
        return false;
    }
    
    normal = vf3::normalize(k0);
    
    return true;
}
