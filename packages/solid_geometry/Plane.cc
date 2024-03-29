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

    if (abs(k) <= relation_tolerance_)
    {
        return Relation::EQUAL;
    }
    else if (k > 0)
    {
        return Relation::POSITIVE;
    }
    else // if (k < 0)
    {
        return Relation::NEGATIVE;
    }
}

Plane::Intersection Plane::
intersection(vector<double> const &particle_position,
             vector<double> const &particle_direction,
             double &distance,
             vector<double> &position) const
{
    vector<double> const k0 = vf3::subtract(origin_,
                                            particle_position);
    double const l0 = vf3::dot(k0,
                               normal_);
    double const l1 = vf3::dot(particle_direction,
                               normal_);
    
    if (abs(l1) <= intersection_tolerance_)
    {
        return Intersection::PARALLEL;
    }

    double const s = l0 / l1;

    if (s > 0)
    {
        distance = s;
    }
    else
    {
        return Intersection::NEGATIVE;
    }
    
    position = vf3::add(particle_position,
                        vf3::multiply(particle_direction,
                                      distance));
    
    return Intersection::INTERSECTS;
}

bool Plane::
normal_direction(vector<double> const &position,
                 vector<double> &normal) const
{
    // Check whether point lies on plane

    if (check_normal_)
    {
        vector<double> const k0 = vf3::subtract(position, origin_);
        double k1 = vf3::magnitude(k0);

        if (vf3::dot(normal_, k0) > normal_tolerance_ * k1)
        {
            return false;
        }
    }
    
    normal = normal_;
    
    return true;
}
