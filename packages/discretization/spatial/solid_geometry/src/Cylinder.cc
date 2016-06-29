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
    // Cross product approach
    // vector<double> const k0 = vf3::cross(direction_,
    //                                      vf3::subtract(particle_position,
    //                                                    origin_));

    // double const r = vf3::magnitude(k0);

    // Dot product approach
    vector<double> const k0 = vf3::subtract(particle_position,
                                            origin_);
    vector<double> const n = vf3::subtract(k0,
                                           vf3::multiply(direction_,
                                                         vf3::dot(direction_,
                                                                  k0)));
    double const r = vf3::magnitude(n);
    double const dr = r - radius_;

    if (abs(dr) <= relation_tolerance_)
    {
        return Relation::EQUAL;
    }
    else if (dr < 0)
    {
        return Relation::INSIDE;
    }
    else // if (dr > 0)
    {
        return Relation::OUTSIDE;
    }
}

Cylinder::Intersection Cylinder::
intersection(vector<double> const &particle_position,
             vector<double> const &particle_direction,
             double &distance,
             vector<double> &position) const
{
    // Cross product approach
    // vector<double> const k0 = vf3::cross(direction_,
    //                                      vf3::subtract(particle_position,
    //                                                    origin_));
    // vector<double> const k1 = vf3::cross(direction_,
    //                                      particle_direction);

    // Dot product approach
    vector<double> const j0 = vf3::subtract(particle_position,
                                            origin_);
    vector<double> const k0 = vf3::subtract(j0,
                                            vf3::multiply(direction_,
                                                          vf3::dot(direction_,
                                                                   j0)));
    vector<double> const k1 = vf3::subtract(particle_direction,
                                            vf3::multiply(direction_,
                                                          vf3::dot(direction_,
                                                                   particle_direction)));
    
    double const l0 = vf3::magnitude_squared(k0) - radius_ * radius_;
    double const l1 = vf3::dot(k0,
                               k1);
    double const l2 = vf3::magnitude_squared(k1);
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
    
    position = vf3::add(particle_position,
                        vf3::multiply(particle_direction,
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

bool Cylinder::
normal_direction(vector<double> const &position,
                 vector<double> &normal) const
{
    vector<double> const k0 = vf3::subtract(position,
                                            origin_);
    vector<double> const n = vf3::subtract(k0,
                                           vf3::multiply(direction_,
                                                         vf3::dot(direction_,
                                                                  k0)));

    if (check_normal_)
    {
        if (abs(vf3::magnitude_squared(n) - radius_ * radius_) > normal_tolerance_ * radius_ * radius_)
        {
            return false;
        }
    }
    
    normal = vf3::normalize(n);
    
    return true;
}
