#include "Cylinder.hh"
#include "Vector_Functions_3D.hh"

#include <cmath>

using namespace std;

namespace vf3 = Vector_Functions_3D;

/* 
   Describes an infinite cylinder with radius "radius" about the axis specified 
   by "direction" that goes through the point "origin." Two formulas are used,
   
   ||x - x0 - \Omega \cdot (x - x0) \Omega||^2 = r^2,
   ||\Omega \times (x - x0)||^2 = r^2,
   
   both of which should be valid with the magnitude applied. The first finds the
   magnitude of the vector pointing from the point x to the center of the cylinder,
   while the second finds the magnitude of the vector perpendicular to \Omega and 
   (x - x0) that runs from the origin to the surface of the cylinder. For the normal 
   direction, the first equation should be used. 
*/
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
    else if (l2 == 0)
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

    if (l3 == 0)
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

    if (abs(vf3::magnitude_squared(n) - radius_ * radius_) > tolerance_)
    {
        return false;
    }
    
    normal = vf3::normalize(n);
    
    return true;
}
