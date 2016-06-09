#include "Line.hh"
#include "Vector_Functions_2D.hh"

#include <cmath>

using namespace std;

namespace vf2 = Vector_Functions_2D;

Line::
Line(Surface_Type surface_type,
     vector<double> const &origin,
     vector<double> const &normal):
    Surface(2, // dimension
            surface_type),
    origin_(origin),
    normal_(normal)
{
    
}

Line::Relation Line::
relation(vector<double> const &particle_position) const
{
    double const k = vf2::dot(normal_,
                              vf2::subtract(particle_position,
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

bool Line::
intersection(vector<double> const &particle_position,
             vector<double> const &particle_direction,
             double &distance,
             vector<double> &position) const
{
    double direction_dot = vf2::dot(particle_direction,
                                    normal_);

    if (direction_dot == 0)
    {
        return false;
    }
    
    distance = vf2::dot(vf2::subtract(origin_,
                                      particle_position), normal_) / direction_dot;

    if (distance <= 0)
    {
        return false;
    }
    
    // double const direction_cross = vf2::cross(direction_,
    //                                           particle_direction);
    
    // if(direction_cross == 0)
    // {
    //     return false;
    // }

    // vector<double> const center_distance = vf2::subtract(origin_,
    //                                                      particle_position);
    
    // double const t = vf2::cross(direction_,
    //                             center_distance) / direction_cross;
    
    // if (t <= 0)
    // {
    //     return false;
    // }
    
    // distance = t;
    
    position = vf2::add(particle_position,
                        vf2::multiply(particle_direction,
                                      distance));
    
    return true;
}

bool Line::
normal_direction(vector<double> const &position,
                 vector<double> &normal) const
{
    // Check whether point lies on line
    if (vf2::dot(normal_,
                 vf2::subtract(position, origin_))
        > tolerance_)
    {
        return false;
    }
    
    normal = normal_;
    
    return true;
}
