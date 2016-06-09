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

Line::Intersection Line::
intersection(vector<double> const &particle_position,
             vector<double> const &particle_direction,
             double &distance,
             vector<double> &position) const
{
    vector<double> const k0 = vf2::subtract(origin_,
                                            particle_position);
    double const l0 = vf2::dot(k0,
                               normal_);
    double const l1 = vf2::dot(particle_direction,
                               normal_);
    
    if (l1 == 0)
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
    
    position = vf2::add(particle_position,
                        vf2::multiply(particle_direction,
                                      distance));
    
    return Intersection::INTERSECTS;
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
