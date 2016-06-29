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
    vector<double> const k0 = vf2::subtract(particle_position,
                                            origin_);
    double const k1 = vf2::magnitude(k0);

    double const k = vf2::dot(normal_, 
                              k0);

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
    
    position = vf2::add(particle_position,
                        vf2::multiply(particle_direction,
                                      distance));
    
    return Intersection::INTERSECTS;
}

bool Line::
normal_direction(vector<double> const &position,
                 vector<double> &normal) const
{
    if (check_normal_)
    {
        vector<double> const k0 = vf2::subtract(position, origin_);
        double const k1 = vf2::magnitude(k0);

        // Check whether point lies on line
        if (vf2::dot(normal_, k0) > normal_tolerance_ * k1)
        {
            return false;
        }
    }
    
    normal = normal_;
    
    return true;
}
