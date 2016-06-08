#include "Line.hh"

#include <cmath>

using namespace std;

namespace // anonymous
{
    double dot(vector<double> const &x,
               vector<double> const &y)
    {
        return x[0] * y[0] + x[1] * y[1];
    }
    
    double cross(vector<double> const &x,
                 vector<double> const &y)
    {
        return x[0] * y[1] - x[1] * y[0];
    }
}

Line::
Line(Surface_Type surface_type,
     vector<double> const &origin,
     vector<double> const &direction,
     vector<double> const &normal):
    Surface(2, // dimension
            surface_type),
    origin_(origin),
    direction_(direction),
    normal_(normal)
{
    
}

Line::Relation Line::
relation(vector<double> const &particle_position) const
{
    vector<double> distance(dimension_);
    
    for (int i = 0; i < dimension_; ++i)
    {
        distance[i] = particle_position[i] - origin_[i];
    }
    
    double k = dot(normal_, distance);

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
    double direction_cross = cross(direction_, particle_direction);
    
    if(direction_cross == 0)
    {
        return false;
    }

    vector<double> center_distance(dimension_);

    for (int i = 0; i < dimension_; ++i)
    {
        center_distance[i] = origin_[i] - particle_position[i];
    }
    
    double t = cross(direction_, center_distance) / direction_cross;

    if (t <= 0)
    {
        return false;
    }
    
    distance = t;
    
    position.resize(dimension_);
    
    for (int i = 0; i < dimension_; ++i)
    {
        position[i] = particle_position[i] + particle_direction[i] * t;
    }

    return true;
}

bool Line::
normal_direction(vector<double> const &position,
                 vector<double> &normal) const
{
    // Check whether point lies on line
    if (abs((position[1] - origin_[1]) * direction_[0]
            - (position[0] - origin_[0]) * direction_[1])
        > tolerance_)
    {
        return false;
    }
    
    normal = normal_;
    
    return true;
}
