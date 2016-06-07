#include "Line.hh"

namespace // anonymous
{
    double cross(vector<double> const &x,
                 vector<double> const &y)
    {
        return x[0] * y[1] - x[1] * y[0];
    }
}

Line::
Line(vector<double> const &origin,
     vector<double> const &direction):
    origin_(origin),
    direction_(direction)
{
    
}

Relation Line::
relation(vector<double> &particle_position) const
{
    if (direction_[0] == 0) // line perpendicular to x axis
    {
        if (particle_position_[0] > origin_[0])
        {
            return Relation::POSITIVE;
        }
        else if (particle_position_[0] < origin_[0])
        {
            return Relation::NEGATIVE;
        }
        else
        {
            return Relation::EQUAL;
        }
    }
    else if (direction_[1] == 0) // line perpendicular to y axis
    {
        if (particle_position_[1] > origin_[1])
        {
            return Relation::POSITIVE;
        }
        else if (particle_position_[1] < origin_[1])
        {
            return Relation::NEGATIVE;
        }
        else
        {
            return Relation::EQUAL;
        }
    }
    
    vector<double> pos_direction = {1, 0};
    vector<double> neg_direction = {-1, 0};
    double distance;
    vector<double> position;
    
    // Check if higher on x axis
    if(intersection(particle_position,
                    pos_direction,
                    distance,
                    position))
    {
        return Relation::POSITIVE;
    }
    else if (distance == 0)
    {
        return Relation::EQUAL;
    }

    // Check if lower on x axis
    if (intersection(particle_position,
                     neg_direction,
                     distance,
                     position))
    {
        return Relation::NEGATIVE;
    }
    else if (distance == 0)
    {
        return Relation::EQUAL;
    }
    
    // Error
    AssertMsg(false, "relation not found");
    
    return Relation::ERROR;
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
