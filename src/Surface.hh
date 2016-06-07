#ifndef Surface_hh
#define Surface_hh

class Surface
{
public:
    
    enum class Relation
    {
        POSITIVE,
        NEGATIVE,
        EQUAL,
        INSIDE = NEGATIVE,
        OUTSIDE = POSITIVE,
        ERROR
    };

    Surface();

private:
    
    virtual Relation relation(vector<double> const &position) = 0;
    
    virtual bool intersection(vector<double> const &particle_position,
                              vector<double> const &particle_direction,
                              double &distance,
                              vector<double> &position) const = 0;
    
}

#endif
