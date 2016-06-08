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

    enum class Surface_Type
    {
        VACUUM_BOUNDARY,
        REFLECTIVE_BOUNDARY,
        INTERNAL
    };
    
    Surface(Surface_Type surface_type);

    virtual Surface_Type surface_type()
    {
        return surface_type_;
    }
    virtual Relation relation(vector<double> const &position) = 0;
    virtual bool intersection(vector<double> const &particle_position,
                              vector<double> const &particle_direction,
                              double &distance,
                              vector<double> &position) const = 0;
    
private:
    
    Surface_Type surface_type_;
}

#endif
