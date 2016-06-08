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
        OUTSIDE = POSITIVE,
        INSIDE = NEGATIVE
    };

    enum class Surface_Type
    {
        VACUUM_BOUNDARY,
        REFLECTIVE_BOUNDARY,
        INTERNAL
    };

    Surface(int dimension_,
            Surface_Type surface_type);

    virtual int dimension()
    {
        return dimension_;
    }
    virtual Surface_Type surface_type()
    {
        return surface_type_;
    }
    virtual Relation relation(vector<double> const &position) = 0;
    virtual bool intersection(vector<double> const &particle_position,
                              vector<double> const &particle_direction,
                              double &distance,
                              vector<double> &position) const = 0;
    virtual bool normal_direction(vector<double> const &position,
                                  vector<double> &normal) const = 0;
    
private:

    int dimension_;
    double tolerance_;
    Surface_Type surface_type_;
}

#endif
