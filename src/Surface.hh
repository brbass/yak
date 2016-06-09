#ifndef Surface_hh
#define Surface_hh

#include <vector>

using std::vector;

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

    enum class Intersection
    {
        INTERSECTS, // has intersection
        PARALLEL, // has no intersection, parallel to surface
        NONE, // has no intersection, negative or positive
        NEGATIVE, // only has negative intersection
        TANGEANT // only intersects at one infinitesimal point
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
    virtual Relation relation(vector<double> const &particle_position) const = 0;
    virtual Intersection intersection(vector<double> const &particle_position,
                                      vector<double> const &particle_direction,
                                      double &distance,
                                      vector<double> &position) const = 0;
    virtual bool normal_direction(vector<double> const &position,
                                  vector<double> &normal) const = 0;
    
protected:

    int dimension_;
    double tolerance_;
    Surface_Type surface_type_;
};

#endif
