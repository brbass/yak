#ifndef Surface_hh
#define Surface_hh

#include <vector>

#include "Check.hh"

using std::vector;

/*
  General class for a solid geometry surface
*/
class Surface
{
public:

    /* Relationship of point to surface */
    enum class Relation
    {
        POSITIVE,
        NEGATIVE,
        EQUAL,
        OUTSIDE = POSITIVE,
        INSIDE = NEGATIVE
    };
    
    /* Type of surface */
    enum class Surface_Type
    {
        VACUUM_BOUNDARY,
        REFLECTIVE_BOUNDARY,
        INTERNAL
    };

    /* Types of intersection of vector and surface */
    enum class Intersection
    {
        INTERSECTS, // has intersection
        PARALLEL, // has no intersection, parallel to surface
        NONE, // has no intersection, negative or positive
        NEGATIVE, // only has negative intersection
        TANGEANT // only intersects at one infinitesimal point
    };

    /* Constructor */
    Surface(int dimension_,
            Surface_Type surface_type);
    
    /* Number of spatial dimensions */
    virtual int dimension()
    {
        return dimension_;
    }

    /* Returns surface type */
    virtual Surface_Type surface_type()
    {
        return surface_type_;
    }

    /* Returns relationship between point and surface */
    virtual Relation relation(vector<double> const &position) const = 0;
    
    /* Returns distance to the surface */
    virtual double distance(vector<double> const &position) const
    {
        AssertMsg(false, "Not yet implemented for this surface");
    }

    /* Type of intersection of streaming particle with surface
       If type is TANGEANT or PARALLEL, the distance and position are
       returned. Otherwise, the distance and position remain unchanged.*/
    virtual Intersection intersection(vector<double> const &initial_position,
                                      vector<double> const &initial_direction,
                                      double &distance,
                                      vector<double> &final_position) const = 0;

    /* Normal direction of surface at a point on the surface */
    virtual bool normal_direction(vector<double> const &position,
                                  vector<double> &normal) const = 0;

    /* Reflected direction */
    virtual bool reflected_direction(vector<double> const &position,
                                     vector<double> const &initial_direction,
                                     vector<double> &final_direction);
    
protected:

    bool check_normal_;
    int dimension_;
    Surface_Type surface_type_;
    double relation_tolerance_;
    double intersection_tolerance_;
    double normal_tolerance_;
};

#endif
