#include "Surface.hh"

Surface::
Surface(int dimension,
        Surface_Type surface_type,
        double tolerance):
    dimension_(dimension),
    surface_type_(surface_type),
    tolerance_(tolerance)
{
};
