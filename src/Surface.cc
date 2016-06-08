#include "Surface.hh"

Surface::
Surface(int dimension,
        Surface_Type surface_type):
    dimension_(dimension),
    tolerance_(1e-12),
    surface_type_(surface_type)
{
};
