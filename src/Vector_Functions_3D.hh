#ifndef Vector_Functions_3D_hh
#define Vector_Functions_3D_hh

#include "Check.hh"

namespace Vector_Functions_3D
{
    int const dim = 3;
    
    template<class T> vector<T> add(vector<T> const &x,
                                    vector<T> const &y)
    {
        return {x[0] + y[0],
                x[1] + y[1],
                x[2] + y[2]};
    }
    template<class T> vector<T> subtract(vector<T> const &x,
                                         vector<T> const &y)
    {
        return {x[0] - y[0],
                x[1] - y[1],
                x[2] - y[2]};
    }
    template<class T> vector<T> multiply(vector<T> const &x,
                                         T const t)
    {
        return {x[0] * t,
                x[1] * t,
                x[2] * t};
    }
    template<class T> T dot(vector<T> const &x,
                            vector<T> const &y)
    {
        return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
    }
    
    template<class T> vector<T> cross(vector<T> const &x,
                                      vector<T> const &y)
    {
        return {x[1] * y[2] - x[2] * y[1],
                x[2] * y[0] - x[0] * y[2],
                x[0] * y[1] - x[1] * y[0]};
    }
    
    template<class T> T magnitude(vector<T> const &x)
    {
        return sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    }
}

#endif
