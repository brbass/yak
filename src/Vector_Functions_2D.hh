#ifndef Vector_Functions_2D_hh
#define Vector_Functions_2D_hh

#include "Check.hh"

namespace Vector_Functions_2D
{
    int const dim = 2;
    
    template<class T> vector<T> add(vector<T> const &x,
                                    vector<T> const &y)
    {
        Check(x.size() == dim);
        Check(y.size() == dim);

        return {x[0] + y[0], x[1] + y[1]};
    }
    template<class T> vector<T> subtract(vector<T> const &x,
                                         vector<T> const &y)
    {
        Check(x.size() == dim);
        Check(y.size() == dim);
        
        return {x[0] - y[0], x[1] - y[1]};
    }
    template<class T> vector<T> multiply(vector<T> const &x,
                                         T const t)
    {
        Check(x.size() == dim);

        return {x[0] * t, x[1] * t};
    }
    template<class T> T dot(vector<T> const &x,
                            vector<T> const &y)
    {
        Check(x.size() == dim);
        Check(y.size() == dim);

        return x[0] * y[0] + x[1] * y[1];
    }
    
    template<class T> T cross(vector<T> const &x,
                              vector<T> const &y)
    {
        Check(x.size() == dim);
        Check(y.size() == dim);

        return x[0] * y[1] - x[1] * y[0];
    }
    
    template<class T> T magnitude(vector<T> const &x)
    {
        Check(x.size() == dim);

        return sqrt(x[0] * x[0] + x[1] * x[1]);
    }
}

#endif
