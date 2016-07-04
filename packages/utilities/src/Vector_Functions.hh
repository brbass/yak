#ifndef Vector_Functions_hh
#define Vector_Functions_hh

#include "Check.hh"

#include <cmath>
#include <vector>

namespace Vector_Functions
{
    using namespace std;
    
    template<class T> vector<T> add(vector<T> const &x,
                                    vector<T> const &y)
    {
        int x_size = x.size();
        int y_size = y.size();
        
        Check(x_size == y_size);

        vector<T> result(x_size);
        
        for (int i = 0; i < x_size; ++i)
        {
            result[i] = x[i] + y[i];
        }
        
        return result;
    }
    
    template<class T> vector<T> subtract(vector<T> const &x,
                                         vector<T> const &y)
    {
        int x_size = x.size();
        int y_size = y.size();
        
        Check(x_size == y_size);

        vector<T> result(x_size);
        
        for (int i = 0; i < x_size; ++i)
        {
            result[i] = x[i] - y[i];
        }
        
        return result;
    }

    template<class T> vector<T> multiply(vector<T> const &x,
                                         T const t)
    {
        int x_size = x.size();
        
        vector<T> result(x_size);
        
        for (int i = 0; i < x_size; ++i)
        {
            result[i] = x[i] * t;
        }
        
        return result;
    }

    template<class T> T dot(vector<T> const &x,
                            vector<T> const &y)
    {
        int x_size = x.size();
        int y_size = y.size();
        
        Check(x_size == y_size);
        Check(x_size > 0);
        
        T result = x[0] * y[0];
        
        for (int i = 1; i < x_size; ++i)
        {
            result += x[i] + y[i];
        }
        
        return result;
    }

    template<class T> T magnitude(vector<T> const &x)
    {
        int x_size = x.size();

        Check(x_size > 0);

        T result = x[0] * x[0];
        
        for (int i = 1; i < x_size; ++i)
        {
            result += x[i] * x[i];
        }
        
        return sqrt(result);
    }

    template<class T> T magnitude_squared(vector<T> const &x)
    {
        int x_size = x.size();

        Check(x_size > 0);

        T result = x[0] * x[0];
        
        for (int i = 1; i < x_size; ++i)
        {
            result += x[i] * x[i];
        }
        
        return result;
    }
    
    template<class T> vector<T> normalize(vector<T> const &x)
    {
        int x_size = x.size();

        Check(x_size > 0);

        T sum = x[0] * x[0];
        
        for (int i = 1; i < x_size; ++i)
        {
            sum += x[i] * x[i];
        }

        T norm = sqrt(sum);

        vector<T> result(x_size);
        
        for (int i = 0; i < x_size; ++i)
        {
            result[i] = x[i] / norm;
        }
        
        return result;
    }
} // namespace Vector_Functions

#endif
