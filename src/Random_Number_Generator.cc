#include "Random_Number_Generator.hh"

#include <chrono>
#include <vector>

using namespace std;

Random_Number_Generator::
Random_Number_Generator(double lower_bound,
                        double upper_bound):
    generator(get_seed()),
    int_distribution(lower_bound, upper_bound),
    real_distribution(lower_bound, upper_bound)
{
}

int Random_Number_Generator::
random_int()
{
    return int_distribution(generator);
}

double Random_Number_Generator::
random_double()
{
    return real_distribution(generator);
}

unsigned Random_Number_Generator::
get_seed()
{
    using namespace chrono;

    duration<int> time = duration_cast<seconds>(high_resolution_clock::now().time_since_epoch());
    
    return time.count();
}

vector<int> Random_Number_Generator::
random_int_vector(unsigned number_of_elements)
{
    vector<int> vec(number_of_elements);
    
    for (unsigned i=0; i<number_of_elements; ++i)
    {
        vec[i] = int_distribution(generator);
    }

    return vec;
}

vector<double> Random_Number_Generator::
random_double_vector(unsigned number_of_elements)
{
    vector<double> vec(number_of_elements);
    
    for (unsigned i=0; i<number_of_elements; ++i)
    {
        vec[i] = real_distribution(generator);
    }

    return vec;
}
