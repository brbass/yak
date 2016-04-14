#ifndef Random_Number_Generator_hh
#define Random_Number_Generator_hh

#include <random>
#include <vector>

using std::vector;

class Random_Number_Generator
{
private:
    
    std::default_random_engine generator;
    std:: uniform_real_distribution<double> real_distribution;
    std:: uniform_int_distribution<int> int_distribution;

    unsigned get_seed();

public:

    Random_Number_Generator(double lower_bound,
                            double upper_bound);
    
    int random_int();
    double random_double();
    vector<int> random_int_vector(unsigned number_of_elements);
    vector<double> random_double_vector(unsigned number_of_elements);
};

#endif
