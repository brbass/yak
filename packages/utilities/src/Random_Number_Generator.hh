#ifndef Random_Number_Generator_hh
#define Random_Number_Generator_hh

#include <random>
#include <vector>

using std::vector;

/*
  Simple class to generate several types or sequences of random numbers
*/
class Random_Number_Generator
{
public:

    // Constructor
    Random_Number_Generator(double lower_bound,
                            double upper_bound);

    // Returns a random integer
    int random_int();

    // Returns a random double
    double random_double();

    // Returns a sequence of random integers
    vector<int> random_int_vector(unsigned number_of_elements);

    // Returns a sequence of random doubles
    vector<double> random_double_vector(unsigned number_of_elements);

private:

    // Get seed for random number generator
    // Uses time since epoch
    unsigned get_seed();
    
    std::default_random_engine generator;
    std:: uniform_real_distribution<double> real_distribution;
    std:: uniform_int_distribution<int> int_distribution;

};

#endif
