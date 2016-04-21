#ifndef Transport_Problem_hh
#define Transport_Problem_hh

#include <memory>
#include "pugixml.hpp"

#include "Solver.hh"

using std::shared_ptr;

class Transport_Problem
{
public:
    
    enum Problem_Type
    {
        STEADY_STATE,
        K_EIGENVALUE,
        TIME_DEPENDENT
    };
    
    Transport_Problem(Problem_Type problem_type,
                      shared_ptr<Solver> solver);
    
    void solve();

    void output(pugi::xml_node &output_node) const;

private:
    
    Problem_Type problem_type_;
    shared_ptr<Solver> solver_;
    
    double k_eigenvalue_;
    vector<double> phi_;
    vector<double> psi_;
};

#endif
