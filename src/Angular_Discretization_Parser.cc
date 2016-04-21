#include "Angular_Discretization_Parser.hh"

#include "Gauss_Legendre_Ordinates.hh"

using namespace std;

Angular_Discretization_Parser::
Angular_Discretization_Parser(pugi::xml_node &input_file):
    Parser(input_file)
{
    pugi::xml_node angular = input_file.child("angular_discretization");
    
    int dimension = child_value<int>(angular, "dimension");
    int number_of_moments = child_value<int>(angular, "number_of_moments");
    int number_of_ordinates = child_value<int>(angular, "number_of_ordinates");
    
    angular_ = make_shared<Gauss_Legendre_Ordinates>(dimension,
                                                     number_of_moments,
                                                     number_of_ordinates);
}
