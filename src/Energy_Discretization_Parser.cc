#include "Energy_Discretization_Parser.hh"

Energy_Discretization_Parser::
Energy_Discretization_Parser(pugi::xml_node &input_file):
    Parser(input_file)
{
    pugi::xml_node energy = input_file.child("energy_discretization");
    
    int number_of_groups = child_value<int>(energy, "number_of_groups");
    
    energy_ = make_shared<Energy_Discretization>(number_of_groups);
}
