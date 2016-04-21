#ifndef Parser_hh
#define Parser_hh

#include <memory>

#include "pugixml.hpp"

#include "XML_Child_Value.hh"

using std::shared_ptr;

template<class T>
class Parser
{
public:

    Parser(pugi::xml_node &input_file):
        input_file_(input_file)
    {
    }

    virtual shared_ptr<T> get_ptr() = 0;
    
protected:
    
    pugi::xml_node &input_file_;
};

#endif
