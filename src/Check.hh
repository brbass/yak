#ifndef Check_hh
#define Check_hh

#include <string>

namespace ch_ns
{
    using std::string;
    
    void check(string condition,
               string file,
               int line);

    void check(string condition,
               string message,
               string file,
               int line);
}

// Check only happens in debug
#ifdef NDEBUG
#  define Check(cond)
#  define CheckMsg(cond, desc)
#else
#  define Check(cond)                                           \
    if (!(cond)) ch_ns::check(#cond, __FILE__, __LINE__)
#  define CheckMsg(cond, desc)                                  \
    if (!(cond)) ch_ns::check(#cond, desc, __FILE__, __LINE__)
#endif

// Assert always happens
#define Assert(cond)                                            \
    if (!(cond)) ch_ns::check(#cond, __FILE__, __LINE__)
#define AssertMsg(cond, desc)                                   \
    if (!(cond)) ch_ns::check(#cond, desc, __FILE__, __LINE__)

#endif
