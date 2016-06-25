set(library_name discretization_spatial)
set(library_dependencies external utilities)

include_directories(${global_include_directories})

file(GLOB src *.cc *.hh)
add_library(${library_name} ${src})
target_link_libraries(${library_dependencies})
install(TARGETS ${library_name} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
