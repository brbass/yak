set(library_name discretization_spatial_solid)
# set(library_dependencies external utilities)

include_directories(${global_include_directories})

file(GLOB src *.cc *.hh)
add_library(${library_name} OBJECT ${src})

# target_link_libraries(${library_name} ${library_dependencies})
# install(TARGETS ${library_name} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
