set(library_name discretization_spatial_rbf)

include_directories(${global_include_directories})

file(GLOB src *.cc *.hh)
add_library(${library_name} OBJECT ${src})
