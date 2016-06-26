set(library_name driver)
set(library_dependencies external utilities discretization data solver)

include_directories(${global_include_directories})

file(GLOB src *.cc *.hh)
add_library(${library_name} ${src})
target_link_libraries(${library_name} ${library_dependencies})
install(TARGETS ${library_name} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
