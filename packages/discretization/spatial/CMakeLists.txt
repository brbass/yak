set(subdirectories fe rbf solid_geometry)

set(library_name discretization_spatial)
set(library_directories fe/src rbf/src solid_geometry/src src)
set(library_include_directories ${global_include_directories} ${global_trilinos_include_directories})
set(library_link_directories ${global_trilinos_link_directories})
set(library_dependencies utilities ${global_trilinos_link_libraries})

set(src "")
foreach(library_directory ${library_directories})
    file(GLOB temp ${library_directory}/*.cc ${library_directory}/*.hh)
    set(src ${src} ${temp})
endforeach()

include_directories(${library_include_directories})
link_directories(${library_link_directories})

add_library(${library_name} ${src})

target_link_libraries(${library_name} ${library_dependencies})

install(TARGETS ${library_name} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

foreach(subdirectory ${subdirectories})
    add_subdirectory(${subdirectory})
endforeach()
