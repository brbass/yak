set(subdirectories fe rbf solid_geometry)

set(library_name discretization_spatial)
set(library_directories src)
set(library_objects discretization_spatial_solid discretization_spatial_fe discretization_spatial_rbf)
set(library_include_directories ${global_include_directories} ${global_trilinos_include_directories})
set(library_link_directories ${global_trilinos_link_directories})
set(library_dependencies external utilities ${global_trilinos_link_libraries})

set(src "")
foreach(library_directory ${library_directories})
    file(GLOB temp ${library_directory}/*.cc ${library_directory}/*.hh)
    set(src ${src} ${temp})
endforeach()
foreach(library_object ${library_objects})
   set(src ${src} $<TARGET_OBJECTS:${library_object}>)
endforeach()

include_directories(${library_include_directories})
link_directories(${library_link_directories})

add_library(${library_name} ${src})

target_link_libraries(${library_name} ${library_dependencies})

foreach(subdirectory ${subdirectories})
    add_subdirectory(${subdirectory})
endforeach()

install(TARGETS ${library_name} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
