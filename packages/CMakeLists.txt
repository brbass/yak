set(subdirectories external utilities solid_geometry energy_discretization angular_discretization spatial_discretization data operator transport preconditioner solver driver yak)

foreach(subdirectory ${subdirectories})
    add_subdirectory(${subdirectory})
endforeach()
