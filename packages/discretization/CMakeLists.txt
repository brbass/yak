set(subdirectories angular energy spatial)

foreach(subdirectory ${subdirectories})
    add_subdirectory(${subdirectory})
endforeach()
