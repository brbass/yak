cmake_minimum_required (VERSION 3.1)

project(Yak)

####################
# Initialize lists #
####################

# Finds all header files recursively in given directories
macro(header_directories return_list)
    file(GLOB_RECURSE header_files packages/*.hh lib/*.hpp)
    set(header_directories "")
    foreach(header_file ${header_files})
        get_filename_component(header_directory ${header_file} DIRECTORY)
        set(header_directories ${header_directories} ${header_directory})
    endforeach()
    list(REMOVE_DUPLICATES header_directories)
    set(${return_list} ${header_directories})
endmacro()

# Lists of include and link values
set(global_include_directories "")
set(global_trilinos_include_directories "")
set(global_trilinos_link_directories "")
set(global_trilinos_link_libraries "")

header_directories(global_include_directories)

####################
# Compiler options #
####################

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "-std=c++11 -o3")
set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/bin)

#############
# GSL Flags #
#############

# find_package(GSL REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -L/home/software/rhel6/gsl/1.15/lib -lgsl -lgslcblas -I/home/software/rhel6/gsl/1.15/include")

#################
# Find Trilinos #
#################

set(trilinosPackages Amesos AztecOO Epetra)

foreach(trilinosPackage ${trilinosPackages})
    # Find package
    find_package(${trilinosPackage})
    
    # Set compiler flags
    set(CMAKE_CXX_COMPILER ${${trilinosPackage}_CXX_COMPILER} )
    set(CMAKE_C_COMPILER ${${trilinosPackage}_C_COMPILER} )
    set(CMAKE_CXX_FLAGS  "${${trilinosPackage}_CXX_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}")
    set(CMAKE_C_FLAGS  "${${trilinosPackage}_C_COMPILER_FLAGS} ${CMAKE_C_FLAGS}")
    
    # Get lists of include and link values
    set(global_trilinos_include_directories ${global_trilinos_include_directories} ${${trilinosPackage}_INCLUDE_DIRS} ${${trilinosPackage}_TPL_INCLUDE_DIRS})
    set(global_trilinos_link_directories ${global_trilinos_link_directories} ${${trilinosPackage}_LIBRARY_DIRS} ${${trilinosPackage}_TPL_LIBRARY_DIRS})
    set(global_trilinos_link_libraries ${global_trilinos_link_libraries} ${${trilinosPackage}_LIBRARIES} ${${trilinosPackage}_TPL_LIBRARIES})
endforeach()

########################
# Include source files #
########################

# Make lists global
set(global_include_directories ${global_include_directories} CACHE INTERNAL "global_include_directories")
set(global_trilinos_directories ${global_trilinos_directories} CACHE INTERNAL "global_trilinos_directories")
set(global_trilinos_libraries ${global_trilinos_libraries} CACHE INTERNAL "global_trilinos_libraries")

set(subdirectories packages)
foreach(subdirectory ${subdirectories})
    add_subdirectory(${subdirectory})
endforeach()
