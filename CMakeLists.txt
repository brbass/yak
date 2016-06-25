cmake_minimum_required (VERSION 3.1)

project(Yak)

####################
# Compiler options #
####################

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "-std=c++11 -o3")
# set(CMAKE_INSTALL_PREFIX /home/brbass/research/yak/bin)

#############
# GSL Flags #
#############

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -L/home/software/rhel6/gsl/1.15/lib -lgsl -lgslcblas -I/home/software/rhel6/gsl/1.15/include")

#################
# Find Trilinos #
#################

set(TrilinosPackages Amesos AztecOO Epetra)

foreach(TrilinosPackage ${TrilinosPackages})
    find_package(${TrilinosPackage})
    set(CMAKE_CXX_COMPILER ${${TrilinosPackage}_CXX_COMPILER} )
    set(CMAKE_C_COMPILER ${${TrilinosPackage}_C_COMPILER} )
    set(CMAKE_CXX_FLAGS  "${${TrilinosPackage}_CXX_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}")
    set(CMAKE_C_FLAGS  "${${TrilinosPackage}_C_COMPILER_FLAGS} ${CMAKE_C_FLAGS}")
    include_directories(${${TrilinosPackage}_INCLUDE_DIRS} ${${TrilinosPackage}_TPL_INCLUDE_DIRS})
    link_directories(${${TrilinosPackage}_LIBRARY_DIRS} ${${TrilinosPackage}_TPL_LIBRARY_DIRS})
endforeach(TrilinosPackage)

########################
# Include source files #
########################

file(GLOB SOURCES "src/*.cc" "lib/*.cpp")
set(MAIN "main/Yak.cc")

include_directories(src lib)

add_executable(yak ${SOURCES} ${MAIN})

foreach(TrilinosPackage ${TrilinosPackages})
    target_link_libraries(yak ${${TrilinosPackage}_LIBRARIES} ${${TrilinosPackage}_TPL_LIBRARIES})
endforeach(TrilinosPackage)

# install(TARGETS yak DESTINATION .)

