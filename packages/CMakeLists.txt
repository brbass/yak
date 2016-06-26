# set(subdirectories data discretization driver operator solver sweep utilities yak)
set(subdirectories external utilities discretization data operator sweep solver driver yak)

foreach(subdirectory ${subdirectories})
    add_subdirectory(${subdirectory})
endforeach()
