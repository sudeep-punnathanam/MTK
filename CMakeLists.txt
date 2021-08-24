cmake_minimum_required(VERSION 3.5)

project(mtk)
enable_language(Fortran)
set(CMAKE_Fortran_COMPILER "/usr/local/bin/gfortran")
set(CMAKE_C_COMPILER "/usr/local/bin/gcc-11")
set(CMAKE_CXX_COMPILER "/usr/local/bin/g++-11")

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-cpp -ffree-line-length-none -ffree-form  -fimplicit-none")
    set(debug "-O0 -fbounds-check")
    set(release "-O1 -ffast-math -funroll-loops")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-stand f08 -free -implicitnone")
    set(debug "-O0 -check bounds")
    set(release "-O3")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
    set(bounds "-C")
endif()

set(CMAKE_Fortran_FLAGS "${dialect}")
set(CMAKE_Fortran_FLAGS_DEBUG "${debug}")
set(CMAKE_Fortran_FLAGS_RELEASE "${release}")

#
# Compile.
#
#file(GLOB_RECURSE sources  src/*.f90)
file(GLOB sources  src/*.f90)
message(STATUS ${CMAKE_Fortran_FLAGS_DEBUG})
message(STATUS ${CMAKE_Fortran_FLAGS_RELEASE})
message(STATUS ${CMAKE_Fortran_FLAGS})
message(STATUS ${CMAKE_Fortran_COMPILER})
message(STATUS ${CMAKE_C_COMPILER})
message(STATUS ${CMAKE_CXX_COMPILER})
add_executable(simulate ${sources})
