cmake_minimum_required(VERSION 3.5)

project(mtk)
enable_language(Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-ffree-form  -fimplicit-none")
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

#set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
#set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${opt}")
set(CMAKE_Fortran_FLAGS_DEBUG "${debug}")
set(CMAKE_Fortran_FLAGS_RELEASE "${release}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect} -cpp -ffree-line-length-none")

#
# Compile.
#
file(GLOB_RECURSE sources  src/*.f90)
message(STATUS ${CMAKE_Fortran_FLAGS_DEBUG})
message(STATUS ${CMAKE_Fortran_FLAGS_RELEASE})
add_executable(simulate ${sources})
