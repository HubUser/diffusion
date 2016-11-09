cmake_minimum_required(VERSION 3.6)
set (PROJECT fort_lib) 
set (SOURCES radbelInput_tet2_kp.doubleQ.f)

project (${PROJECT} Fortran)
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
        "Choose the type of build, options are: None Debug Release."
        FORCE)
endif (NOT CMAKE_BUILD_TYPE)
# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
# gfortran
set (CMAKE_Fortran_FLAGS_RELEASE "-Wall -frecursive -O3 -fPIC")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-Wall -Og -pg -fbacktrace -fPIC")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
    # ifort (untested)
    set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -Og -g")
elseif (Fortran_COMPILER_NAME MATCHES "g77")
    # g77
    set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
    message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
    message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
    message ("No optimized Fortran compiler flags are known, we just try -O2...")
    set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-Og -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

add_library(${PROJECT} ${SOURCES})
