cmake_minimum_required(VERSION 3.6)

set (PROJECT Diffusion) 
set (SOURCES main.cpp
             radbel.cpp
)

#set (HEADERS headers/lab1.h
#             headers/controller.h
#)
project (${PROJECT} C CXX Fortran)
include (CMakeAddFortranSubdirectory)
include(FortranCInterface)
find_package(PythonLibs)
find_package(Boost)
FortranCInterface_HEADER(FC.h MACRO_NAMESPACE FC_ SYMBOLS localdtet2dst alpha)
#set (CMAKE_BUILD_TYPE Debug) 
set (CMAKE_BUILD_TYPE Release) 

include_directories (${PYTHON_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
#include_directories (${PYTHON_INCLUDE_DIRS})

#link_directories (${CMAKE_PREFIX_PATH}/src/plugins/platforms/xcb/xcb-static)

#add_definitions (${QT_DEFINITIONS} -Werror -Wall -fPIC)

cmake_add_fortran_subdirectory ("src_fortran/" 
    PROJECT fort_lib 
    ARCHIVE_DIR "src_fortran/lib"
    RUNTIME_DIR "src_fortran/bin"
    LIBRARIES fort_lib
    LINK_LIBRARIES
    LINK_LIBS fort_lib
    NO_EXTERNAL_INSTALL)

add_library (${PROJECT} MODULE ${SOURCES})
#add_library (${PROJECT} SHARED ${SOURCES})
set_property(TARGET ${PROJECT} PROPERTY CXX_STANDARD 11)
target_link_libraries (${PROJECT} fort_lib)
