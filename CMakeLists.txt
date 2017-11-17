cmake_minimum_required(VERSION 2.6)

PROJECT(ZRSTD C CXX)

SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/" ${CMAKE_MODULE_PATH})
#set(CMAKE_VERBOSE_MAKEFILE on )

# Where to build shared libs
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Where to build exes
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Global linker options
SET(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-as-needed")

# Set build type (Debug, Release, RelWithDebInfo)
# If the user specifies -DCMAKE_BUILD_TYPE on the command line, take their definition
# and dump it in the cache along with proper documentation, otherwise set CMAKE_BUILD_TYPE
# to Debug prior to calling PROJECT()
IF (DEFINED CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of
	build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug
	Release RelWithDebInfo MinSizeRel.")
ELSE ()
    SET(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build,
	options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release
	RelWithDebInfo MinSizeRel.")
ENDIF ()

# Compile options
ADD_DEFINITIONS(-std=c++11 -O2 -g -ffloat-store -pthread)
#ADD_DEFINITIONS(-std=c++11 -O3 -g -fopenmp -w -ffloat-store)
# turn no return-value warnings into errors
ADD_DEFINITIONS(-w)
#ADD_DEFINITIONS(-Wall)
#ADD_DEFINITIONS(-Werror=return-type)
#ADD_DEFINITIONS(-Wno-reorder)

# Dependencies
FIND_PACKAGE(Armadillo REQUIRED)
FIND_PACKAGE(Threads REQUIRED)
FIND_PACKAGE(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif ()


FIND_PACKAGE(Boost REQUIRED)

FIND_PACKAGE(CUDA REQUIRED)

ADD_SUBDIRECTORY(ZREval)