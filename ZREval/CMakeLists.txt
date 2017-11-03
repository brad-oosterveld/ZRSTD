cmake_minimum_required(VERSION 2.6)
#set(CMAKE_VERBOSE_MAKEFILE ON)
PROJECT(ZREval)

include_directories(./cudacompare)

SET(SOURCES
        ZREval.cpp
        MFCCs.cpp
        Preprocessing.cpp
        ComparisonMethod.cpp
        Cluster.cpp
        TranscriptionUnit.cpp
        )

SET(CUDA_SOURCES
        ./cudacompare/cudacompare.cu
        )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,--no-as-needed ") # need this because of a gcc bug or something
set(GENCODES "${GENCODES} -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=sm_35 -gencode arch=compute_37,code=sm_37 -gencode arch=compute_50,code=sm_50 -gencode arch=compute_52,code=sm_52")
set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -std=c++11 -O3  -g -G -DDEBUG=6 ${GENCODES}") # need this because of a gcc bug or something

add_executable(ZREval main.cpp ${SOURCES})

cuda_add_library(cudacompare ${CUDA_SOURCES})

target_link_libraries(ZREval
        armadillo
        cudacompare
        )