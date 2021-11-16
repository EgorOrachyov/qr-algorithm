cmake_minimum_required(VERSION 3.19)
project(qr_algo_omp)

set(CMAKE_CXX_STANDARD 17)
add_executable(qr_algo src/qr_algorithm.cpp)

target_compile_options(qr_algo PUBLIC /openmp)

if (${CMAKE_BUILD_TYPE} STREQUAL Release)
    message(STATUS "Configure flags for Release build")
    target_compile_options(qr_algo PUBLIC /O2)
    target_compile_options(qr_algo PUBLIC /fp:precise)
endif ()