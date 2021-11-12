cmake_minimum_required(VERSION 3.19)
project(qr_decomposition)

function(set_compile_options target)
    target_compile_options(${target} PUBLIC /openmp)
endfunction()

set(CMAKE_CXX_STANDARD 17)
set(QR_SOURCES src/main.cpp src/math.hpp src/qr_omp.hpp src/qr_seq.hpp)
add_executable(qr_algo ${QR_SOURCES})
set_compile_options(qr_algo)
