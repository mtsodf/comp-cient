cmake_minimum_required(VERSION 2.8)

SET( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -mkl" )
add_executable(Laplaciano Laplaciano.c solvers.c)
add_executable(LaplacianoNaoLinear LaplacianoNaoLinear.c solvers.c)