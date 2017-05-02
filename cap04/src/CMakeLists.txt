cmake_minimum_required(VERSION 3.4)

project (RayleightCoeffIter)

# this one is important
SET(CMAKE_SYSTEM_NAME Generic)

include(CMakeForceCompiler)
CMAKE_FORCE_C_COMPILER(icc "Intel C Compiler")
CMAKE_FORCE_CXX_COMPILER(icpc "Intel C++ Compiler")
CMAKE_FORCE_Fortran_COMPILER(ifort "Intel Fortran Compiler")

SET (CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -mkl)

add_executable(RayleightCoeffIter RayleightCoeffIter.cpp)
