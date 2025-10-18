cmake_minimum_required(VERSION 3.15)
project(LatticeQCD LANGUAGES Fortran)

option(ENABLE_MPI "Enable MPI parallelization" ON)
option(ENABLE_OPENMP "Enable OpenMP parallelization" ON)

# Enable Fortran preprocessing for conditional compilation
set(CMAKE_Fortran_PREPROCESS ON)

# Prefer Release by default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

add_subdirectory(src)
