set(SRC
  kinds.f90
  parallel_env.f90
  params.f90
  lattice_geometry.f90
  sun_group.f90
  gauge_field.f90
  fermion_wilson.f90
  fermion_staggered.f90
  hamiltonian.f90
  main.f90
)

add_executable(lqcd ${SRC})

# Compiler options
if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  target_compile_options(lqcd PRIVATE -Wall -Wextra -ffree-line-length-none)
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  target_compile_options(lqcd PRIVATE -warn all)
endif()

# OpenMP
if(ENABLE_OPENMP)
  find_package(OpenMP COMPONENTS Fortran)
  if(OpenMP_Fortran_FOUND)
    target_link_libraries(lqcd PRIVATE OpenMP::OpenMP_Fortran)
    target_compile_definitions(lqcd PRIVATE USE_OPENMP)
  endif()
endif()

# MPI
if(ENABLE_MPI)
  find_package(MPI COMPONENTS Fortran)
  if(MPI_Fortran_FOUND)
    target_link_libraries(lqcd PRIVATE MPI::MPI_Fortran)
    target_compile_definitions(lqcd PRIVATE USE_MPI)
  endif()
endif()

# Ensure preprocessing is enabled
set_property(TARGET lqcd PROPERTY Fortran_PREPROCESS ON)
