cmake_minimum_required(VERSION 3.0)

project(md)

find_package(MPI REQUIRED)
include_directories(${MPI_Fortran_INCLUDE_PATH})

enable_language(Fortran)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_Fortran_COMPILER mpiifort)
set(CMAKE_Fortran_FLAGS "-g -traceback -cpp -lmkl_intel_lp64 -lmkl_sequential -lmkl_core")

set(source_files
  src/multi_diag.F90
)

add_library(md ${source_files})
target_link_libraries(md ${MPI_Fortran_LIBRARIES})

add_executable(test.exe test/main.F90)
target_link_libraries(test.exe md ${MPI_Fortran_LIBRARIES})