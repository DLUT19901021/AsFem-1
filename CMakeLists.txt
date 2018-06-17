cmake_minimum_required(VERSION 3.8)
project(AsFem)

set(CMAKE_CXX_STANDARD 14)

# For debug mode
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -fopenmp -O2")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")

# For release mode
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fopenmp -Ofast -march=corei7-avx")
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

# For PETSC
set(petsc "/home/by/Programs/petsc3.9.2")
include_directories("${petsc}/include")
link_libraries("${petsc}/lib/libpetsc.so")

# For libmpi
set(mpi "/opt/moose/mpich-3.2/gcc-7.3.0")
link_libraries("${mpi}/lib/libmpi.so")
include_directories("${mpi}/include")

add_executable(ASFEM src/main.cpp src/Mesh/MeshBase.h src/Mesh/Mesh1D.cpp src/Mesh/Mesh1D.h src/Welcome.h src/Mesh/Mesh2D.cpp src/Mesh/Mesh2D.h src/StringUtils/StringUtils.h src/StringUtils/RemoveSpace.cpp src/StringUtils/StrToLower.cpp src/StringUtils/StrVecToLower.cpp src/StringUtils/StrToUpper.cpp src/StringUtils/StrVecToUpper.cpp src/StringUtils/SplitStr.cpp src/StringUtils/SplitNum.cpp src/StringUtils/SplitNumAfter.cpp src/Mesh/Mesh.cpp src/Mesh/Mesh.h)
