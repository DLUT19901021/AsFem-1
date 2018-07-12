cmake_minimum_required(VERSION 3.8)
project(AsFem)

set(CMAKE_CXX_STANDARD 14)

# For debug mode
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -fopenmp -O2")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")

# For release mode
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fopenmp -Ofast -march=corei7-avx")
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

###################################################
### Please edit the PETSc and MPI library path  ###
###  based on your own installtion              ###
###################################################
# For PETSC
set(petsc "/home/by/Programs/petsc3.9.2")
include_directories("${petsc}/include")
link_libraries("${petsc}/lib/libpetsc.so")
# For libmpi
set(mpi "/opt/moose/mpich-3.2/gcc-7.3.0")
link_libraries("${mpi}/lib/libmpi.so")
include_directories("${mpi}/include")
###################################################


##################################################
##################################################
### Do not edit the following lines !!!        ###
### Do not edit the following lines !!!        ###
### Do not edit the following lines !!!        ###
##################################################
##################################################
# For Welcome header file and main.cpp
set(HEADER include/Welcome.h)
set(SRC src/main.cpp)

# For string utils 
#set(HEADER ${HEADER} include/StringUtils/StringUtils.h)
#set(SRC ${SRC} src/StringUtils/RemoveSpace.cpp src/StringUtils/SplitNum.cpp)
#set(SRC ${SRC} src/StringUtils/SplitNumAfter.cpp src/StringUtils/SplitStr.cpp)
#set(SRC ${SRC} src/StringUtils/StrToLower.cpp src/StringUtils/StrToUpper.cpp)
#set(SRC ${SRC} src/StringUtils/StrVecToLower.cpp src/StringUtils/StrVecToUpper.cpp)
#set(SRC ${SRC} src/StringUtils/IsBracketMatch.cpp src/StringUtils/GotoLine.cpp src/main.cpp)

# For mesh utils 
set(HEADER ${HEADER} include/Mesh/Mesh.h)
set(SRC ${SRC} src/Mesh/Mesh.cpp src/Mesh/CreateMesh.cpp)

add_executable(ASFEM ${HEADER} ${SRC})

