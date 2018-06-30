cmake_minimum_required(VERSION 3.8)
project(AsFem)

set(CMAKE_CXX_STANDARD 14)

# For debug mode
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -fopenmp -O2")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")

# For release mode
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fopenmp -Ofast -march=corei7-avx")
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# For PETSC
set(petsc "/home/by/Programs/petsc3.9.2")
include_directories("${petsc}/include")
link_libraries("${petsc}/lib/libpetsc.so")

# For libmpi
set(mpi "/opt/moose/mpich-3.2/gcc-7.3.0")
link_libraries("${mpi}/lib/libmpi.so")
include_directories("${mpi}/include")


# For Welcome header file and main.cpp
set(HEADER include/Welcome.h)
set(SRC src/main.cpp)

# For string utils 
set(HEADER ${HEADER} include/StringUtils/StringUtils.h)
set(SRC ${SRC} src/StringUtils/RemoveSpace.cpp src/StringUtils/SplitNum.cpp)
set(SRC ${SRC} src/StringUtils/SplitNumAfter.cpp src/StringUtils/SplitStr.cpp)
set(SRC ${SRC} src/StringUtils/StrToLower.cpp src/StringUtils/StrToUpper.cpp)
set(SRC ${SRC} src/StringUtils/StrVecToLower.cpp src/StringUtils/StrVecToUpper.cpp)
set(SRC ${SRC} src/StringUtils/IsBracketMatch.cpp src/StringUtils/GotoLine.cpp)

# For mesh utils 
set(HEADER ${HEADER} include/Mesh/MeshBase.h include/Mesh/Mesh1D.h include/Mesh/Mesh2D.h)
set(HEADER ${HEADER} include/Mesh/Mesh.h)

set(SRC ${SRC} src/Mesh/Mesh1D.cpp src/Mesh/Mesh2D.cpp)
set(SRC ${SRC} src/Mesh/Mesh.cpp)
# For gauss point
set(HEADER ${HEADER} include/FE/Qrule.h)
set(SRC ${SRC} src/FE/Qrule.cpp)

# For DofHandler
set(HEADER ${HEADER} include/DofHandler/DofHandler.h)
set(SRC ${SRC} src/DofHandler/DofHandler.cpp)

# For Equation System
set(HEADER ${HEADER} include/EquationSystem/EquationSystem.h)
set(SRC ${SRC} src/EquationSystem/EquationSystem.cpp)
# For SolverSystem
set(HEADER ${HEADER} include/SolverSystem/SolverSystem.h)
set(SRC ${SRC} src/SolverSystem/SolverSystem.cpp)

# For InputSystem
set(HEADER ${HEADER} include/InputSystem/InputSystem.h)
set(SRC ${SRC} src/InputSystem/InputSystem.cpp src/InputSystem/ReadDofsName.cpp)

#--> for mesh
set(SRC ${SRC} src/InputSystem/ReadMesh.cpp)

# For FESystem
set(HEADER ${HEADER} include/FESystem/FESystem.h)
set(SRC ${SRC} src/FESystem/FESystem.cpp include/BCSystem/BCSystem.h)
set(HEADER ${HEADER} include/FE/FE.h)
set(SRC ${SRC} src/FE/FE.cpp)

# For Kernels
set(HEADER ${HEADER} include/Kernels/Kernel.h include/Kernels/BoundaryKernel.h)

set(HEADER ${HEADER} include/Kernels/Laplace.h)
set(SRC ${SRC} src/Kernels/Laplace.cpp)


add_executable(ASFEM ${HEADER} ${SRC} )
