cmake_minimum_required(VERSION 3.19)
project(cppprojct)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER mpic++)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/usr/lib/eigen")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")

# VERSIONS:
# OpenMP Unchecked Version
# Eigen Unchecked Version
# g++ 7.50
# Ubuntu 18.04
# OpenMPI 2.1.1
# Boost 1.77.0

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Looking for OpenMP
find_package(OpenMP REQUIRED)

# Looking for MPI
find_package(MPI REQUIRED)

# Looking for Boost: https://cmake.org/cmake/help/latest/module/FindBoost.html
set(Boost_NO_BOOST_CMAKE    ON)
find_package(Boost 1.77.0 REQUIRED COMPONENTS mpi graph_parallel system serialization)
include_directories(${Boost_INCLUDE_DIRS})

# OPTIMIZATION
# Please use -02 or -03 for compilation in production

# Run with mpirun or mpiexec
add_executable(cppprojct main.cpp Utils/timers.cpp Utils/timers.h GraphClasses/ErdosRenyiGraph.cpp GraphClasses/ErdosRenyiGraph.h Utils/error.cpp Utils/error.h GraphClasses/GeneralGraph.cpp GraphClasses/GeneralGraph.h GraphClasses/CliqueGraph.h GraphClasses/CliqueGraph.cpp Simulation/Simulation.cpp  Tests/graph-test-init.cpp Tests/graph-test-init.h  Utils/reproductibility.cpp Utils/reproductibility.h Utils/adequate_synchronization.h Utils/adequate_synchronization.cpp GraphClasses/RingGraph.h GraphClasses/RingGraph.cpp  Tests/graph-test-singlestep-evolution.h Tests/graph-test-singlestep-evolution.cpp Utils/global_standard_messages.h Solvers/GeneralSolver.h Solvers/GeneralSolver.cpp Utils/memory_management.h Utils/memory_management.cpp Tests/solvers-test-init.h Tests/solvers-test-init.cpp DifferentialEquations/GeneralDifferentialEquation.h DifferentialEquations/GeneralDifferentialEquation.cpp Solvers/EulerSolver.h Solvers/EulerSolver.cpp Utils/differential_equations_aux.h Utils/differential_equations_aux.cpp DifferentialEquations/NoiselessKuramoto.h DifferentialEquations/NoiselessKuramoto.cpp Solvers/RungeKuttaSolver.cpp Solvers/RungeKuttaSolver.h GraphClasses/GraphFunctions.h GraphClasses/GraphFunctions.cpp Utils/typed_combinations.h Utils/msleep.h Utils/msleep.cpp)
target_link_libraries(cppprojct PRIVATE OpenMP::OpenMP_CXX ${Boost_LIBRARIES})
