cmake_minimum_required(VERSION 3.19)
project(cppprojct)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_C_COMPILER mpicc)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/usr/lib/eigen")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")

# VERSIONS:
# OpenMP Unchecked Version
# Eigen Unchecked Version
# g++ 7.50
# Ubuntu 18.04
# OpenMPI 2.1.1
# Boost 1.77.0

# ----------Build Type----------
#set(CMAKE_BUILD_TYPE Debug) # DEBUG SHOULD BE ON FOR PROFILING
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O1") # DEBUG SHOULD REQUIRE 00, PRODUCTION 03
# MacOS should additionally set -fno-omit-frame-pointer

# this is to set verbosity :-)
#add_definitions(-DVERBOSE)


# Looking for OpenMP
find_package(OpenMP REQUIRED)

# Looking for MPI
find_package(MPI REQUIRED)

# Looking for Boost: https://cmake.org/cmake/help/latest/module/FindBoost.html
set(Boost_NO_BOOST_CMAKE    ON)
find_package(Boost 1.77.0 REQUIRED COMPONENTS mpi graph_parallel system serialization)
include_directories(${Boost_INCLUDE_DIRS})


#add_definitions(-DVERBOSE)
#                                                                       (max)Nthreads                 Nprocs
# Run with  "mpiexec -x NNODES=15 -x TEST=1 -x SEED=1234 -x OMP_THREAD_LIMIT=5  -x  OMP_NESTED=true -n 3  cmake-build-debug/cppprojct"
add_executable(cppprojct main.cpp Utils/timers.cpp Utils/timers.h GraphClasses/ErdosRenyiGraph.cpp GraphClasses/ErdosRenyiGraph.h Utils/error.cpp Utils/error.h GraphClasses/GeneralGraph.cpp GraphClasses/GeneralGraph.h GraphClasses/CliqueGraph.h GraphClasses/CliqueGraph.cpp Simulation/Simulation.cpp  Tests/graph-test-init.cpp Tests/graph-test-init.h  Utils/reproductibility.cpp Utils/reproductibility.h Utils/adequate_synchronization.h Utils/adequate_synchronization.cpp GraphClasses/RingGraph.h GraphClasses/RingGraph.cpp  Tests/graph-test-singlestep-evolution.h Tests/graph-test-singlestep-evolution.cpp Utils/global_standard_messages.h Solvers/GeneralSolver.h Solvers/GeneralSolver.cpp Utils/memory_management.h Utils/memory_management.cpp Tests/solvers-test-init.h Tests/solvers-test-init.cpp DifferentialEquations/GeneralDifferentialEquation.h DifferentialEquations/GeneralDifferentialEquation.cpp Solvers/EulerSolver.h Solvers/EulerSolver.cpp Utils/differential_equations_aux.h Utils/differential_equations_aux.cpp DifferentialEquations/NoiselessKuramoto.h DifferentialEquations/NoiselessKuramoto.cpp Solvers/RungeKuttaSolver.cpp Solvers/RungeKuttaSolver.h  Utils/typed_combinations.h Utils/msleep.h Utils/msleep.cpp GraphClasses/GraphFunctions.h GraphClasses/GraphFunctions.cpp Communication/CommunicationFunctions.cpp Communication/CommunicationFunctions.h Utils/HelperClasses.h Utils/HelperClasses.cpp Utils/print_init.h Utils/print_warnings.h Utils/print_init.cpp Utils/print_warnings.cpp Utils/parallel_sanitizer.h Tests/long-singlestep-run.h Tests/long-singlestep-run.cpp macros/macros.h Utils/certify_compliant.h Utils/certify_compliant.cpp)
target_link_libraries(cppprojct PRIVATE OpenMP::OpenMP_CXX ${Boost_LIBRARIES})
