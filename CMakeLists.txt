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
add_executable(cppprojct main.cpp utils/timers.cpp utils/timers.h GraphClasses/ErdosRenyiGraph.cpp GraphClasses/ErdosRenyiGraph.h utils/error.cpp utils/error.h GraphClasses/GeneralGraph.cpp GraphClasses/GeneralGraph.h GraphClasses/CliqueGraph.h GraphClasses/CliqueGraph.cpp Simulation/Simulation.cpp tests/graph-test.cpp tests/graph-test.h utils/reproductibility.cpp utils/reproductibility.h utils/adequate_synchronization.h utils/adequate_synchronization.cpp GraphClasses/RingGraph.h GraphClasses/RingGraph.cpp)
target_link_libraries(cppprojct PRIVATE OpenMP::OpenMP_CXX ${Boost_LIBRARIES})

