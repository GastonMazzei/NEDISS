cmake_minimum_required(VERSION 3.19)
project(cppprojct)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_COMPILER mpic++)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/usr/lib/eigen")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L/usr/lib/eigen")


find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

add_executable(cppprojct main.cpp utils/timers.cpp utils/timers.h classes/graph.cpp classes/graph.h utils/error.cpp utils/error.h)
target_link_libraries(cppprojct PRIVATE OpenMP::OpenMP_CXX)

