# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.19

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/clion-2021.1.1/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /opt/clion-2021.1.1/bin/cmake/linux/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/m4zz31/cppprojct

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/m4zz31/cppprojct/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/cppprojct.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/cppprojct.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/cppprojct.dir/flags.make

CMakeFiles/cppprojct.dir/main.cpp.o: CMakeFiles/cppprojct.dir/flags.make
CMakeFiles/cppprojct.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/m4zz31/cppprojct/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/cppprojct.dir/main.cpp.o"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cppprojct.dir/main.cpp.o -c /home/m4zz31/cppprojct/main.cpp

CMakeFiles/cppprojct.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cppprojct.dir/main.cpp.i"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/m4zz31/cppprojct/main.cpp > CMakeFiles/cppprojct.dir/main.cpp.i

CMakeFiles/cppprojct.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cppprojct.dir/main.cpp.s"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/m4zz31/cppprojct/main.cpp -o CMakeFiles/cppprojct.dir/main.cpp.s

CMakeFiles/cppprojct.dir/Utils/timers.cpp.o: CMakeFiles/cppprojct.dir/flags.make
CMakeFiles/cppprojct.dir/Utils/timers.cpp.o: ../Utils/timers.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/m4zz31/cppprojct/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/cppprojct.dir/Utils/timers.cpp.o"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cppprojct.dir/Utils/timers.cpp.o -c /home/m4zz31/cppprojct/Utils/timers.cpp

CMakeFiles/cppprojct.dir/Utils/timers.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cppprojct.dir/Utils/timers.cpp.i"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/m4zz31/cppprojct/Utils/timers.cpp > CMakeFiles/cppprojct.dir/Utils/timers.cpp.i

CMakeFiles/cppprojct.dir/Utils/timers.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cppprojct.dir/Utils/timers.cpp.s"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/m4zz31/cppprojct/Utils/timers.cpp -o CMakeFiles/cppprojct.dir/Utils/timers.cpp.s

CMakeFiles/cppprojct.dir/GraphClasses/ErdosRenyiGraph.cpp.o: CMakeFiles/cppprojct.dir/flags.make
CMakeFiles/cppprojct.dir/GraphClasses/ErdosRenyiGraph.cpp.o: ../GraphClasses/ErdosRenyiGraph.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/m4zz31/cppprojct/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/cppprojct.dir/GraphClasses/ErdosRenyiGraph.cpp.o"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cppprojct.dir/GraphClasses/ErdosRenyiGraph.cpp.o -c /home/m4zz31/cppprojct/GraphClasses/ErdosRenyiGraph.cpp

CMakeFiles/cppprojct.dir/GraphClasses/ErdosRenyiGraph.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cppprojct.dir/GraphClasses/ErdosRenyiGraph.cpp.i"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/m4zz31/cppprojct/GraphClasses/ErdosRenyiGraph.cpp > CMakeFiles/cppprojct.dir/GraphClasses/ErdosRenyiGraph.cpp.i

CMakeFiles/cppprojct.dir/GraphClasses/ErdosRenyiGraph.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cppprojct.dir/GraphClasses/ErdosRenyiGraph.cpp.s"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/m4zz31/cppprojct/GraphClasses/ErdosRenyiGraph.cpp -o CMakeFiles/cppprojct.dir/GraphClasses/ErdosRenyiGraph.cpp.s

CMakeFiles/cppprojct.dir/Utils/error.cpp.o: CMakeFiles/cppprojct.dir/flags.make
CMakeFiles/cppprojct.dir/Utils/error.cpp.o: ../Utils/error.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/m4zz31/cppprojct/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/cppprojct.dir/Utils/error.cpp.o"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cppprojct.dir/Utils/error.cpp.o -c /home/m4zz31/cppprojct/Utils/error.cpp

CMakeFiles/cppprojct.dir/Utils/error.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cppprojct.dir/Utils/error.cpp.i"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/m4zz31/cppprojct/Utils/error.cpp > CMakeFiles/cppprojct.dir/Utils/error.cpp.i

CMakeFiles/cppprojct.dir/Utils/error.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cppprojct.dir/Utils/error.cpp.s"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/m4zz31/cppprojct/Utils/error.cpp -o CMakeFiles/cppprojct.dir/Utils/error.cpp.s

CMakeFiles/cppprojct.dir/GraphClasses/GeneralGraph.cpp.o: CMakeFiles/cppprojct.dir/flags.make
CMakeFiles/cppprojct.dir/GraphClasses/GeneralGraph.cpp.o: ../GraphClasses/GeneralGraph.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/m4zz31/cppprojct/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/cppprojct.dir/GraphClasses/GeneralGraph.cpp.o"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cppprojct.dir/GraphClasses/GeneralGraph.cpp.o -c /home/m4zz31/cppprojct/GraphClasses/GeneralGraph.cpp

CMakeFiles/cppprojct.dir/GraphClasses/GeneralGraph.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cppprojct.dir/GraphClasses/GeneralGraph.cpp.i"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/m4zz31/cppprojct/GraphClasses/GeneralGraph.cpp > CMakeFiles/cppprojct.dir/GraphClasses/GeneralGraph.cpp.i

CMakeFiles/cppprojct.dir/GraphClasses/GeneralGraph.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cppprojct.dir/GraphClasses/GeneralGraph.cpp.s"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/m4zz31/cppprojct/GraphClasses/GeneralGraph.cpp -o CMakeFiles/cppprojct.dir/GraphClasses/GeneralGraph.cpp.s

CMakeFiles/cppprojct.dir/GraphClasses/CliqueGraph.cpp.o: CMakeFiles/cppprojct.dir/flags.make
CMakeFiles/cppprojct.dir/GraphClasses/CliqueGraph.cpp.o: ../GraphClasses/CliqueGraph.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/m4zz31/cppprojct/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/cppprojct.dir/GraphClasses/CliqueGraph.cpp.o"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cppprojct.dir/GraphClasses/CliqueGraph.cpp.o -c /home/m4zz31/cppprojct/GraphClasses/CliqueGraph.cpp

CMakeFiles/cppprojct.dir/GraphClasses/CliqueGraph.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cppprojct.dir/GraphClasses/CliqueGraph.cpp.i"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/m4zz31/cppprojct/GraphClasses/CliqueGraph.cpp > CMakeFiles/cppprojct.dir/GraphClasses/CliqueGraph.cpp.i

CMakeFiles/cppprojct.dir/GraphClasses/CliqueGraph.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cppprojct.dir/GraphClasses/CliqueGraph.cpp.s"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/m4zz31/cppprojct/GraphClasses/CliqueGraph.cpp -o CMakeFiles/cppprojct.dir/GraphClasses/CliqueGraph.cpp.s

CMakeFiles/cppprojct.dir/Simulation/Simulation.cpp.o: CMakeFiles/cppprojct.dir/flags.make
CMakeFiles/cppprojct.dir/Simulation/Simulation.cpp.o: ../Simulation/Simulation.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/m4zz31/cppprojct/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/cppprojct.dir/Simulation/Simulation.cpp.o"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cppprojct.dir/Simulation/Simulation.cpp.o -c /home/m4zz31/cppprojct/Simulation/Simulation.cpp

CMakeFiles/cppprojct.dir/Simulation/Simulation.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cppprojct.dir/Simulation/Simulation.cpp.i"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/m4zz31/cppprojct/Simulation/Simulation.cpp > CMakeFiles/cppprojct.dir/Simulation/Simulation.cpp.i

CMakeFiles/cppprojct.dir/Simulation/Simulation.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cppprojct.dir/Simulation/Simulation.cpp.s"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/m4zz31/cppprojct/Simulation/Simulation.cpp -o CMakeFiles/cppprojct.dir/Simulation/Simulation.cpp.s

CMakeFiles/cppprojct.dir/Tests/graph-test-init.cpp.o: CMakeFiles/cppprojct.dir/flags.make
CMakeFiles/cppprojct.dir/Tests/graph-test-init.cpp.o: ../Tests/graph-test-init.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/m4zz31/cppprojct/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/cppprojct.dir/Tests/graph-test-init.cpp.o"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cppprojct.dir/Tests/graph-test-init.cpp.o -c /home/m4zz31/cppprojct/Tests/graph-test-init.cpp

CMakeFiles/cppprojct.dir/Tests/graph-test-init.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cppprojct.dir/Tests/graph-test-init.cpp.i"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/m4zz31/cppprojct/Tests/graph-test-init.cpp > CMakeFiles/cppprojct.dir/Tests/graph-test-init.cpp.i

CMakeFiles/cppprojct.dir/Tests/graph-test-init.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cppprojct.dir/Tests/graph-test-init.cpp.s"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/m4zz31/cppprojct/Tests/graph-test-init.cpp -o CMakeFiles/cppprojct.dir/Tests/graph-test-init.cpp.s

CMakeFiles/cppprojct.dir/Utils/reproductibility.cpp.o: CMakeFiles/cppprojct.dir/flags.make
CMakeFiles/cppprojct.dir/Utils/reproductibility.cpp.o: ../Utils/reproductibility.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/m4zz31/cppprojct/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/cppprojct.dir/Utils/reproductibility.cpp.o"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cppprojct.dir/Utils/reproductibility.cpp.o -c /home/m4zz31/cppprojct/Utils/reproductibility.cpp

CMakeFiles/cppprojct.dir/Utils/reproductibility.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cppprojct.dir/Utils/reproductibility.cpp.i"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/m4zz31/cppprojct/Utils/reproductibility.cpp > CMakeFiles/cppprojct.dir/Utils/reproductibility.cpp.i

CMakeFiles/cppprojct.dir/Utils/reproductibility.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cppprojct.dir/Utils/reproductibility.cpp.s"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/m4zz31/cppprojct/Utils/reproductibility.cpp -o CMakeFiles/cppprojct.dir/Utils/reproductibility.cpp.s

CMakeFiles/cppprojct.dir/Utils/adequate_synchronization.cpp.o: CMakeFiles/cppprojct.dir/flags.make
CMakeFiles/cppprojct.dir/Utils/adequate_synchronization.cpp.o: ../Utils/adequate_synchronization.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/m4zz31/cppprojct/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/cppprojct.dir/Utils/adequate_synchronization.cpp.o"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cppprojct.dir/Utils/adequate_synchronization.cpp.o -c /home/m4zz31/cppprojct/Utils/adequate_synchronization.cpp

CMakeFiles/cppprojct.dir/Utils/adequate_synchronization.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cppprojct.dir/Utils/adequate_synchronization.cpp.i"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/m4zz31/cppprojct/Utils/adequate_synchronization.cpp > CMakeFiles/cppprojct.dir/Utils/adequate_synchronization.cpp.i

CMakeFiles/cppprojct.dir/Utils/adequate_synchronization.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cppprojct.dir/Utils/adequate_synchronization.cpp.s"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/m4zz31/cppprojct/Utils/adequate_synchronization.cpp -o CMakeFiles/cppprojct.dir/Utils/adequate_synchronization.cpp.s

CMakeFiles/cppprojct.dir/GraphClasses/RingGraph.cpp.o: CMakeFiles/cppprojct.dir/flags.make
CMakeFiles/cppprojct.dir/GraphClasses/RingGraph.cpp.o: ../GraphClasses/RingGraph.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/m4zz31/cppprojct/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object CMakeFiles/cppprojct.dir/GraphClasses/RingGraph.cpp.o"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cppprojct.dir/GraphClasses/RingGraph.cpp.o -c /home/m4zz31/cppprojct/GraphClasses/RingGraph.cpp

CMakeFiles/cppprojct.dir/GraphClasses/RingGraph.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cppprojct.dir/GraphClasses/RingGraph.cpp.i"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/m4zz31/cppprojct/GraphClasses/RingGraph.cpp > CMakeFiles/cppprojct.dir/GraphClasses/RingGraph.cpp.i

CMakeFiles/cppprojct.dir/GraphClasses/RingGraph.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cppprojct.dir/GraphClasses/RingGraph.cpp.s"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/m4zz31/cppprojct/GraphClasses/RingGraph.cpp -o CMakeFiles/cppprojct.dir/GraphClasses/RingGraph.cpp.s

CMakeFiles/cppprojct.dir/Tests/graph-test-singlestep-evolution.cpp.o: CMakeFiles/cppprojct.dir/flags.make
CMakeFiles/cppprojct.dir/Tests/graph-test-singlestep-evolution.cpp.o: ../Tests/graph-test-singlestep-evolution.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/m4zz31/cppprojct/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object CMakeFiles/cppprojct.dir/Tests/graph-test-singlestep-evolution.cpp.o"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cppprojct.dir/Tests/graph-test-singlestep-evolution.cpp.o -c /home/m4zz31/cppprojct/Tests/graph-test-singlestep-evolution.cpp

CMakeFiles/cppprojct.dir/Tests/graph-test-singlestep-evolution.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cppprojct.dir/Tests/graph-test-singlestep-evolution.cpp.i"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/m4zz31/cppprojct/Tests/graph-test-singlestep-evolution.cpp > CMakeFiles/cppprojct.dir/Tests/graph-test-singlestep-evolution.cpp.i

CMakeFiles/cppprojct.dir/Tests/graph-test-singlestep-evolution.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cppprojct.dir/Tests/graph-test-singlestep-evolution.cpp.s"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/m4zz31/cppprojct/Tests/graph-test-singlestep-evolution.cpp -o CMakeFiles/cppprojct.dir/Tests/graph-test-singlestep-evolution.cpp.s

CMakeFiles/cppprojct.dir/Solvers/GeneralSolvers.cpp.o: CMakeFiles/cppprojct.dir/flags.make
CMakeFiles/cppprojct.dir/Solvers/GeneralSolvers.cpp.o: ../Solvers/GeneralSolvers.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/m4zz31/cppprojct/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building CXX object CMakeFiles/cppprojct.dir/Solvers/GeneralSolvers.cpp.o"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cppprojct.dir/Solvers/GeneralSolvers.cpp.o -c /home/m4zz31/cppprojct/Solvers/GeneralSolvers.cpp

CMakeFiles/cppprojct.dir/Solvers/GeneralSolvers.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cppprojct.dir/Solvers/GeneralSolvers.cpp.i"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/m4zz31/cppprojct/Solvers/GeneralSolvers.cpp > CMakeFiles/cppprojct.dir/Solvers/GeneralSolvers.cpp.i

CMakeFiles/cppprojct.dir/Solvers/GeneralSolvers.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cppprojct.dir/Solvers/GeneralSolvers.cpp.s"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/m4zz31/cppprojct/Solvers/GeneralSolvers.cpp -o CMakeFiles/cppprojct.dir/Solvers/GeneralSolvers.cpp.s

CMakeFiles/cppprojct.dir/Utils/memory_management.cpp.o: CMakeFiles/cppprojct.dir/flags.make
CMakeFiles/cppprojct.dir/Utils/memory_management.cpp.o: ../Utils/memory_management.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/m4zz31/cppprojct/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Building CXX object CMakeFiles/cppprojct.dir/Utils/memory_management.cpp.o"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cppprojct.dir/Utils/memory_management.cpp.o -c /home/m4zz31/cppprojct/Utils/memory_management.cpp

CMakeFiles/cppprojct.dir/Utils/memory_management.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cppprojct.dir/Utils/memory_management.cpp.i"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/m4zz31/cppprojct/Utils/memory_management.cpp > CMakeFiles/cppprojct.dir/Utils/memory_management.cpp.i

CMakeFiles/cppprojct.dir/Utils/memory_management.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cppprojct.dir/Utils/memory_management.cpp.s"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/m4zz31/cppprojct/Utils/memory_management.cpp -o CMakeFiles/cppprojct.dir/Utils/memory_management.cpp.s

CMakeFiles/cppprojct.dir/Tests/solvers-test-init.cpp.o: CMakeFiles/cppprojct.dir/flags.make
CMakeFiles/cppprojct.dir/Tests/solvers-test-init.cpp.o: ../Tests/solvers-test-init.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/m4zz31/cppprojct/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_15) "Building CXX object CMakeFiles/cppprojct.dir/Tests/solvers-test-init.cpp.o"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cppprojct.dir/Tests/solvers-test-init.cpp.o -c /home/m4zz31/cppprojct/Tests/solvers-test-init.cpp

CMakeFiles/cppprojct.dir/Tests/solvers-test-init.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cppprojct.dir/Tests/solvers-test-init.cpp.i"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/m4zz31/cppprojct/Tests/solvers-test-init.cpp > CMakeFiles/cppprojct.dir/Tests/solvers-test-init.cpp.i

CMakeFiles/cppprojct.dir/Tests/solvers-test-init.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cppprojct.dir/Tests/solvers-test-init.cpp.s"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/m4zz31/cppprojct/Tests/solvers-test-init.cpp -o CMakeFiles/cppprojct.dir/Tests/solvers-test-init.cpp.s

CMakeFiles/cppprojct.dir/DifferentialEquations/GeneralDifferentialEquation.cpp.o: CMakeFiles/cppprojct.dir/flags.make
CMakeFiles/cppprojct.dir/DifferentialEquations/GeneralDifferentialEquation.cpp.o: ../DifferentialEquations/GeneralDifferentialEquation.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/m4zz31/cppprojct/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_16) "Building CXX object CMakeFiles/cppprojct.dir/DifferentialEquations/GeneralDifferentialEquation.cpp.o"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cppprojct.dir/DifferentialEquations/GeneralDifferentialEquation.cpp.o -c /home/m4zz31/cppprojct/DifferentialEquations/GeneralDifferentialEquation.cpp

CMakeFiles/cppprojct.dir/DifferentialEquations/GeneralDifferentialEquation.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cppprojct.dir/DifferentialEquations/GeneralDifferentialEquation.cpp.i"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/m4zz31/cppprojct/DifferentialEquations/GeneralDifferentialEquation.cpp > CMakeFiles/cppprojct.dir/DifferentialEquations/GeneralDifferentialEquation.cpp.i

CMakeFiles/cppprojct.dir/DifferentialEquations/GeneralDifferentialEquation.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cppprojct.dir/DifferentialEquations/GeneralDifferentialEquation.cpp.s"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/m4zz31/cppprojct/DifferentialEquations/GeneralDifferentialEquation.cpp -o CMakeFiles/cppprojct.dir/DifferentialEquations/GeneralDifferentialEquation.cpp.s

CMakeFiles/cppprojct.dir/Solvers/EulerSolver.cpp.o: CMakeFiles/cppprojct.dir/flags.make
CMakeFiles/cppprojct.dir/Solvers/EulerSolver.cpp.o: ../Solvers/EulerSolver.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/m4zz31/cppprojct/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_17) "Building CXX object CMakeFiles/cppprojct.dir/Solvers/EulerSolver.cpp.o"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cppprojct.dir/Solvers/EulerSolver.cpp.o -c /home/m4zz31/cppprojct/Solvers/EulerSolver.cpp

CMakeFiles/cppprojct.dir/Solvers/EulerSolver.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cppprojct.dir/Solvers/EulerSolver.cpp.i"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/m4zz31/cppprojct/Solvers/EulerSolver.cpp > CMakeFiles/cppprojct.dir/Solvers/EulerSolver.cpp.i

CMakeFiles/cppprojct.dir/Solvers/EulerSolver.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cppprojct.dir/Solvers/EulerSolver.cpp.s"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/m4zz31/cppprojct/Solvers/EulerSolver.cpp -o CMakeFiles/cppprojct.dir/Solvers/EulerSolver.cpp.s

CMakeFiles/cppprojct.dir/Utils/differential_equations_aux.cpp.o: CMakeFiles/cppprojct.dir/flags.make
CMakeFiles/cppprojct.dir/Utils/differential_equations_aux.cpp.o: ../Utils/differential_equations_aux.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/m4zz31/cppprojct/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_18) "Building CXX object CMakeFiles/cppprojct.dir/Utils/differential_equations_aux.cpp.o"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cppprojct.dir/Utils/differential_equations_aux.cpp.o -c /home/m4zz31/cppprojct/Utils/differential_equations_aux.cpp

CMakeFiles/cppprojct.dir/Utils/differential_equations_aux.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cppprojct.dir/Utils/differential_equations_aux.cpp.i"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/m4zz31/cppprojct/Utils/differential_equations_aux.cpp > CMakeFiles/cppprojct.dir/Utils/differential_equations_aux.cpp.i

CMakeFiles/cppprojct.dir/Utils/differential_equations_aux.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cppprojct.dir/Utils/differential_equations_aux.cpp.s"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/m4zz31/cppprojct/Utils/differential_equations_aux.cpp -o CMakeFiles/cppprojct.dir/Utils/differential_equations_aux.cpp.s

CMakeFiles/cppprojct.dir/DifferentialEquations/NoiselessKuramoto.cpp.o: CMakeFiles/cppprojct.dir/flags.make
CMakeFiles/cppprojct.dir/DifferentialEquations/NoiselessKuramoto.cpp.o: ../DifferentialEquations/NoiselessKuramoto.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/m4zz31/cppprojct/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_19) "Building CXX object CMakeFiles/cppprojct.dir/DifferentialEquations/NoiselessKuramoto.cpp.o"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cppprojct.dir/DifferentialEquations/NoiselessKuramoto.cpp.o -c /home/m4zz31/cppprojct/DifferentialEquations/NoiselessKuramoto.cpp

CMakeFiles/cppprojct.dir/DifferentialEquations/NoiselessKuramoto.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cppprojct.dir/DifferentialEquations/NoiselessKuramoto.cpp.i"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/m4zz31/cppprojct/DifferentialEquations/NoiselessKuramoto.cpp > CMakeFiles/cppprojct.dir/DifferentialEquations/NoiselessKuramoto.cpp.i

CMakeFiles/cppprojct.dir/DifferentialEquations/NoiselessKuramoto.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cppprojct.dir/DifferentialEquations/NoiselessKuramoto.cpp.s"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/m4zz31/cppprojct/DifferentialEquations/NoiselessKuramoto.cpp -o CMakeFiles/cppprojct.dir/DifferentialEquations/NoiselessKuramoto.cpp.s

# Object files for target cppprojct
cppprojct_OBJECTS = \
"CMakeFiles/cppprojct.dir/main.cpp.o" \
"CMakeFiles/cppprojct.dir/Utils/timers.cpp.o" \
"CMakeFiles/cppprojct.dir/GraphClasses/ErdosRenyiGraph.cpp.o" \
"CMakeFiles/cppprojct.dir/Utils/error.cpp.o" \
"CMakeFiles/cppprojct.dir/GraphClasses/GeneralGraph.cpp.o" \
"CMakeFiles/cppprojct.dir/GraphClasses/CliqueGraph.cpp.o" \
"CMakeFiles/cppprojct.dir/Simulation/Simulation.cpp.o" \
"CMakeFiles/cppprojct.dir/Tests/graph-test-init.cpp.o" \
"CMakeFiles/cppprojct.dir/Utils/reproductibility.cpp.o" \
"CMakeFiles/cppprojct.dir/Utils/adequate_synchronization.cpp.o" \
"CMakeFiles/cppprojct.dir/GraphClasses/RingGraph.cpp.o" \
"CMakeFiles/cppprojct.dir/Tests/graph-test-singlestep-evolution.cpp.o" \
"CMakeFiles/cppprojct.dir/Solvers/GeneralSolvers.cpp.o" \
"CMakeFiles/cppprojct.dir/Utils/memory_management.cpp.o" \
"CMakeFiles/cppprojct.dir/Tests/solvers-test-init.cpp.o" \
"CMakeFiles/cppprojct.dir/DifferentialEquations/GeneralDifferentialEquation.cpp.o" \
"CMakeFiles/cppprojct.dir/Solvers/EulerSolver.cpp.o" \
"CMakeFiles/cppprojct.dir/Utils/differential_equations_aux.cpp.o" \
"CMakeFiles/cppprojct.dir/DifferentialEquations/NoiselessKuramoto.cpp.o"

# External object files for target cppprojct
cppprojct_EXTERNAL_OBJECTS =

cppprojct: CMakeFiles/cppprojct.dir/main.cpp.o
cppprojct: CMakeFiles/cppprojct.dir/Utils/timers.cpp.o
cppprojct: CMakeFiles/cppprojct.dir/GraphClasses/ErdosRenyiGraph.cpp.o
cppprojct: CMakeFiles/cppprojct.dir/Utils/error.cpp.o
cppprojct: CMakeFiles/cppprojct.dir/GraphClasses/GeneralGraph.cpp.o
cppprojct: CMakeFiles/cppprojct.dir/GraphClasses/CliqueGraph.cpp.o
cppprojct: CMakeFiles/cppprojct.dir/Simulation/Simulation.cpp.o
cppprojct: CMakeFiles/cppprojct.dir/Tests/graph-test-init.cpp.o
cppprojct: CMakeFiles/cppprojct.dir/Utils/reproductibility.cpp.o
cppprojct: CMakeFiles/cppprojct.dir/Utils/adequate_synchronization.cpp.o
cppprojct: CMakeFiles/cppprojct.dir/GraphClasses/RingGraph.cpp.o
cppprojct: CMakeFiles/cppprojct.dir/Tests/graph-test-singlestep-evolution.cpp.o
cppprojct: CMakeFiles/cppprojct.dir/Solvers/GeneralSolvers.cpp.o
cppprojct: CMakeFiles/cppprojct.dir/Utils/memory_management.cpp.o
cppprojct: CMakeFiles/cppprojct.dir/Tests/solvers-test-init.cpp.o
cppprojct: CMakeFiles/cppprojct.dir/DifferentialEquations/GeneralDifferentialEquation.cpp.o
cppprojct: CMakeFiles/cppprojct.dir/Solvers/EulerSolver.cpp.o
cppprojct: CMakeFiles/cppprojct.dir/Utils/differential_equations_aux.cpp.o
cppprojct: CMakeFiles/cppprojct.dir/DifferentialEquations/NoiselessKuramoto.cpp.o
cppprojct: CMakeFiles/cppprojct.dir/build.make
cppprojct: /usr/lib/x86_64-linux-gnu/libboost_mpi.so
cppprojct: /usr/lib/x86_64-linux-gnu/libboost_graph_parallel.so
cppprojct: /usr/lib/x86_64-linux-gnu/libboost_system.so
cppprojct: /usr/lib/x86_64-linux-gnu/libboost_serialization.so
cppprojct: /usr/lib/gcc/x86_64-linux-gnu/7/libgomp.so
cppprojct: /usr/lib/x86_64-linux-gnu/libpthread.so
cppprojct: CMakeFiles/cppprojct.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/m4zz31/cppprojct/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_20) "Linking CXX executable cppprojct"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/cppprojct.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/cppprojct.dir/build: cppprojct

.PHONY : CMakeFiles/cppprojct.dir/build

CMakeFiles/cppprojct.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/cppprojct.dir/cmake_clean.cmake
.PHONY : CMakeFiles/cppprojct.dir/clean

CMakeFiles/cppprojct.dir/depend:
	cd /home/m4zz31/cppprojct/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/m4zz31/cppprojct /home/m4zz31/cppprojct /home/m4zz31/cppprojct/cmake-build-debug /home/m4zz31/cppprojct/cmake-build-debug /home/m4zz31/cppprojct/cmake-build-debug/CMakeFiles/cppprojct.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/cppprojct.dir/depend

