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

CMakeFiles/cppprojct.dir/utils/timers.cpp.o: CMakeFiles/cppprojct.dir/flags.make
CMakeFiles/cppprojct.dir/utils/timers.cpp.o: ../utils/timers.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/m4zz31/cppprojct/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/cppprojct.dir/utils/timers.cpp.o"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cppprojct.dir/utils/timers.cpp.o -c /home/m4zz31/cppprojct/utils/timers.cpp

CMakeFiles/cppprojct.dir/utils/timers.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cppprojct.dir/utils/timers.cpp.i"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/m4zz31/cppprojct/utils/timers.cpp > CMakeFiles/cppprojct.dir/utils/timers.cpp.i

CMakeFiles/cppprojct.dir/utils/timers.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cppprojct.dir/utils/timers.cpp.s"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/m4zz31/cppprojct/utils/timers.cpp -o CMakeFiles/cppprojct.dir/utils/timers.cpp.s

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

CMakeFiles/cppprojct.dir/utils/error.cpp.o: CMakeFiles/cppprojct.dir/flags.make
CMakeFiles/cppprojct.dir/utils/error.cpp.o: ../utils/error.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/m4zz31/cppprojct/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/cppprojct.dir/utils/error.cpp.o"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cppprojct.dir/utils/error.cpp.o -c /home/m4zz31/cppprojct/utils/error.cpp

CMakeFiles/cppprojct.dir/utils/error.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cppprojct.dir/utils/error.cpp.i"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/m4zz31/cppprojct/utils/error.cpp > CMakeFiles/cppprojct.dir/utils/error.cpp.i

CMakeFiles/cppprojct.dir/utils/error.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cppprojct.dir/utils/error.cpp.s"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/m4zz31/cppprojct/utils/error.cpp -o CMakeFiles/cppprojct.dir/utils/error.cpp.s

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

CMakeFiles/cppprojct.dir/tests/graph-test-init.cpp.o: CMakeFiles/cppprojct.dir/flags.make
CMakeFiles/cppprojct.dir/tests/graph-test-init.cpp.o: ../tests/graph-test-init.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/m4zz31/cppprojct/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/cppprojct.dir/tests/graph-test-init.cpp.o"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cppprojct.dir/tests/graph-test-init.cpp.o -c /home/m4zz31/cppprojct/tests/graph-test-init.cpp

CMakeFiles/cppprojct.dir/tests/graph-test-init.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cppprojct.dir/tests/graph-test-init.cpp.i"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/m4zz31/cppprojct/tests/graph-test-init.cpp > CMakeFiles/cppprojct.dir/tests/graph-test-init.cpp.i

CMakeFiles/cppprojct.dir/tests/graph-test-init.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cppprojct.dir/tests/graph-test-init.cpp.s"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/m4zz31/cppprojct/tests/graph-test-init.cpp -o CMakeFiles/cppprojct.dir/tests/graph-test-init.cpp.s

CMakeFiles/cppprojct.dir/utils/reproductibility.cpp.o: CMakeFiles/cppprojct.dir/flags.make
CMakeFiles/cppprojct.dir/utils/reproductibility.cpp.o: ../utils/reproductibility.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/m4zz31/cppprojct/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/cppprojct.dir/utils/reproductibility.cpp.o"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cppprojct.dir/utils/reproductibility.cpp.o -c /home/m4zz31/cppprojct/utils/reproductibility.cpp

CMakeFiles/cppprojct.dir/utils/reproductibility.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cppprojct.dir/utils/reproductibility.cpp.i"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/m4zz31/cppprojct/utils/reproductibility.cpp > CMakeFiles/cppprojct.dir/utils/reproductibility.cpp.i

CMakeFiles/cppprojct.dir/utils/reproductibility.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cppprojct.dir/utils/reproductibility.cpp.s"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/m4zz31/cppprojct/utils/reproductibility.cpp -o CMakeFiles/cppprojct.dir/utils/reproductibility.cpp.s

CMakeFiles/cppprojct.dir/utils/adequate_synchronization.cpp.o: CMakeFiles/cppprojct.dir/flags.make
CMakeFiles/cppprojct.dir/utils/adequate_synchronization.cpp.o: ../utils/adequate_synchronization.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/m4zz31/cppprojct/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/cppprojct.dir/utils/adequate_synchronization.cpp.o"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cppprojct.dir/utils/adequate_synchronization.cpp.o -c /home/m4zz31/cppprojct/utils/adequate_synchronization.cpp

CMakeFiles/cppprojct.dir/utils/adequate_synchronization.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cppprojct.dir/utils/adequate_synchronization.cpp.i"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/m4zz31/cppprojct/utils/adequate_synchronization.cpp > CMakeFiles/cppprojct.dir/utils/adequate_synchronization.cpp.i

CMakeFiles/cppprojct.dir/utils/adequate_synchronization.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cppprojct.dir/utils/adequate_synchronization.cpp.s"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/m4zz31/cppprojct/utils/adequate_synchronization.cpp -o CMakeFiles/cppprojct.dir/utils/adequate_synchronization.cpp.s

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

CMakeFiles/cppprojct.dir/tests/graph-test-singlestep-evolution.cpp.o: CMakeFiles/cppprojct.dir/flags.make
CMakeFiles/cppprojct.dir/tests/graph-test-singlestep-evolution.cpp.o: ../tests/graph-test-singlestep-evolution.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/m4zz31/cppprojct/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object CMakeFiles/cppprojct.dir/tests/graph-test-singlestep-evolution.cpp.o"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cppprojct.dir/tests/graph-test-singlestep-evolution.cpp.o -c /home/m4zz31/cppprojct/tests/graph-test-singlestep-evolution.cpp

CMakeFiles/cppprojct.dir/tests/graph-test-singlestep-evolution.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cppprojct.dir/tests/graph-test-singlestep-evolution.cpp.i"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/m4zz31/cppprojct/tests/graph-test-singlestep-evolution.cpp > CMakeFiles/cppprojct.dir/tests/graph-test-singlestep-evolution.cpp.i

CMakeFiles/cppprojct.dir/tests/graph-test-singlestep-evolution.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cppprojct.dir/tests/graph-test-singlestep-evolution.cpp.s"
	mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/m4zz31/cppprojct/tests/graph-test-singlestep-evolution.cpp -o CMakeFiles/cppprojct.dir/tests/graph-test-singlestep-evolution.cpp.s

# Object files for target cppprojct
cppprojct_OBJECTS = \
"CMakeFiles/cppprojct.dir/main.cpp.o" \
"CMakeFiles/cppprojct.dir/utils/timers.cpp.o" \
"CMakeFiles/cppprojct.dir/GraphClasses/ErdosRenyiGraph.cpp.o" \
"CMakeFiles/cppprojct.dir/utils/error.cpp.o" \
"CMakeFiles/cppprojct.dir/GraphClasses/GeneralGraph.cpp.o" \
"CMakeFiles/cppprojct.dir/GraphClasses/CliqueGraph.cpp.o" \
"CMakeFiles/cppprojct.dir/Simulation/Simulation.cpp.o" \
"CMakeFiles/cppprojct.dir/tests/graph-test-init.cpp.o" \
"CMakeFiles/cppprojct.dir/utils/reproductibility.cpp.o" \
"CMakeFiles/cppprojct.dir/utils/adequate_synchronization.cpp.o" \
"CMakeFiles/cppprojct.dir/GraphClasses/RingGraph.cpp.o" \
"CMakeFiles/cppprojct.dir/tests/graph-test-singlestep-evolution.cpp.o"

# External object files for target cppprojct
cppprojct_EXTERNAL_OBJECTS =

cppprojct: CMakeFiles/cppprojct.dir/main.cpp.o
cppprojct: CMakeFiles/cppprojct.dir/utils/timers.cpp.o
cppprojct: CMakeFiles/cppprojct.dir/GraphClasses/ErdosRenyiGraph.cpp.o
cppprojct: CMakeFiles/cppprojct.dir/utils/error.cpp.o
cppprojct: CMakeFiles/cppprojct.dir/GraphClasses/GeneralGraph.cpp.o
cppprojct: CMakeFiles/cppprojct.dir/GraphClasses/CliqueGraph.cpp.o
cppprojct: CMakeFiles/cppprojct.dir/Simulation/Simulation.cpp.o
cppprojct: CMakeFiles/cppprojct.dir/tests/graph-test-init.cpp.o
cppprojct: CMakeFiles/cppprojct.dir/utils/reproductibility.cpp.o
cppprojct: CMakeFiles/cppprojct.dir/utils/adequate_synchronization.cpp.o
cppprojct: CMakeFiles/cppprojct.dir/GraphClasses/RingGraph.cpp.o
cppprojct: CMakeFiles/cppprojct.dir/tests/graph-test-singlestep-evolution.cpp.o
cppprojct: CMakeFiles/cppprojct.dir/build.make
cppprojct: /usr/lib/x86_64-linux-gnu/libboost_mpi.so
cppprojct: /usr/lib/x86_64-linux-gnu/libboost_graph_parallel.so
cppprojct: /usr/lib/x86_64-linux-gnu/libboost_system.so
cppprojct: /usr/lib/x86_64-linux-gnu/libboost_serialization.so
cppprojct: /usr/lib/gcc/x86_64-linux-gnu/7/libgomp.so
cppprojct: /usr/lib/x86_64-linux-gnu/libpthread.so
cppprojct: CMakeFiles/cppprojct.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/m4zz31/cppprojct/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Linking CXX executable cppprojct"
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

