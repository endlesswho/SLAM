# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/home/ncslabtb/SLAM_SYSU/slam-sysu (copy)"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/ncslabtb/SLAM_SYSU/slam-sysu (copy)"

# Include any dependencies generated for this target.
include CMakeFiles/slam.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/slam.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/slam.dir/flags.make

CMakeFiles/slam.dir/main.cpp.o: CMakeFiles/slam.dir/flags.make
CMakeFiles/slam.dir/main.cpp.o: main.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report "/home/ncslabtb/SLAM_SYSU/slam-sysu (copy)/CMakeFiles" $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/slam.dir/main.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/slam.dir/main.cpp.o -c "/home/ncslabtb/SLAM_SYSU/slam-sysu (copy)/main.cpp"

CMakeFiles/slam.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/slam.dir/main.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E "/home/ncslabtb/SLAM_SYSU/slam-sysu (copy)/main.cpp" > CMakeFiles/slam.dir/main.cpp.i

CMakeFiles/slam.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/slam.dir/main.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S "/home/ncslabtb/SLAM_SYSU/slam-sysu (copy)/main.cpp" -o CMakeFiles/slam.dir/main.cpp.s

CMakeFiles/slam.dir/main.cpp.o.requires:
.PHONY : CMakeFiles/slam.dir/main.cpp.o.requires

CMakeFiles/slam.dir/main.cpp.o.provides: CMakeFiles/slam.dir/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/slam.dir/build.make CMakeFiles/slam.dir/main.cpp.o.provides.build
.PHONY : CMakeFiles/slam.dir/main.cpp.o.provides

CMakeFiles/slam.dir/main.cpp.o.provides.build: CMakeFiles/slam.dir/main.cpp.o

CMakeFiles/slam.dir/icp.cpp.o: CMakeFiles/slam.dir/flags.make
CMakeFiles/slam.dir/icp.cpp.o: icp.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report "/home/ncslabtb/SLAM_SYSU/slam-sysu (copy)/CMakeFiles" $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/slam.dir/icp.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/slam.dir/icp.cpp.o -c "/home/ncslabtb/SLAM_SYSU/slam-sysu (copy)/icp.cpp"

CMakeFiles/slam.dir/icp.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/slam.dir/icp.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E "/home/ncslabtb/SLAM_SYSU/slam-sysu (copy)/icp.cpp" > CMakeFiles/slam.dir/icp.cpp.i

CMakeFiles/slam.dir/icp.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/slam.dir/icp.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S "/home/ncslabtb/SLAM_SYSU/slam-sysu (copy)/icp.cpp" -o CMakeFiles/slam.dir/icp.cpp.s

CMakeFiles/slam.dir/icp.cpp.o.requires:
.PHONY : CMakeFiles/slam.dir/icp.cpp.o.requires

CMakeFiles/slam.dir/icp.cpp.o.provides: CMakeFiles/slam.dir/icp.cpp.o.requires
	$(MAKE) -f CMakeFiles/slam.dir/build.make CMakeFiles/slam.dir/icp.cpp.o.provides.build
.PHONY : CMakeFiles/slam.dir/icp.cpp.o.provides

CMakeFiles/slam.dir/icp.cpp.o.provides.build: CMakeFiles/slam.dir/icp.cpp.o

# Object files for target slam
slam_OBJECTS = \
"CMakeFiles/slam.dir/main.cpp.o" \
"CMakeFiles/slam.dir/icp.cpp.o"

# External object files for target slam
slam_EXTERNAL_OBJECTS =

bin/slam: CMakeFiles/slam.dir/main.cpp.o
bin/slam: CMakeFiles/slam.dir/icp.cpp.o
bin/slam: CMakeFiles/slam.dir/build.make
bin/slam: /usr/local/lib/libpointmatcher.so
bin/slam: /usr/lib/x86_64-linux-gnu/libboost_thread.so
bin/slam: /usr/lib/x86_64-linux-gnu/libboost_filesystem.so
bin/slam: /usr/lib/x86_64-linux-gnu/libboost_system.so
bin/slam: /usr/lib/x86_64-linux-gnu/libboost_program_options.so
bin/slam: /usr/lib/x86_64-linux-gnu/libboost_date_time.so
bin/slam: /usr/lib/x86_64-linux-gnu/libboost_chrono.so
bin/slam: /usr/lib/x86_64-linux-gnu/libpthread.so
bin/slam: /home/ncslabtb/ethzasl_icp_mapping/libnabo/build/libnabo.a
bin/slam: lib/libdbscan.a
bin/slam: lib/libransac.a
bin/slam: CMakeFiles/slam.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable bin/slam"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/slam.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/slam.dir/build: bin/slam
.PHONY : CMakeFiles/slam.dir/build

CMakeFiles/slam.dir/requires: CMakeFiles/slam.dir/main.cpp.o.requires
CMakeFiles/slam.dir/requires: CMakeFiles/slam.dir/icp.cpp.o.requires
.PHONY : CMakeFiles/slam.dir/requires

CMakeFiles/slam.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/slam.dir/cmake_clean.cmake
.PHONY : CMakeFiles/slam.dir/clean

CMakeFiles/slam.dir/depend:
	cd "/home/ncslabtb/SLAM_SYSU/slam-sysu (copy)" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/ncslabtb/SLAM_SYSU/slam-sysu (copy)" "/home/ncslabtb/SLAM_SYSU/slam-sysu (copy)" "/home/ncslabtb/SLAM_SYSU/slam-sysu (copy)" "/home/ncslabtb/SLAM_SYSU/slam-sysu (copy)" "/home/ncslabtb/SLAM_SYSU/slam-sysu (copy)/CMakeFiles/slam.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/slam.dir/depend
