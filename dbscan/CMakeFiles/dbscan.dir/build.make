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
CMAKE_SOURCE_DIR = /home/ncslabtb/SLAM_SYSU/SLAM

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ncslabtb/SLAM_SYSU/SLAM

# Include any dependencies generated for this target.
include dbscan/CMakeFiles/dbscan.dir/depend.make

# Include the progress variables for this target.
include dbscan/CMakeFiles/dbscan.dir/progress.make

# Include the compile flags for this target's objects.
include dbscan/CMakeFiles/dbscan.dir/flags.make

dbscan/CMakeFiles/dbscan.dir/dbscan.cpp.o: dbscan/CMakeFiles/dbscan.dir/flags.make
dbscan/CMakeFiles/dbscan.dir/dbscan.cpp.o: dbscan/dbscan.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/ncslabtb/SLAM_SYSU/SLAM/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object dbscan/CMakeFiles/dbscan.dir/dbscan.cpp.o"
	cd /home/ncslabtb/SLAM_SYSU/SLAM/dbscan && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/dbscan.dir/dbscan.cpp.o -c /home/ncslabtb/SLAM_SYSU/SLAM/dbscan/dbscan.cpp

dbscan/CMakeFiles/dbscan.dir/dbscan.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dbscan.dir/dbscan.cpp.i"
	cd /home/ncslabtb/SLAM_SYSU/SLAM/dbscan && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/ncslabtb/SLAM_SYSU/SLAM/dbscan/dbscan.cpp > CMakeFiles/dbscan.dir/dbscan.cpp.i

dbscan/CMakeFiles/dbscan.dir/dbscan.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dbscan.dir/dbscan.cpp.s"
	cd /home/ncslabtb/SLAM_SYSU/SLAM/dbscan && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/ncslabtb/SLAM_SYSU/SLAM/dbscan/dbscan.cpp -o CMakeFiles/dbscan.dir/dbscan.cpp.s

dbscan/CMakeFiles/dbscan.dir/dbscan.cpp.o.requires:
.PHONY : dbscan/CMakeFiles/dbscan.dir/dbscan.cpp.o.requires

dbscan/CMakeFiles/dbscan.dir/dbscan.cpp.o.provides: dbscan/CMakeFiles/dbscan.dir/dbscan.cpp.o.requires
	$(MAKE) -f dbscan/CMakeFiles/dbscan.dir/build.make dbscan/CMakeFiles/dbscan.dir/dbscan.cpp.o.provides.build
.PHONY : dbscan/CMakeFiles/dbscan.dir/dbscan.cpp.o.provides

dbscan/CMakeFiles/dbscan.dir/dbscan.cpp.o.provides.build: dbscan/CMakeFiles/dbscan.dir/dbscan.cpp.o

# Object files for target dbscan
dbscan_OBJECTS = \
"CMakeFiles/dbscan.dir/dbscan.cpp.o"

# External object files for target dbscan
dbscan_EXTERNAL_OBJECTS =

lib/libdbscan.a: dbscan/CMakeFiles/dbscan.dir/dbscan.cpp.o
lib/libdbscan.a: dbscan/CMakeFiles/dbscan.dir/build.make
lib/libdbscan.a: dbscan/CMakeFiles/dbscan.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library ../lib/libdbscan.a"
	cd /home/ncslabtb/SLAM_SYSU/SLAM/dbscan && $(CMAKE_COMMAND) -P CMakeFiles/dbscan.dir/cmake_clean_target.cmake
	cd /home/ncslabtb/SLAM_SYSU/SLAM/dbscan && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/dbscan.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
dbscan/CMakeFiles/dbscan.dir/build: lib/libdbscan.a
.PHONY : dbscan/CMakeFiles/dbscan.dir/build

dbscan/CMakeFiles/dbscan.dir/requires: dbscan/CMakeFiles/dbscan.dir/dbscan.cpp.o.requires
.PHONY : dbscan/CMakeFiles/dbscan.dir/requires

dbscan/CMakeFiles/dbscan.dir/clean:
	cd /home/ncslabtb/SLAM_SYSU/SLAM/dbscan && $(CMAKE_COMMAND) -P CMakeFiles/dbscan.dir/cmake_clean.cmake
.PHONY : dbscan/CMakeFiles/dbscan.dir/clean

dbscan/CMakeFiles/dbscan.dir/depend:
	cd /home/ncslabtb/SLAM_SYSU/SLAM && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ncslabtb/SLAM_SYSU/SLAM /home/ncslabtb/SLAM_SYSU/SLAM/dbscan /home/ncslabtb/SLAM_SYSU/SLAM /home/ncslabtb/SLAM_SYSU/SLAM/dbscan /home/ncslabtb/SLAM_SYSU/SLAM/dbscan/CMakeFiles/dbscan.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : dbscan/CMakeFiles/dbscan.dir/depend

