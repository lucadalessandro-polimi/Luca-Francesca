# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

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
CMAKE_COMMAND = /u/sw/toolchains/gcc-glibc/11.2.0/base/bin/cmake

# The command to remove a file.
RM = /u/sw/toolchains/gcc-glibc/11.2.0/base/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/lucadalessandro/Luca-Francesca/test

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/lucadalessandro/Luca-Francesca/test

# Include any dependencies generated for this target.
include CMakeFiles/fdapde_test.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/fdapde_test.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/fdapde_test.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/fdapde_test.dir/flags.make

CMakeFiles/fdapde_test.dir/main.cpp.o: CMakeFiles/fdapde_test.dir/flags.make
CMakeFiles/fdapde_test.dir/main.cpp.o: main.cpp
CMakeFiles/fdapde_test.dir/main.cpp.o: CMakeFiles/fdapde_test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lucadalessandro/Luca-Francesca/test/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/fdapde_test.dir/main.cpp.o"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/fdapde_test.dir/main.cpp.o -MF CMakeFiles/fdapde_test.dir/main.cpp.o.d -o CMakeFiles/fdapde_test.dir/main.cpp.o -c /home/lucadalessandro/Luca-Francesca/test/main.cpp

CMakeFiles/fdapde_test.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fdapde_test.dir/main.cpp.i"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/lucadalessandro/Luca-Francesca/test/main.cpp > CMakeFiles/fdapde_test.dir/main.cpp.i

CMakeFiles/fdapde_test.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fdapde_test.dir/main.cpp.s"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/lucadalessandro/Luca-Francesca/test/main.cpp -o CMakeFiles/fdapde_test.dir/main.cpp.s

# Object files for target fdapde_test
fdapde_test_OBJECTS = \
"CMakeFiles/fdapde_test.dir/main.cpp.o"

# External object files for target fdapde_test
fdapde_test_EXTERNAL_OBJECTS =

fdapde_test: CMakeFiles/fdapde_test.dir/main.cpp.o
fdapde_test: CMakeFiles/fdapde_test.dir/build.make
fdapde_test: lib/libgtest_main.a
fdapde_test: lib/libgtest.a
fdapde_test: CMakeFiles/fdapde_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/lucadalessandro/Luca-Francesca/test/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable fdapde_test"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/fdapde_test.dir/link.txt --verbose=$(VERBOSE)
	/u/sw/toolchains/gcc-glibc/11.2.0/base/bin/cmake -D TEST_TARGET=fdapde_test -D TEST_EXECUTABLE=/home/lucadalessandro/Luca-Francesca/test/fdapde_test -D TEST_EXECUTOR= -D TEST_WORKING_DIR=/home/lucadalessandro/Luca-Francesca/test -D TEST_EXTRA_ARGS= -D TEST_PROPERTIES= -D TEST_PREFIX= -D TEST_SUFFIX= -D NO_PRETTY_TYPES=FALSE -D NO_PRETTY_VALUES=FALSE -D TEST_LIST=fdapde_test_TESTS -D CTEST_FILE=/home/lucadalessandro/Luca-Francesca/test/fdapde_test[1]_tests.cmake -D TEST_DISCOVERY_TIMEOUT=5 -D TEST_XML_OUTPUT_DIR= -P /u/sw/toolchains/gcc-glibc/11.2.0/base/share/cmake-3.20/Modules/GoogleTestAddTests.cmake

# Rule to build all files generated by this target.
CMakeFiles/fdapde_test.dir/build: fdapde_test
.PHONY : CMakeFiles/fdapde_test.dir/build

CMakeFiles/fdapde_test.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/fdapde_test.dir/cmake_clean.cmake
.PHONY : CMakeFiles/fdapde_test.dir/clean

CMakeFiles/fdapde_test.dir/depend:
	cd /home/lucadalessandro/Luca-Francesca/test && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/lucadalessandro/Luca-Francesca/test /home/lucadalessandro/Luca-Francesca/test /home/lucadalessandro/Luca-Francesca/test /home/lucadalessandro/Luca-Francesca/test /home/lucadalessandro/Luca-Francesca/test/CMakeFiles/fdapde_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/fdapde_test.dir/depend

