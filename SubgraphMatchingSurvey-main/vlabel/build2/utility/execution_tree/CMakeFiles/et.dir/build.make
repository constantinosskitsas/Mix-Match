# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_SOURCE_DIR = /home/konstantinos/SubgraphMatchingSurvey-main/SubgraphMatchingSurvey-main/vlabel

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/konstantinos/SubgraphMatchingSurvey-main/SubgraphMatchingSurvey-main/vlabel/build2

# Include any dependencies generated for this target.
include utility/execution_tree/CMakeFiles/et.dir/depend.make

# Include the progress variables for this target.
include utility/execution_tree/CMakeFiles/et.dir/progress.make

# Include the compile flags for this target's objects.
include utility/execution_tree/CMakeFiles/et.dir/flags.make

utility/execution_tree/CMakeFiles/et.dir/execution_tree_node.cpp.o: utility/execution_tree/CMakeFiles/et.dir/flags.make
utility/execution_tree/CMakeFiles/et.dir/execution_tree_node.cpp.o: ../utility/execution_tree/execution_tree_node.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/konstantinos/SubgraphMatchingSurvey-main/SubgraphMatchingSurvey-main/vlabel/build2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object utility/execution_tree/CMakeFiles/et.dir/execution_tree_node.cpp.o"
	cd /home/konstantinos/SubgraphMatchingSurvey-main/SubgraphMatchingSurvey-main/vlabel/build2/utility/execution_tree && /home/konstantinos/miniconda3/bin/x86_64-conda-linux-gnu-c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/et.dir/execution_tree_node.cpp.o -c /home/konstantinos/SubgraphMatchingSurvey-main/SubgraphMatchingSurvey-main/vlabel/utility/execution_tree/execution_tree_node.cpp

utility/execution_tree/CMakeFiles/et.dir/execution_tree_node.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/et.dir/execution_tree_node.cpp.i"
	cd /home/konstantinos/SubgraphMatchingSurvey-main/SubgraphMatchingSurvey-main/vlabel/build2/utility/execution_tree && /home/konstantinos/miniconda3/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/konstantinos/SubgraphMatchingSurvey-main/SubgraphMatchingSurvey-main/vlabel/utility/execution_tree/execution_tree_node.cpp > CMakeFiles/et.dir/execution_tree_node.cpp.i

utility/execution_tree/CMakeFiles/et.dir/execution_tree_node.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/et.dir/execution_tree_node.cpp.s"
	cd /home/konstantinos/SubgraphMatchingSurvey-main/SubgraphMatchingSurvey-main/vlabel/build2/utility/execution_tree && /home/konstantinos/miniconda3/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/konstantinos/SubgraphMatchingSurvey-main/SubgraphMatchingSurvey-main/vlabel/utility/execution_tree/execution_tree_node.cpp -o CMakeFiles/et.dir/execution_tree_node.cpp.s

utility/execution_tree/CMakeFiles/et.dir/execution_tree_node.cpp.o.requires:

.PHONY : utility/execution_tree/CMakeFiles/et.dir/execution_tree_node.cpp.o.requires

utility/execution_tree/CMakeFiles/et.dir/execution_tree_node.cpp.o.provides: utility/execution_tree/CMakeFiles/et.dir/execution_tree_node.cpp.o.requires
	$(MAKE) -f utility/execution_tree/CMakeFiles/et.dir/build.make utility/execution_tree/CMakeFiles/et.dir/execution_tree_node.cpp.o.provides.build
.PHONY : utility/execution_tree/CMakeFiles/et.dir/execution_tree_node.cpp.o.provides

utility/execution_tree/CMakeFiles/et.dir/execution_tree_node.cpp.o.provides.build: utility/execution_tree/CMakeFiles/et.dir/execution_tree_node.cpp.o


utility/execution_tree/CMakeFiles/et.dir/execution_tree.cpp.o: utility/execution_tree/CMakeFiles/et.dir/flags.make
utility/execution_tree/CMakeFiles/et.dir/execution_tree.cpp.o: ../utility/execution_tree/execution_tree.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/konstantinos/SubgraphMatchingSurvey-main/SubgraphMatchingSurvey-main/vlabel/build2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object utility/execution_tree/CMakeFiles/et.dir/execution_tree.cpp.o"
	cd /home/konstantinos/SubgraphMatchingSurvey-main/SubgraphMatchingSurvey-main/vlabel/build2/utility/execution_tree && /home/konstantinos/miniconda3/bin/x86_64-conda-linux-gnu-c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/et.dir/execution_tree.cpp.o -c /home/konstantinos/SubgraphMatchingSurvey-main/SubgraphMatchingSurvey-main/vlabel/utility/execution_tree/execution_tree.cpp

utility/execution_tree/CMakeFiles/et.dir/execution_tree.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/et.dir/execution_tree.cpp.i"
	cd /home/konstantinos/SubgraphMatchingSurvey-main/SubgraphMatchingSurvey-main/vlabel/build2/utility/execution_tree && /home/konstantinos/miniconda3/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/konstantinos/SubgraphMatchingSurvey-main/SubgraphMatchingSurvey-main/vlabel/utility/execution_tree/execution_tree.cpp > CMakeFiles/et.dir/execution_tree.cpp.i

utility/execution_tree/CMakeFiles/et.dir/execution_tree.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/et.dir/execution_tree.cpp.s"
	cd /home/konstantinos/SubgraphMatchingSurvey-main/SubgraphMatchingSurvey-main/vlabel/build2/utility/execution_tree && /home/konstantinos/miniconda3/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/konstantinos/SubgraphMatchingSurvey-main/SubgraphMatchingSurvey-main/vlabel/utility/execution_tree/execution_tree.cpp -o CMakeFiles/et.dir/execution_tree.cpp.s

utility/execution_tree/CMakeFiles/et.dir/execution_tree.cpp.o.requires:

.PHONY : utility/execution_tree/CMakeFiles/et.dir/execution_tree.cpp.o.requires

utility/execution_tree/CMakeFiles/et.dir/execution_tree.cpp.o.provides: utility/execution_tree/CMakeFiles/et.dir/execution_tree.cpp.o.requires
	$(MAKE) -f utility/execution_tree/CMakeFiles/et.dir/build.make utility/execution_tree/CMakeFiles/et.dir/execution_tree.cpp.o.provides.build
.PHONY : utility/execution_tree/CMakeFiles/et.dir/execution_tree.cpp.o.provides

utility/execution_tree/CMakeFiles/et.dir/execution_tree.cpp.o.provides.build: utility/execution_tree/CMakeFiles/et.dir/execution_tree.cpp.o


utility/execution_tree/CMakeFiles/et.dir/execution_tree_generator.cpp.o: utility/execution_tree/CMakeFiles/et.dir/flags.make
utility/execution_tree/CMakeFiles/et.dir/execution_tree_generator.cpp.o: ../utility/execution_tree/execution_tree_generator.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/konstantinos/SubgraphMatchingSurvey-main/SubgraphMatchingSurvey-main/vlabel/build2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object utility/execution_tree/CMakeFiles/et.dir/execution_tree_generator.cpp.o"
	cd /home/konstantinos/SubgraphMatchingSurvey-main/SubgraphMatchingSurvey-main/vlabel/build2/utility/execution_tree && /home/konstantinos/miniconda3/bin/x86_64-conda-linux-gnu-c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/et.dir/execution_tree_generator.cpp.o -c /home/konstantinos/SubgraphMatchingSurvey-main/SubgraphMatchingSurvey-main/vlabel/utility/execution_tree/execution_tree_generator.cpp

utility/execution_tree/CMakeFiles/et.dir/execution_tree_generator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/et.dir/execution_tree_generator.cpp.i"
	cd /home/konstantinos/SubgraphMatchingSurvey-main/SubgraphMatchingSurvey-main/vlabel/build2/utility/execution_tree && /home/konstantinos/miniconda3/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/konstantinos/SubgraphMatchingSurvey-main/SubgraphMatchingSurvey-main/vlabel/utility/execution_tree/execution_tree_generator.cpp > CMakeFiles/et.dir/execution_tree_generator.cpp.i

utility/execution_tree/CMakeFiles/et.dir/execution_tree_generator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/et.dir/execution_tree_generator.cpp.s"
	cd /home/konstantinos/SubgraphMatchingSurvey-main/SubgraphMatchingSurvey-main/vlabel/build2/utility/execution_tree && /home/konstantinos/miniconda3/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/konstantinos/SubgraphMatchingSurvey-main/SubgraphMatchingSurvey-main/vlabel/utility/execution_tree/execution_tree_generator.cpp -o CMakeFiles/et.dir/execution_tree_generator.cpp.s

utility/execution_tree/CMakeFiles/et.dir/execution_tree_generator.cpp.o.requires:

.PHONY : utility/execution_tree/CMakeFiles/et.dir/execution_tree_generator.cpp.o.requires

utility/execution_tree/CMakeFiles/et.dir/execution_tree_generator.cpp.o.provides: utility/execution_tree/CMakeFiles/et.dir/execution_tree_generator.cpp.o.requires
	$(MAKE) -f utility/execution_tree/CMakeFiles/et.dir/build.make utility/execution_tree/CMakeFiles/et.dir/execution_tree_generator.cpp.o.provides.build
.PHONY : utility/execution_tree/CMakeFiles/et.dir/execution_tree_generator.cpp.o.provides

utility/execution_tree/CMakeFiles/et.dir/execution_tree_generator.cpp.o.provides.build: utility/execution_tree/CMakeFiles/et.dir/execution_tree_generator.cpp.o


# Object files for target et
et_OBJECTS = \
"CMakeFiles/et.dir/execution_tree_node.cpp.o" \
"CMakeFiles/et.dir/execution_tree.cpp.o" \
"CMakeFiles/et.dir/execution_tree_generator.cpp.o"

# External object files for target et
et_EXTERNAL_OBJECTS =

utility/execution_tree/libet.so: utility/execution_tree/CMakeFiles/et.dir/execution_tree_node.cpp.o
utility/execution_tree/libet.so: utility/execution_tree/CMakeFiles/et.dir/execution_tree.cpp.o
utility/execution_tree/libet.so: utility/execution_tree/CMakeFiles/et.dir/execution_tree_generator.cpp.o
utility/execution_tree/libet.so: utility/execution_tree/CMakeFiles/et.dir/build.make
utility/execution_tree/libet.so: utility/execution_tree/CMakeFiles/et.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/konstantinos/SubgraphMatchingSurvey-main/SubgraphMatchingSurvey-main/vlabel/build2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX shared library libet.so"
	cd /home/konstantinos/SubgraphMatchingSurvey-main/SubgraphMatchingSurvey-main/vlabel/build2/utility/execution_tree && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/et.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
utility/execution_tree/CMakeFiles/et.dir/build: utility/execution_tree/libet.so

.PHONY : utility/execution_tree/CMakeFiles/et.dir/build

utility/execution_tree/CMakeFiles/et.dir/requires: utility/execution_tree/CMakeFiles/et.dir/execution_tree_node.cpp.o.requires
utility/execution_tree/CMakeFiles/et.dir/requires: utility/execution_tree/CMakeFiles/et.dir/execution_tree.cpp.o.requires
utility/execution_tree/CMakeFiles/et.dir/requires: utility/execution_tree/CMakeFiles/et.dir/execution_tree_generator.cpp.o.requires

.PHONY : utility/execution_tree/CMakeFiles/et.dir/requires

utility/execution_tree/CMakeFiles/et.dir/clean:
	cd /home/konstantinos/SubgraphMatchingSurvey-main/SubgraphMatchingSurvey-main/vlabel/build2/utility/execution_tree && $(CMAKE_COMMAND) -P CMakeFiles/et.dir/cmake_clean.cmake
.PHONY : utility/execution_tree/CMakeFiles/et.dir/clean

utility/execution_tree/CMakeFiles/et.dir/depend:
	cd /home/konstantinos/SubgraphMatchingSurvey-main/SubgraphMatchingSurvey-main/vlabel/build2 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/konstantinos/SubgraphMatchingSurvey-main/SubgraphMatchingSurvey-main/vlabel /home/konstantinos/SubgraphMatchingSurvey-main/SubgraphMatchingSurvey-main/vlabel/utility/execution_tree /home/konstantinos/SubgraphMatchingSurvey-main/SubgraphMatchingSurvey-main/vlabel/build2 /home/konstantinos/SubgraphMatchingSurvey-main/SubgraphMatchingSurvey-main/vlabel/build2/utility/execution_tree /home/konstantinos/SubgraphMatchingSurvey-main/SubgraphMatchingSurvey-main/vlabel/build2/utility/execution_tree/CMakeFiles/et.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : utility/execution_tree/CMakeFiles/et.dir/depend

