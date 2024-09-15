# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build

# Include any dependencies generated for this target.
include test/CMakeFiles/zblat3.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include test/CMakeFiles/zblat3.dir/compiler_depend.make

# Include the progress variables for this target.
include test/CMakeFiles/zblat3.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/zblat3.dir/flags.make

test/CMakeFiles/zblat3.dir/zblat3.f.o: test/CMakeFiles/zblat3.dir/flags.make
test/CMakeFiles/zblat3.dir/zblat3.f.o: /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/test/zblat3.f
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object test/CMakeFiles/zblat3.dir/zblat3.f.o"
	cd /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/test && /usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/test/zblat3.f -o CMakeFiles/zblat3.dir/zblat3.f.o

test/CMakeFiles/zblat3.dir/zblat3.f.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/zblat3.dir/zblat3.f.i"
	cd /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/test && /usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/test/zblat3.f > CMakeFiles/zblat3.dir/zblat3.f.i

test/CMakeFiles/zblat3.dir/zblat3.f.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/zblat3.dir/zblat3.f.s"
	cd /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/test && /usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/test/zblat3.f -o CMakeFiles/zblat3.dir/zblat3.f.s

# Object files for target zblat3
zblat3_OBJECTS = \
"CMakeFiles/zblat3.dir/zblat3.f.o"

# External object files for target zblat3
zblat3_EXTERNAL_OBJECTS =

test/zblat3: test/CMakeFiles/zblat3.dir/zblat3.f.o
test/zblat3: test/CMakeFiles/zblat3.dir/build.make
test/zblat3: lib/libopenblas.a
test/zblat3: test/CMakeFiles/zblat3.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking Fortran executable zblat3"
	cd /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/zblat3.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/zblat3.dir/build: test/zblat3
.PHONY : test/CMakeFiles/zblat3.dir/build

test/CMakeFiles/zblat3.dir/clean:
	cd /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/test && $(CMAKE_COMMAND) -P CMakeFiles/zblat3.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/zblat3.dir/clean

test/CMakeFiles/zblat3.dir/depend:
	cd /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/test /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/test /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/test/CMakeFiles/zblat3.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/zblat3.dir/depend

