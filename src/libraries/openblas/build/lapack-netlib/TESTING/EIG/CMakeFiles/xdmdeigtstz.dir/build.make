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
include lapack-netlib/TESTING/EIG/CMakeFiles/xdmdeigtstz.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include lapack-netlib/TESTING/EIG/CMakeFiles/xdmdeigtstz.dir/compiler_depend.make

# Include the progress variables for this target.
include lapack-netlib/TESTING/EIG/CMakeFiles/xdmdeigtstz.dir/progress.make

# Include the compile flags for this target's objects.
include lapack-netlib/TESTING/EIG/CMakeFiles/xdmdeigtstz.dir/flags.make

lapack-netlib/TESTING/EIG/CMakeFiles/xdmdeigtstz.dir/zchkdmd.f90.o: lapack-netlib/TESTING/EIG/CMakeFiles/xdmdeigtstz.dir/flags.make
lapack-netlib/TESTING/EIG/CMakeFiles/xdmdeigtstz.dir/zchkdmd.f90.o: /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/lapack-netlib/TESTING/EIG/zchkdmd.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object lapack-netlib/TESTING/EIG/CMakeFiles/xdmdeigtstz.dir/zchkdmd.f90.o"
	cd /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/lapack-netlib/TESTING/EIG && /usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/lapack-netlib/TESTING/EIG/zchkdmd.f90 -o CMakeFiles/xdmdeigtstz.dir/zchkdmd.f90.o

lapack-netlib/TESTING/EIG/CMakeFiles/xdmdeigtstz.dir/zchkdmd.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/xdmdeigtstz.dir/zchkdmd.f90.i"
	cd /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/lapack-netlib/TESTING/EIG && /usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/lapack-netlib/TESTING/EIG/zchkdmd.f90 > CMakeFiles/xdmdeigtstz.dir/zchkdmd.f90.i

lapack-netlib/TESTING/EIG/CMakeFiles/xdmdeigtstz.dir/zchkdmd.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/xdmdeigtstz.dir/zchkdmd.f90.s"
	cd /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/lapack-netlib/TESTING/EIG && /usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/lapack-netlib/TESTING/EIG/zchkdmd.f90 -o CMakeFiles/xdmdeigtstz.dir/zchkdmd.f90.s

# Object files for target xdmdeigtstz
xdmdeigtstz_OBJECTS = \
"CMakeFiles/xdmdeigtstz.dir/zchkdmd.f90.o"

# External object files for target xdmdeigtstz
xdmdeigtstz_EXTERNAL_OBJECTS =

lapack-netlib/TESTING/EIG/xdmdeigtstz: lapack-netlib/TESTING/EIG/CMakeFiles/xdmdeigtstz.dir/zchkdmd.f90.o
lapack-netlib/TESTING/EIG/xdmdeigtstz: lapack-netlib/TESTING/EIG/CMakeFiles/xdmdeigtstz.dir/build.make
lapack-netlib/TESTING/EIG/xdmdeigtstz: lib/libopenblas.a
lapack-netlib/TESTING/EIG/xdmdeigtstz: lapack-netlib/TESTING/EIG/CMakeFiles/xdmdeigtstz.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking Fortran executable xdmdeigtstz"
	cd /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/lapack-netlib/TESTING/EIG && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/xdmdeigtstz.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
lapack-netlib/TESTING/EIG/CMakeFiles/xdmdeigtstz.dir/build: lapack-netlib/TESTING/EIG/xdmdeigtstz
.PHONY : lapack-netlib/TESTING/EIG/CMakeFiles/xdmdeigtstz.dir/build

lapack-netlib/TESTING/EIG/CMakeFiles/xdmdeigtstz.dir/clean:
	cd /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/lapack-netlib/TESTING/EIG && $(CMAKE_COMMAND) -P CMakeFiles/xdmdeigtstz.dir/cmake_clean.cmake
.PHONY : lapack-netlib/TESTING/EIG/CMakeFiles/xdmdeigtstz.dir/clean

lapack-netlib/TESTING/EIG/CMakeFiles/xdmdeigtstz.dir/depend:
	cd /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/lapack-netlib/TESTING/EIG /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/lapack-netlib/TESTING/EIG /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/lapack-netlib/TESTING/EIG/CMakeFiles/xdmdeigtstz.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lapack-netlib/TESTING/EIG/CMakeFiles/xdmdeigtstz.dir/depend

