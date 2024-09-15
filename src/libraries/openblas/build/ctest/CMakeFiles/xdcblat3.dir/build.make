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
include ctest/CMakeFiles/xdcblat3.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include ctest/CMakeFiles/xdcblat3.dir/compiler_depend.make

# Include the progress variables for this target.
include ctest/CMakeFiles/xdcblat3.dir/progress.make

# Include the compile flags for this target's objects.
include ctest/CMakeFiles/xdcblat3.dir/flags.make

ctest/CMakeFiles/xdcblat3.dir/c_dblat3.f.o: ctest/CMakeFiles/xdcblat3.dir/flags.make
ctest/CMakeFiles/xdcblat3.dir/c_dblat3.f.o: /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/ctest/c_dblat3.f
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object ctest/CMakeFiles/xdcblat3.dir/c_dblat3.f.o"
	cd /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/ctest && /usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/ctest/c_dblat3.f -o CMakeFiles/xdcblat3.dir/c_dblat3.f.o

ctest/CMakeFiles/xdcblat3.dir/c_dblat3.f.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/xdcblat3.dir/c_dblat3.f.i"
	cd /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/ctest && /usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/ctest/c_dblat3.f > CMakeFiles/xdcblat3.dir/c_dblat3.f.i

ctest/CMakeFiles/xdcblat3.dir/c_dblat3.f.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/xdcblat3.dir/c_dblat3.f.s"
	cd /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/ctest && /usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/ctest/c_dblat3.f -o CMakeFiles/xdcblat3.dir/c_dblat3.f.s

ctest/CMakeFiles/xdcblat3.dir/c_dblas3.c.o: ctest/CMakeFiles/xdcblat3.dir/flags.make
ctest/CMakeFiles/xdcblat3.dir/c_dblas3.c.o: /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/ctest/c_dblas3.c
ctest/CMakeFiles/xdcblat3.dir/c_dblas3.c.o: ctest/CMakeFiles/xdcblat3.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object ctest/CMakeFiles/xdcblat3.dir/c_dblas3.c.o"
	cd /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/ctest && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT ctest/CMakeFiles/xdcblat3.dir/c_dblas3.c.o -MF CMakeFiles/xdcblat3.dir/c_dblas3.c.o.d -o CMakeFiles/xdcblat3.dir/c_dblas3.c.o -c /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/ctest/c_dblas3.c

ctest/CMakeFiles/xdcblat3.dir/c_dblas3.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/xdcblat3.dir/c_dblas3.c.i"
	cd /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/ctest && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/ctest/c_dblas3.c > CMakeFiles/xdcblat3.dir/c_dblas3.c.i

ctest/CMakeFiles/xdcblat3.dir/c_dblas3.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/xdcblat3.dir/c_dblas3.c.s"
	cd /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/ctest && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/ctest/c_dblas3.c -o CMakeFiles/xdcblat3.dir/c_dblas3.c.s

ctest/CMakeFiles/xdcblat3.dir/c_d3chke.c.o: ctest/CMakeFiles/xdcblat3.dir/flags.make
ctest/CMakeFiles/xdcblat3.dir/c_d3chke.c.o: /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/ctest/c_d3chke.c
ctest/CMakeFiles/xdcblat3.dir/c_d3chke.c.o: ctest/CMakeFiles/xdcblat3.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object ctest/CMakeFiles/xdcblat3.dir/c_d3chke.c.o"
	cd /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/ctest && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT ctest/CMakeFiles/xdcblat3.dir/c_d3chke.c.o -MF CMakeFiles/xdcblat3.dir/c_d3chke.c.o.d -o CMakeFiles/xdcblat3.dir/c_d3chke.c.o -c /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/ctest/c_d3chke.c

ctest/CMakeFiles/xdcblat3.dir/c_d3chke.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/xdcblat3.dir/c_d3chke.c.i"
	cd /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/ctest && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/ctest/c_d3chke.c > CMakeFiles/xdcblat3.dir/c_d3chke.c.i

ctest/CMakeFiles/xdcblat3.dir/c_d3chke.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/xdcblat3.dir/c_d3chke.c.s"
	cd /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/ctest && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/ctest/c_d3chke.c -o CMakeFiles/xdcblat3.dir/c_d3chke.c.s

ctest/CMakeFiles/xdcblat3.dir/auxiliary.c.o: ctest/CMakeFiles/xdcblat3.dir/flags.make
ctest/CMakeFiles/xdcblat3.dir/auxiliary.c.o: /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/ctest/auxiliary.c
ctest/CMakeFiles/xdcblat3.dir/auxiliary.c.o: ctest/CMakeFiles/xdcblat3.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object ctest/CMakeFiles/xdcblat3.dir/auxiliary.c.o"
	cd /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/ctest && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT ctest/CMakeFiles/xdcblat3.dir/auxiliary.c.o -MF CMakeFiles/xdcblat3.dir/auxiliary.c.o.d -o CMakeFiles/xdcblat3.dir/auxiliary.c.o -c /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/ctest/auxiliary.c

ctest/CMakeFiles/xdcblat3.dir/auxiliary.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/xdcblat3.dir/auxiliary.c.i"
	cd /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/ctest && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/ctest/auxiliary.c > CMakeFiles/xdcblat3.dir/auxiliary.c.i

ctest/CMakeFiles/xdcblat3.dir/auxiliary.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/xdcblat3.dir/auxiliary.c.s"
	cd /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/ctest && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/ctest/auxiliary.c -o CMakeFiles/xdcblat3.dir/auxiliary.c.s

ctest/CMakeFiles/xdcblat3.dir/c_xerbla.c.o: ctest/CMakeFiles/xdcblat3.dir/flags.make
ctest/CMakeFiles/xdcblat3.dir/c_xerbla.c.o: /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/ctest/c_xerbla.c
ctest/CMakeFiles/xdcblat3.dir/c_xerbla.c.o: ctest/CMakeFiles/xdcblat3.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object ctest/CMakeFiles/xdcblat3.dir/c_xerbla.c.o"
	cd /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/ctest && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT ctest/CMakeFiles/xdcblat3.dir/c_xerbla.c.o -MF CMakeFiles/xdcblat3.dir/c_xerbla.c.o.d -o CMakeFiles/xdcblat3.dir/c_xerbla.c.o -c /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/ctest/c_xerbla.c

ctest/CMakeFiles/xdcblat3.dir/c_xerbla.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/xdcblat3.dir/c_xerbla.c.i"
	cd /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/ctest && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/ctest/c_xerbla.c > CMakeFiles/xdcblat3.dir/c_xerbla.c.i

ctest/CMakeFiles/xdcblat3.dir/c_xerbla.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/xdcblat3.dir/c_xerbla.c.s"
	cd /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/ctest && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/ctest/c_xerbla.c -o CMakeFiles/xdcblat3.dir/c_xerbla.c.s

ctest/CMakeFiles/xdcblat3.dir/constant.c.o: ctest/CMakeFiles/xdcblat3.dir/flags.make
ctest/CMakeFiles/xdcblat3.dir/constant.c.o: /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/ctest/constant.c
ctest/CMakeFiles/xdcblat3.dir/constant.c.o: ctest/CMakeFiles/xdcblat3.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building C object ctest/CMakeFiles/xdcblat3.dir/constant.c.o"
	cd /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/ctest && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT ctest/CMakeFiles/xdcblat3.dir/constant.c.o -MF CMakeFiles/xdcblat3.dir/constant.c.o.d -o CMakeFiles/xdcblat3.dir/constant.c.o -c /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/ctest/constant.c

ctest/CMakeFiles/xdcblat3.dir/constant.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/xdcblat3.dir/constant.c.i"
	cd /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/ctest && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/ctest/constant.c > CMakeFiles/xdcblat3.dir/constant.c.i

ctest/CMakeFiles/xdcblat3.dir/constant.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/xdcblat3.dir/constant.c.s"
	cd /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/ctest && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/ctest/constant.c -o CMakeFiles/xdcblat3.dir/constant.c.s

# Object files for target xdcblat3
xdcblat3_OBJECTS = \
"CMakeFiles/xdcblat3.dir/c_dblat3.f.o" \
"CMakeFiles/xdcblat3.dir/c_dblas3.c.o" \
"CMakeFiles/xdcblat3.dir/c_d3chke.c.o" \
"CMakeFiles/xdcblat3.dir/auxiliary.c.o" \
"CMakeFiles/xdcblat3.dir/c_xerbla.c.o" \
"CMakeFiles/xdcblat3.dir/constant.c.o"

# External object files for target xdcblat3
xdcblat3_EXTERNAL_OBJECTS =

ctest/xdcblat3: ctest/CMakeFiles/xdcblat3.dir/c_dblat3.f.o
ctest/xdcblat3: ctest/CMakeFiles/xdcblat3.dir/c_dblas3.c.o
ctest/xdcblat3: ctest/CMakeFiles/xdcblat3.dir/c_d3chke.c.o
ctest/xdcblat3: ctest/CMakeFiles/xdcblat3.dir/auxiliary.c.o
ctest/xdcblat3: ctest/CMakeFiles/xdcblat3.dir/c_xerbla.c.o
ctest/xdcblat3: ctest/CMakeFiles/xdcblat3.dir/constant.c.o
ctest/xdcblat3: ctest/CMakeFiles/xdcblat3.dir/build.make
ctest/xdcblat3: lib/libopenblas.a
ctest/xdcblat3: ctest/CMakeFiles/xdcblat3.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Linking Fortran executable xdcblat3"
	cd /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/ctest && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/xdcblat3.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
ctest/CMakeFiles/xdcblat3.dir/build: ctest/xdcblat3
.PHONY : ctest/CMakeFiles/xdcblat3.dir/build

ctest/CMakeFiles/xdcblat3.dir/clean:
	cd /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/ctest && $(CMAKE_COMMAND) -P CMakeFiles/xdcblat3.dir/cmake_clean.cmake
.PHONY : ctest/CMakeFiles/xdcblat3.dir/clean

ctest/CMakeFiles/xdcblat3.dir/depend:
	cd /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/ctest /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/ctest /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/ctest/CMakeFiles/xdcblat3.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : ctest/CMakeFiles/xdcblat3.dir/depend

