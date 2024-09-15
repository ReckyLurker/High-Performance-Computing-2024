# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# compile ASM with /usr/bin/cc
# compile C with /usr/bin/cc
ASM_DEFINES = -DUSE_GEMM3M

ASM_INCLUDES = -I/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src -I/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build

ASM_FLAGS =   -DHAVE_C11 -Wall -m64 -mavx2 -mavx -msse -msse2 -msse3 -mssse3 -msse4.1 -DF_INTERFACE_GFORT -fPIC -DSMALL_MATRIX_OPT -DSMP_SERVER -DNO_WARMUP -DMAX_CPU_NUMBER=8 -DMAX_PARALLEL_NUMBER=1 -DMAX_STACK_ALLOC=2048 -DNO_AFFINITY -DVERSION="\"0.3.28\""

C_DEFINES = -DUSE_GEMM3M

C_INCLUDES = -I/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src -I/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build

C_FLAGS =   -DHAVE_C11 -Wall -m64 -mavx2 -mavx -msse -msse2 -msse3 -mssse3 -msse4.1 -DF_INTERFACE_GFORT -fPIC -DSMALL_MATRIX_OPT -DSMP_SERVER -DNO_WARMUP -DMAX_CPU_NUMBER=8 -DMAX_PARALLEL_NUMBER=1 -DMAX_STACK_ALLOC=2048 -DNO_AFFINITY -DVERSION="\"0.3.28\"" -DBUILD_SINGLE -DBUILD_DOUBLE -DBUILD_COMPLEX -DBUILD_COMPLEX16

# Custom options: kernel/CMakeFiles/kernel.dir/CMakeFiles/srot_k.c.o_OPTIONS = -mfma

# Custom options: kernel/CMakeFiles/kernel.dir/CMakeFiles/drot_k.c.o_OPTIONS = -mfma

