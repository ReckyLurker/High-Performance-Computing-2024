# High-Performance-Computing-2024 - src

Contains the source code for linear solvers and sparse storage formats. The code uses Armadillo and OpenBLAS libraries for calculations. Currently limited to sequential instructions, will be upgraded to allow 
parallel instructions based on OpenMPI. 

1. matrix - Contains some overloads for the matrix formats offered by Armadillo 
2. solvers - Contains source code for linear solvers based on OpenBLAS and Armadillo