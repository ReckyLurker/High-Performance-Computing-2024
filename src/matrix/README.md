# High-Performance-Computing-2024

This folder contains definitions for different matrix storage formats. The formats currently implemented: 

1. Dense Matrix (based on `std::vector<std::vector>`, requires $\mathcal{O}(N^2)$ storage space. 

2. Sparse Matrices: 
    - Compressed Sparse Row (CSR) / Compressed Row Storage (CRS)
    - Compressed Sparse Column (CSC) / Compressed Column Storage (CCS)
    - Coordinate List (COO) 
    - Diagonal Storage 
    - Dictonary of Keys (DOK) 
 
Overloads were developed for the storage formats offered by Armadillo. Currently limited to sequential machines. 