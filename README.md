# High-Performance-Computing-2024
Codes for the assignments of the course High Performance Computing (CD61002) taken in Autumn 24-25

The "src" folder contains codes for linear solvers and some common overloads for dense and sparse matrices based on the "armadillo" and "openblas" libraries. To run the code on your machine, openblas and armadillo must be compiled from source. 

## Steps to build OpenBLAS and Armadillo from source:

1. Run `install_openblas.sh` to compile openblas from source. 
    - Please allow sudo access for installation.
    - Run `sudo bash install_openblas.sh` 

2. Run `install_armadillo.sh` to compile armadillo from source.
    - Please allow sudo access for installation. 
    - Run `sudo bash install_armadillo.sh` 

The above steps will perform a default installation for the armadillo and openblas libraries. Currently, this setup is only limited to Linux based OS. 