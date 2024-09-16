/*
    Jacobi Solver. 
    Depends on Armadillo and OpenBLAS. 
*/

#include <iostream> 
#include <armadillo> 
#include "preconditioner.hpp"

arma::Col<double> jacobi(const arma::SpMat<double>& A, const arma::Col<double>& B){

    // Apply Diagonal Preconditioning // 
    diagonal_preconditioner()
}
