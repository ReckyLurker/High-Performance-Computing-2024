#include <iostream> 
#include <armadillo>

#include "matrix.hpp" // Batchwise loading of sparse matrices
#include "diagonal_preconditioning.hpp" // Diagonal Preconditioning of sparse and dense matrices 
#include "jacobi.hpp" // Jacobi iterative method for sparse and dense matrices 

int main(void){

    arma::SpMat<double> A = loadCoordAscii("../matrix.txt");  // If given in sparse format (CSC)
    // arma::Mat<double> A; // If given in dense format
    // A.load("../matrix.txt");
    // Load source vector // 
    arma::Col<double> B;
    B.load("../sourcevec.txt", arma::raw_ascii); // Source Vector will always be in raw ascii form // 
    
    diagonal_preconditioning(A, B); // Diagonal (Jacobi) Preconditioning to reduce condition number    

    arma::Col<double> x_jacobi = jacobi(A, B, arma::zeros(A.n_rows), 1e-15, 50000); // Jacobi Method
}