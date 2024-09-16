#ifndef DIAGONAL_PRECONDITIONING_HPP 
#define DIAGONAL_PRECONDITIONING_HPP

/*
    Source Code for Diagonal Preconditioning 
*/
#include <armadillo> 
#include <string>

#define RED "\033[31m" // Change Terminal color to Red 
#define YELLOW "\033[33m" // Change Terminal color to Yellow 
#define NC "\033[0m" // Reset terminal color 

#define VSMALL 1e-12 // Small Regularization term 

// For Dense Matrix //
void diagonal_preconditioning(arma::Mat<double>& A, arma::Col<double>& B){
    
    if(A.n_rows != A.n_cols){
        std::string error = std::string(RED) + std::string("Diagonal Preconditiong only for square matrices") + std::string(NC) + std::string("\n");
        throw std::runtime_error(error);
    }

    arma::Col<double> diag_inv = 1.0/(A.diag() + VSMALL); // Inverse Diagonal

    // Diagonal Preconditioning O(N^2) // 
    for(std::size_t i = 0; i < A.n_rows; ++i){
        for(std::size_t j = 0; j < A.n_cols; ++j){
            A(i,j) *= diag_inv(i); 
        }
        B(i) *= diag_inv(i);
    }
}

// For Sparse Matrix // 
void diagonal_preconditioning(arma::SpMat<double>& A, arma::Col<double>& B){
    if(A.n_rows != A.n_cols){
        std::string error = std::string(RED) + std::string("Diagonal Preconditiong only for square matrices") + std::string(NC) + std::string("\n");
        throw std::runtime_error(error);
    }

    arma::Col<double> diag_inv = 1.0/(A.diag() + VSMALL); // Inverse Diagonal

    // Diagonal Preconditioning - Sparse Matrix O(N) ops // 

    // Iterate through non-zero elements of A 
    for(arma::SpMat<double>::const_iterator itA = A.begin() ; itA != A.end() ; ++itA){
        std::size_t row_id = itA.row();
        std::size_t col_id = itA.col();
        A(row_id, col_id) *= diag_inv(row_id);
        B(row_id) *= diag_inv(row_id);
    }
}

#endif