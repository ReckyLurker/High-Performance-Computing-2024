#ifndef PRECONDITIONER_HPP 
#define PRECONDITIONER_HPP 

/*
    Implementation for different preconditioners. 
*/

#include <armadillo> 

void diagonal_preconditioner(const arma::Mat<double>& A){
    arma::Col<double> diagonal = A.diag();
    arma::Col<bool> zero_values = (diagonal == 0);

    size_t n_rows = A.n_rows; 
    size_t n_cols = A.n_cols;

    for(size_t i = 0; i < n_rows ; ++i){
        if(zero_values(i)) continue; // Skip zero values 
        double d_inv = diagonal(i);
        for(size_t j = 0; j < n_cols; ++j){
            
        }
    }
}

#endif 