#ifndef JACOBI_HPP
#define JACOBI_HPP 

/*
    Defines Jacobi Iterative Solver
*/

#include <armadillo> 
#include <string> 
#include <iostream>

#define RED "\033[31m" // Change Terminal color to Red 
#define GREEN "\033[32m" // Change Terminal color to Green 
#define YELLOW "\033[33m" // Change Terminal color to Yellow 
#define NC "\033[0m" // Reset terminal color 

#define VSMALL 1e-12 // Small term for regularization
#define BATCH_SIZE 5 // Print after this many iterations

// Jacobi Iterative Method for Dense matrices // 
arma::Col<double> jacobi(const arma::Mat<double>& A, const arma::Col<double>& B, arma::Col<double> x = arma::Col<double>(), const double tol = 1e-6, const int max_iter = 1000) {
    std::size_t n_rows = A.n_rows; 
    std::size_t n_cols = A.n_cols;

    if(n_rows != n_cols){
        std::string error = std::string(RED) + std::string("Jacobi only for Square Matrices") + std::string(NC) + std::string("\n");
        throw std::runtime_error(error);
    }

    if(x.n_rows == 0){
        x = arma::zeros(n_rows);
    }

    arma::Col<double> D_inv = 1.0/(A.diag() + VSMALL);
    arma::Mat<double> R = A - arma::diagmat(A.diag()); // R = L + U 
    
    int iter = 0;
    double residual = 1e6; 

    arma::Col<double> x_new = x;
    int last_iter = 0; 

    while(residual > tol && iter <= max_iter){
        x_new = D_inv % (B - R*x); // Update the solution vector 
        double residual = arma::norm(x_new - x, "inf");
        if(residual < tol){
            x = x_new; 
            std::cout << GREEN << "Solution Converged" << NC << std::endl;
            return x;
        }
        if(iter >= last_iter + BATCH_SIZE){
            last_iter = iter; 
            std::cout << YELLOW << "Jacobi Method -> Residual: " << residual << " Iterations: " << iter << NC << "\n";
        }
        x = x_new; 
        ++iter;
    }

    std::cout << RED << "Maximum Iterations reached." << NC << std::endl;
    return x;
}

// Jacobi method for sparse matrices // 
arma::Col<double> jacobi(const arma::SpMat<double>& A, const arma::Col<double>& B, arma::Col<double> x = arma::Col<double>(), const double tol = 1e-6, const int max_iter = 1000) {
    std::size_t n_rows = A.n_rows; 
    std::size_t n_cols = A.n_cols;

    if(n_rows != n_cols){
        std::string error = std::string(RED) + std::string("Jacobi only for Square Matrices") + std::string(NC) + std::string("\n");
        throw std::runtime_error(error);
    }

    if(x.n_rows == 0){
        x = arma::zeros(n_rows);
    }

    arma::Col<double> D_inv = 1.0/(A.diag() + VSMALL);
    arma::SpMat<double> R = A; // R = L + U 
    R.diag() = arma::zeros(A.n_rows);

    int iter = 0;
    double residual = 1e6; 

    arma::Col<double> x_new = x;
    int last_iter = 0; 

    while(residual > tol && iter <= max_iter){
        x_new = D_inv % (B - R*x); // Update the solution vector 
        double residual = arma::norm(x_new - x, "inf");
        if(residual < tol){
            x = x_new; 
            std::cout << GREEN << "Jacobi Method -> Final Residual: " << residual << " iterations: " << iter << NC << "\n";
            std::cout << GREEN << "Solution Converged" << NC << std::endl;
            return x;
        }
        if(iter >= last_iter + BATCH_SIZE){
            last_iter = iter; 
            std::cout << YELLOW << "Jacobi Method -> Residual: " << residual << " Iterations: " << iter << NC << "\n";
        }
        x = x_new; 
        ++iter;
    }

    std::cout << RED << "Maximum Iterations reached." << NC << std::endl;
    return x;
}

#endif