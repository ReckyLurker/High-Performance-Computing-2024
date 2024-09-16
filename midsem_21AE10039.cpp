#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <sstream>

#define RED "\033[31m" // Change Terminal color to Red 
#define GREEN "\033[32m" // Change Terminal color to Green 
#define YELLOW "\033[33m" // Change Terminal color to Yellow 
#define NC "\033[0m" // Reset terminal color 

#define VSMALL 1e-12 // Small term for regularization
#define BATCH_SIZE 5 // Print after this many iterations

// Function Definitions - to be defined after main // 
std::vector<std::vector<double>> load(const char* filename, std::size_t& n_rows, std::size_t& n_cols);
std::vector<double> loadvector(const char* filename);
void diagonal_preconditioning(std::vector<std::vector<double>>& A, std::vector<double>& B);
double norm(const std::vector<double>& x, const std::vector<double>& x0); 
std::vector<double> jacobi(const std::vector<std::vector<double>>& A, const std::vector<double>& B, std::vector<double> x0, const double tol , const int max_iter); 
std::vector<double> gauss_seidel(const std::vector<std::vector<double>>& A, const std::vector<double>& B, std::vector<double> x0, const double tol , const int max_iter); 
std::vector<double> gauss_seidel_sor(const std::vector<std::vector<double>>& A, const std::vector<double>& B, const double omega, std::vector<double> x0, const double tol, const int max_iter);


int main() {

    // Loading into the space // 
    std::size_t n_rows = 0, n_cols = 0;
    std::vector<std::vector<double>> A = load("matrix.txt", n_rows, n_cols); // Load Dense Matrix into space // 
    std::vector<double> B = loadvector("sourcevec.txt"); 
    
    diagonal_preconditioning(A,B);

    // Solve Using Different Iterative Methods // 
    std::vector<double> x0(B.size(), 0.0);
    std::vector<double> x_jacobi = jacobi(A, B, x0, 1e-6, 1000); // jacobi(A, B, initial_guess, tol, max_iter);
    std::vector<double> x_gs = gauss_seidel(A, B, x0, 1e-6, 1000); // gauss_seidel(A, B, initial_guess, tol, max_iter)
    std::vector<double> x_gs_sor = gauss_seidel_sor(A, B, 1.5, x0, 1e-6, 1000); // gauss_seidel_sor(A, B, omega, initial_guess, tol, max_iter)
    return 0;
}


// Functions // 

// Load Dense Matrix // 
std::vector<std::vector<double>> load(const char* filename, std::size_t& n_rows, std::size_t& n_cols) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::string error = std::string(RED) + std::string(filename) + std::string(" couldn't be opened") + std::string(NC) + std::string("\n"); 
        throw std::runtime_error(error);
    }

    std::string line;
    std::vector<double> row; 
    double value; 

    n_rows = 0;
    n_cols = 0;
    std::vector<std::vector<double>> data; 

    while(std::getline(file, line)){
        std::istringstream iss(line);
        while(iss >> value){
            row.push_back(value);
        }
        if(n_cols == 0){
            n_cols = row.size();
        }
        if(n_cols != row.size()){
            std::string error = std::string(RED) + std::string("Incompatible File") + std::string(NC) + std::string("\n");
            throw std::runtime_error(error);
        }
        n_rows++;
        data.push_back(row);
        row.clear();
    }
    file.close();
    return data; // Update the Matrix 
}

// Load Source Vector from ASCII file // 
std::vector<double> loadvector(const char* filename){
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::string error = std::string(RED) + std::string(filename) + std::string(" couldn't be opened") + std::string(NC) + std::string("\n"); 
        throw std::runtime_error(error);
    }
 
    std::string line;
    std::vector<double> row; 
    double value; 

    while(std::getline(file, line)){
        std::istringstream iss(line);
        while(iss >> value){
            row.push_back(value);
        }
    }
    return row; // The source vector 
}

// Perform Diagonal Preconditioning to reduce condition number // 
void diagonal_preconditioning(std::vector<std::vector<double>>& A, std::vector<double>& B){
    std::size_t n_rows = B.size();
    if(A.size() == 0){
        std::string error = std::string(RED) + std::string("Empty Matrix A") + std::string(NC) + std::string("\n");
        throw std::runtime_error(error);
    }
    std::size_t n_cols = A[0].size();
    if(n_rows != n_cols){
        std::string error = std::string(RED) + std::string("Diagonal Preconditioning only for square matrices") + std::string(NC) + std::string("\n");
        throw std::runtime_error(error);
    }

    double D_inv; 
    for(std::size_t i = 0; i<n_rows; ++i){
        D_inv = 1.0/(A[i][i] + VSMALL);
        B[i] *= D_inv;
        for(std::size_t j = 0; j < n_cols ; ++j){
            A[i][j] *= D_inv;
        }
    }
}

// Returns the second norm of the vector // 
double norm(const std::vector<double>& x, const std::vector<double>& x0){
    if(x.size() != x0.size()){
        std::string error = std::string(RED) + std::string("Incompatible Vectors") + std::string(NC) + std::string("\n");
        throw std::runtime_error(error);
    }
    double sum = 0.0;
    for(std::size_t i = 0; i < x.size(); ++i){
        sum += (x[i] - x0[i])*(x[i] - x0[i]);
    }
    return std::sqrt(sum);
}

// Jacobi Iterative Method // 
std::vector<double> jacobi(const std::vector<std::vector<double>>& A, const std::vector<double>& B, std::vector<double> x0 = std::vector<double>(), const double tol = 1e-6, const int max_iter = 1000) {
    std::size_t n_rows = A.size();
    if(n_rows == 0){
        std::string error = std::string(RED) + std::string("Empty Matrix A") + std::string(NC) + std::string("\n");
        throw std::runtime_error(error);
    } 
    std::size_t n_cols = A[0].size();

    if(n_rows != n_cols){
        std::string error = std::string(RED) + std::string("Jacobi only for Square Matrices") + std::string(NC) + std::string("\n");
        throw std::runtime_error(error);
    }

    if(x0.size() == 0){
        x0 = std::vector<double>(n_rows, 0.0);
    }

    std::vector<double> D_inv(n_rows, 0.0);
    std::vector<std::vector<double>> R(n_rows, std::vector<double>(n_cols, 0.0)); // R = L + U 
    
    // Pre-process //
    for(std::size_t i = 0; i < n_rows ; ++i){
        D_inv[i] = 1.0/(A[i][i] + VSMALL);
        R[i][i] = -A[i][i];
        for(std::size_t j = 0; j < n_cols; ++j){
            R[i][j] += A[i][j];
        }
    }

    std::vector<double> x = x0;
    double residual = 1e6; 
    int iter = 0;
    int last_iter = 0; 

    while(residual > tol && iter <= max_iter){
        for(std::size_t i = 0; i < n_rows; ++i){
            x[i] = D_inv[i] * B[i];
            for(std::size_t j = 0 ; j < n_cols ; ++j){
                x[i] -= D_inv[i] * R[i][j] * x0[j]; // Matrix Vector Product 
            }
        }
        double residual = norm(x, x0);
        if(residual < tol){
            x0 = x; // Copy the new iteration to the old iteration
            std::cout << GREEN << "Jacobi Method -> Final Residual: " << residual << " Iterations: " << iter << NC << "\n";
            std::cout << GREEN << "Solution Converged" << NC << std::endl;
            return x;
        }
        if(iter >= last_iter + BATCH_SIZE){
            last_iter = iter; 
            std::cout << YELLOW << "Jacobi Method -> Residual: " << residual << " Iterations: " << iter << NC << "\n";
        }
        x0 = x; // Copy the new iteration to the old iteration 
        ++iter;
    }
    std::cout << RED << "Maximum Iterations reached." << NC << std::endl;
    return x;
}

// Gauss Seidel Iteration // 
std::vector<double> gauss_seidel(const std::vector<std::vector<double>>& A, const std::vector<double>& B, std::vector<double> x0 = std::vector<double>(), const double tol = 1e-6, const int max_iter = 1000){
    std::size_t n_rows = A.size();
    if(n_rows == 0){
        std::string error = std::string(RED) + std::string("Empty Matrix A") + std::string(NC) + std::string("\n");
        throw std::runtime_error(error);
    } 
    std::size_t n_cols = A[0].size();

    if(n_rows != n_cols){
        std::string error = std::string(RED) + std::string("Gauss Seidel only for Square Matrices") + std::string(NC) + std::string("\n");
        throw std::runtime_error(error);
    }

    if(x0.size() == 0){
        x0 = std::vector<double>(n_rows, 0.0);
    }

    std::vector<double> D_inv(n_rows, 0.0);
    
    // Pre-process //
    for(std::size_t i = 0; i < n_rows ; ++i){
        D_inv[i] = 1.0/(A[i][i] + VSMALL);
    }

    std::vector<double> x = x0;
    double residual = 1e6; 
    int iter = 0;
    int last_iter = 0; 

    while(residual > tol && iter <= max_iter){
        for(std::size_t i = 0; i < n_rows; ++i){
            x[i] = D_inv[i] * B[i];
            for(std::size_t j = 0 ; j < i ; ++j){
                x[i] -= D_inv[i] * A[i][j] * x[j]; // Matrix Vector Product (Lower Triangle)
            }
            for(std::size_t j = i+1; j < n_cols ; ++j){
                x[i] -= D_inv[i] * A[i][j] * x0[j]; // Matrix Vector Product (Upper Triangle)
            }
        }
        double residual = norm(x, x0);
        if(residual < tol){
            x0 = x; // Copy the new iteration to the old iteration
            std::cout << GREEN << "Gauss-Seidel Method -> Final Residual: " << residual << " Iterations: " << iter << NC << "\n";
            std::cout << GREEN << "Solution Converged" << NC << std::endl;
            return x;
        }
        if(iter >= last_iter + BATCH_SIZE){
            last_iter = iter; 
            std::cout << YELLOW << "Gauss-Seidel Method -> Residual: " << residual << " Iterations: " << iter << NC << "\n";
        }
        x0 = x; // Copy the new iteration to the old iteration 
        ++iter;
    }
    std::cout << RED << "Maximum Iterations reached." << NC << std::endl;
    return x;
}

// Gauss Seidel Successive Over Relaxation (SOR) Iteration // 
std::vector<double> gauss_seidel_sor(const std::vector<std::vector<double>>& A, const std::vector<double>& B, const double omega = 1.0, std::vector<double> x0 = std::vector<double>(), const double tol = 1e-6, const int max_iter = 1000){
    std::size_t n_rows = A.size();
    if(n_rows == 0){
        std::string error = std::string(RED) + std::string("Empty Matrix A") + std::string(NC) + std::string("\n");
        throw std::runtime_error(error);
    } 
    std::size_t n_cols = A[0].size();

    if(n_rows != n_cols){
        std::string error = std::string(RED) + std::string("Gauss Seidel only for Square Matrices") + std::string(NC) + std::string("\n");
        throw std::runtime_error(error);
    }

    if(x0.size() == 0){
        x0 = std::vector<double>(n_rows, 0.0);
    }

    std::vector<double> D_inv(n_rows, 0.0);
    
    // Pre-process //
    for(std::size_t i = 0; i < n_rows ; ++i){
        D_inv[i] = 1.0/(A[i][i] + VSMALL);
    }

    std::vector<double> x = x0;
    double residual = 1e6; 
    int iter = 0;
    int last_iter = 0; 

    double one_minus_omega = 1-omega; 
    while(residual > tol && iter <= max_iter){
        for(std::size_t i = 0; i < n_rows; ++i){
            x[i] = one_minus_omega * x0[i] + omega * D_inv[i] * B[i];
            for(std::size_t j = 0 ; j < i ; ++j){
                x[i] -= omega * D_inv[i] * A[i][j] * x[j]; // Matrix Vector Product (Lower Triangle)
            }
            for(std::size_t j = i+1; j < n_cols ; ++j){
                x[i] -= omega * D_inv[i] * A[i][j] * x0[j]; // Matrix Vector Product (Upper Triangle)
            }
        }
        double residual = norm(x, x0);
        if(residual < tol){
            x0 = x; // Copy the new iteration to the old iteration
            std::cout << GREEN << "Gauss-Seidel Method -> Final Residual: " << residual << " Iterations: " << iter << NC << "\n";
            std::cout << GREEN << "Solution Converged" << NC << std::endl;
            return x;
        }
        if(iter >= last_iter + BATCH_SIZE){
            last_iter = iter; 
            std::cout << YELLOW << "Gauss-Seidel Method -> Residual: " << residual << " Iterations: " << iter << NC << "\n";
        }
        x0 = x; // Copy the new iteration to the old iteration 
        ++iter;
    }
    std::cout << RED << "Maximum Iterations reached." << NC << std::endl;
    return x;
}