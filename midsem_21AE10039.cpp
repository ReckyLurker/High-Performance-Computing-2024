#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <unordered_map>

#define RED "\e[0;91m" // Change Terminal color to Red 
#define GREEN "\e[0;92m" // Change Terminal color to Green 
#define YELLOW "\e[0;93m" // Change Terminal color to Yellow 
#define NC "\e[0m" // Reset terminal color 

#define VSMALL 1e-12 // Small term for regularization
#define BATCH_SIZE 5 // Print after this many iterations (for iterative methods)

// Function Definitions - to be defined after main // 
std::vector<std::vector<double>> load(const char* filename, std::size_t& n_rows, std::size_t& n_cols);
void loadSparseMatrix(const char* filename, std::vector<double>& A_values, std::vector<std::size_t>& A_rows, std::vector<std::size_t>& A_cols, std::size_t& n_rows, std::size_t& n_cols);
std::vector<double> loadvector(const char* filename);
void diagonal_preconditioning(std::vector<std::vector<double>>& A, std::vector<double>& B);
void diagonal_preconditioning(std::vector<double>& A_values, std::vector<double>& B, const std::vector<std::size_t>& A_rows, const std::vector<std::size_t>& A_cols);
double norm(const std::vector<double>& x, const std::vector<double>& x0); 
std::vector<double> jacobi(const std::vector<std::vector<double>>& A, const std::vector<double>& B, std::vector<double> x0, const double tol , const int max_iter); 
std::vector<double> gauss_seidel(const std::vector<std::vector<double>>& A, const std::vector<double>& B, std::vector<double> x0, const double tol , const int max_iter); 
std::vector<double> gauss_seidel_sor(const std::vector<std::vector<double>>& A, const std::vector<double>& B, const double omega, std::vector<double> x0, const double tol, const int max_iter);
std::vector<double> cg(const std::vector<std::vector<double>>& A, const std::vector<double>& B, std::vector<double> x0, const double tol, const int max_iter);
std::vector<double> minres(const std::vector<std::vector<double>>& A, const std::vector<double>& B, std::vector<double> x0, const double tol, const int max_iter);
std::vector<double> steepest_descent(const std::vector<std::vector<double>>& A, const std::vector<double>& B, std::vector<double> x0, const double tol, const int max_iter);
std::vector<double> bicgstab(const std::vector<std::vector<double>>& A, const std::vector<double>& B, std::vector<double> x0, const double tol, const int max_iter);


int main() {

    std::size_t n_rows = 0, n_cols = 0;

    // Uncomment for loading dense matrix // 
    std::vector<std::vector<double>> A = load("matrix.txt", n_rows, n_cols); 

    // Uncomment for loading sparse matrix // 
    
    // std::vector<double> A_values;
    // std::vector<std::size_t> A_rows, A_cols; 
    // loadSparseMatrix("matrix.txt", A_values, A_rows, A_cols, n_rows, n_cols);

    // Load source vector //
    std::vector<double> B = loadvector("sourcevec.txt"); 
    
    diagonal_preconditioning(A,B); // For dense matrices 

    // diagonal_preconditioning(A_values, B, A_rows, A_cols); // For sparse matrices 

    // Solve Dense Linear System Using Different Iterative Methods // 
    std::vector<double> x0(B.size(), 0.0);
    std::vector<double> x0_temp = x0;
    std::vector<double> x_jacobi = jacobi(A, B, x0, 1e-8, 1000); // jacobi(A, B, initial_guess, tol, max_iter);
    x0 = x0_temp;
    std::vector<double> x_gs = gauss_seidel(A, B, x0, 1e-6, 1000); // gauss_seidel(A, B, initial_guess, tol, max_iter)
    x0 = x0_temp;
    std::vector<double> x_gs_sor = gauss_seidel_sor(A, B, 1.5, x0, 1e-6, 1000); // gauss_seidel_sor(A, B, omega, initial_guess, tol, max_iter)
    x0 = x0_temp;
    std::vector<double> x_cg = cg(A, B, x0, 1e-6, 1000); // cg(A, B, initial_guess, tol, max_iter)
    x0 = x0_temp;
    std::vector<double> x_minres = minres(A, B, x0, 1e-6, 1000); // minres(A, B, initial_guess, tol, max_iter)
    x0 = x0_temp;
    std::vector<double> x_steepest_descent = steepest_descent(A, B, x0, 1e-8, 1000); // steepest_descent(A, B, initial_guess, tol, max_iter)
    x0 = x0_temp;
    std::vector<double> x_bicgstab = bicgstab(A, B, x0, 1e-8, 1000); // bicgstab(A, B, initial_guess, tol, max_iter)
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

// Load Sparse Matrix in CSC format // 
void loadSparseMatrix(const char* filename, std::vector<double>& A_values, std::vector<std::size_t>& A_rows, std::vector<std::size_t>& A_cols, std::size_t& n_rows, std::size_t& n_cols){
    std::ifstream file(filename);
    
    std::size_t row, col; 
    double value;

    while(file >> row >> col >> value){
        A_rows.push_back(row);
        A_cols.push_back(col);
        A_values.push_back(value);
    }

    // Not correct theoretically, but should work for all matrices if complete diagonal provided //
    n_rows = *std::max_element(A_rows.begin(), A_rows.end());
    n_cols = *std::max_element(A_cols.begin(), A_cols.end());
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

// Perform diagonal preconditioning for sparse matrix in CSC format // 
void diagonal_preconditioning(std::vector<double>& A_values, std::vector<double>& B, const std::vector<std::size_t>& A_rows, const std::vector<std::size_t>& A_cols) {
    // Find the indices of diagonal elements // 
    std::unordered_map<std::size_t, std::size_t> diag_ids;  
    std::vector<double> diagonal; 
    for(std::size_t i = 0; i < A_rows.size() ; ++i){
        if(A_rows[i] == A_cols[i]){
            diag_ids.insert({A_rows[i], i}); // store for future use 
            diagonal.push_back(A_values[i]);
        }
    }

    // Invert Diagonal and scale down source vector B // 
    std::vector<double> diag_inv(diagonal.size(), 0.0); 
    for(std::size_t i = 0; i < diagonal.size(); ++i){
        diag_inv[i] = 1.0/(diagonal[i] + VSMALL);
    }

    // Scale down A // 
    for(std::size_t i = 0; i < A_rows.size(); ++i){
        std::size_t idx = diag_ids[A_rows[i]]; // Index of this diagonal - Not cache friendly
        A_values[i] *= diag_inv[idx];
    }

    // Scale down B // 
    for(std::size_t i = 0; i < A_rows.size(); ++i){
        std::size_t idx = diag_ids[i]; // Index of the diagonal 
        B[i] *= diag_inv[idx];
    }
}

// Returns the second norm of the vector subtraction x - x0 // 
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
        x0 = x; // Copy the new iteration to the old iteration 
        ++iter;
        if(iter >= last_iter + BATCH_SIZE){
            last_iter = iter; 
            std::cout << YELLOW << "Jacobi Method -> Residual: " << residual << " Iterations: " << iter << NC << "\n";
        }
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
        x0 = x; // Copy the new iteration to the old iteration 
        ++iter;
        if(iter >= last_iter + BATCH_SIZE){
            last_iter = iter; 
            std::cout << YELLOW << "Gauss-Seidel Method -> Residual: " << residual << " Iterations: " << iter << NC << "\n";
        }
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
        std::string error = std::string(RED) + std::string("Gauss Seidel SOR only for Square Matrices") + std::string(NC) + std::string("\n");
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
            std::cout << GREEN << "Gauss-Seidel SOR Method -> Final Residual: " << residual << " Iterations: " << iter << NC << "\n";
            std::cout << GREEN << "Solution Converged" << NC << std::endl;
            return x;
        }
        x0 = x; // Copy the new iteration to the old iteration 
        ++iter;
        if(iter >= last_iter + BATCH_SIZE){
            last_iter = iter; 
            std::cout << YELLOW << "Gauss-Seidel SOR Method -> Residual: " << residual << " Iterations: " << iter << NC << "\n";
        }
    }
    std::cout << RED << "Maximum Iterations reached." << NC << std::endl;
    return x;
}

// Krylov SUBSPACE METHODS // 

// Matrix-Vector Multiplication // 
std::vector<double> matvecmul(const std::vector<std::vector<double>>& A, const std::vector<double> B){
    std::size_t n_rows = A.size();
    if(n_rows == 0){
        std::string error = std::string(RED) + std::string("Empty Matrix A") + std::string(NC) + std::string("\n");
        throw std::runtime_error(error);
    } 
    std::size_t n_cols = A[0].size();
    if(n_cols != B.size()){
        std::string error = std::string(RED) + std::string("Incompatible Matrix Vector Product") + std::string(NC) + std::string("\n");
        throw std::runtime_error(error);
    }
    std::vector<double> res(n_rows, 0.0);
    for(std::size_t i = 0; i < n_rows; ++i){
        for(std::size_t j = 0 ; j < n_cols ; ++j){
            res[i] += A[i][j] * B[j]; // Matrix vector product
        }
    }
    return res; 
}

// Dot Product (Vector) // 
double dot(const std::vector<double>& x, const std::vector<double>& y){
    if(x.size() != y.size()){
        std::string error = std::string(RED) + std::string("Incompatible vectors for dot Product") + std::string(NC) + std::string("\n");
        throw std::runtime_error(error);
    }
    double sum = 0.0;
    for(std::size_t i = 0; i < x.size(); ++i){
        sum += x[i]*y[i];
    }
    return sum;
}

// Dot Product (Matrix Vector) // 
std::vector<double> dot(const std::vector<std::vector<double>>& A, const std::vector<double>& B){
    std::size_t n_rows = A.size();
    if(n_rows == 0){
        std::string error = std::string(RED) + std::string("Empty Matrix A") + std::string(NC) + std::string("\n");
        throw std::runtime_error(error);
    }
    if(n_rows != B.size()){
        std::string error = std::string(RED) + std::string("Incompatible matrix-vector for dot Product") + std::string(NC) + std::string("\n");
        throw std::runtime_error(error);
    }
    std::size_t n_cols = A[0].size();

    std::vector<double> res(n_cols, 0.0);
    for(std::size_t i = 0; i < n_rows; ++i){
        for(std::size_t j = 0; j < n_cols ; ++j){
            res[j] += A[i][j] * B[j]; // Cache Friendly access
        }
    }
    return res;
}

// Add (vector + vector) // 
std::vector<double> add(const std::vector<double>& x, const std::vector<double>& y){
    if(x.size() != y.size()){
        std::string error = std::string(RED) + std::string("Incompatible vectors for addition") + std::string(NC) + std::string("\n");
        throw std::runtime_error(error);
    }
    std::vector<double> res(x.size());
    for(std::size_t i = 0; i < x.size(); ++i){
        res[i] = x[i] + y[i];
    }
    return res;
}

// Subtract (vector - vector) // 
std::vector<double> subtract(const std::vector<double>& x, const std::vector<double>& y){
    if(x.size() != y.size()){
        std::string error = std::string(RED) + std::string("Incompatible vectors for subtraction") + std::string(NC) + std::string("\n");
        throw std::runtime_error(error);
    }
    std::vector<double> res(x.size());
    for(std::size_t i = 0; i < x.size(); ++i){
        res[i] = x[i] - y[i];
    }
    return res;
}

// Scale a vector by a scalar // 
std::vector<double> scale(const std::vector<double>& x, const double& alpha){
    std::vector<double> res(x.size());
    for(std::size_t i = 0; i < x.size(); ++i){
        res[i] = x[i] * alpha;
    }
    return res;
}

// Conjugate Gradient Iterative Method // 
std::vector<double> cg(const std::vector<std::vector<double>>& A, const std::vector<double>& B, std::vector<double> x0 = std::vector<double>(), const double tol = 1e-6, const int max_iter = 1000) {
    std::size_t n_rows = A.size();
    if(n_rows == 0){
        std::string error = std::string(RED) + std::string("Empty Matrix A") + std::string(NC) + std::string("\n");
        throw std::runtime_error(error);
    } 
    std::size_t n_cols = A[0].size();

    if(n_rows != n_cols){
        std::string error = std::string(RED) + std::string("Conjugate Gradient only for Square Matrices") + std::string(NC) + std::string("\n");
        throw std::runtime_error(error);
    }

    if(x0.size() == 0){
        x0 = std::vector<double>(n_rows, 0.0);
    }

    std::vector<double> x = x0; 
    std::vector<double> res = subtract(B, matvecmul(A,x));
    std::vector<double> p = res; 
    std::vector<double> Ap = matvecmul(A,p);

    int iter = 0; 
    int last_iter = 0; 

    double rs_old = dot(res,res);
    double residual = 1e6;

    while(residual > tol && iter <= max_iter){
        Ap = matvecmul(A, p);
        rs_old = dot(res,res);
        double alpha = rs_old/(dot(Ap,p) + VSMALL);
        x = add(x, scale(p, alpha));
        res = subtract(res, scale(Ap, alpha));
        residual = norm(x,x0);
        x0 = x;
        if(residual < tol){
            std::cout << GREEN << "Conjugate Gradient Method -> Final Residual: " << residual << " Iterations: " << iter << NC << "\n";
            std::cout << GREEN << "Solution Converged" << NC << std::endl;
            return x;
        }
        p = add(res, scale(p, dot(res,res) / (rs_old + VSMALL)));
        ++iter;
        if(iter >= last_iter + BATCH_SIZE){
            last_iter = iter; 
            std::cout << YELLOW << "Conjugate Gradient Method -> Residual: " << residual << " Iterations: " << iter << NC << "\n";
        } 
    }
    std::cout << RED << "Maximum Iterations reached or solution couldn't converge." << NC << std::endl;
    return x;
}

// Minimum Residual Iterations (MINRES) // 
std::vector<double> minres(const std::vector<std::vector<double>>& A, const std::vector<double>& B, std::vector<double> x0 = std::vector<double>(), const double tol = 1e-6, const int max_iter = 1000) {
    std::size_t n_rows = A.size();
    if(n_rows == 0){
        std::string error = std::string(RED) + std::string("Empty Matrix A") + std::string(NC) + std::string("\n");
        throw std::runtime_error(error);
    } 
    std::size_t n_cols = A[0].size();

    if(n_rows != n_cols){
        std::string error = std::string(RED) + std::string("Minimum Residual only for Square Matrices") + std::string(NC) + std::string("\n");
        throw std::runtime_error(error);
    }

    if(x0.size() == 0){
        x0 = std::vector<double>(n_rows, 0.0);
    }

    std::vector<double> x = x0; 
    std::vector<double> res = subtract(B, matvecmul(A,x));
    std::vector<double> p = matvecmul(A,res);

    int iter = 0; 
    int last_iter = 0;
    double residual = 1e6; 

    while(residual > tol && iter <= max_iter){
        double alpha = dot(p,res) / (dot(p,p) + VSMALL);
        x = add(x, scale(res, alpha));
        residual = norm(x,x0);
        x0 = x;
        res = subtract(res, scale(p, alpha));
        p = matvecmul(A, res);

        if(residual < tol){
            std::cout << GREEN << "Minimum Residual Method -> Final Residual: " << residual << " Iterations: " << iter << NC << "\n";
            std::cout << GREEN << "Solution Converged" << NC << std::endl;
            return x;
        }
        ++iter;
        if(iter >= last_iter + BATCH_SIZE){
            last_iter = iter; 
            std::cout << YELLOW << "Minimum Residual Method -> Residual: " << residual << " Iterations: " << iter << NC << "\n";
        } 
    }
    std::cout << RED << "Maximum Iterations reached or solution couldn't converge." << NC << std::endl;
    return x;
}

// Steepest Descent Iteration // 
std::vector<double> steepest_descent(const std::vector<std::vector<double>>& A, const std::vector<double>& B, std::vector<double> x0 = std::vector<double>(), const double tol = 1e-6, const int max_iter = 1000){
    std::size_t n_rows = A.size();
    if(n_rows == 0){
        std::string error = std::string(RED) + std::string("Empty Matrix A") + std::string(NC) + std::string("\n");
        throw std::runtime_error(error);
    } 
    std::size_t n_cols = A[0].size();

    if(n_rows != n_cols){
        std::string error = std::string(RED) + std::string("Steepest Descent only for Square Matrices") + std::string(NC) + std::string("\n");
        throw std::runtime_error(error);
    }

    if(x0.size() == 0){
        x0 = std::vector<double>(n_rows, 0.0);
    }

    std::vector<double> x = x0; 
    std::vector<double> res = subtract(B, matvecmul(A,x));
    std::vector<double> v = dot(A,res); // equivalent to A.t * res
    std::vector<double> Av = matvecmul(A,v);

    int iter = 0; 
    int last_iter = 0;
    double residual = 1e6; 

    while(residual > tol && iter <= max_iter){
        v = dot(A,res);
        Av = matvecmul(A,v);
        double alpha = dot(v,v)/(dot(Av,Av) + VSMALL);
        x = add(x, scale(v, alpha));
        residual = norm(x,x0);
        x0 = x;
        res = subtract(res, scale(Av, alpha));

        if(residual < tol){
            std::cout << GREEN << "Steepest Descent Method -> Final Residual: " << residual << " Iterations: " << iter << NC << "\n";
            std::cout << GREEN << "Solution Converged" << NC << std::endl;
            return x;
        }
        ++iter;
        if(iter >= last_iter + BATCH_SIZE){
            last_iter = iter; 
            std::cout << YELLOW << "Steepest Descent Method -> Residual: " << residual << " Iterations: " << iter << NC << "\n";
        } 
    }
    std::cout << RED << "Maximum Iterations reached or solution couldn't converge." << NC << std::endl;
    return x;
}

// BiCGSTAB Iteration // 
std::vector<double> bicgstab(const std::vector<std::vector<double>>& A, const std::vector<double>& B, std::vector<double> x0 = std::vector<double>(), const double tol = 1e-6, const int max_iter = 1000){
    std::size_t n_rows = A.size();
    if(n_rows == 0){
        std::string error = std::string(RED) + std::string("Empty Matrix A") + std::string(NC) + std::string("\n");
        throw std::runtime_error(error);
    } 
    std::size_t n_cols = A[0].size();

    if(n_rows != n_cols){
        std::string error = std::string(RED) + std::string("BICGSTAB only for Square Matrices") + std::string(NC) + std::string("\n");
        throw std::runtime_error(error);
    }

    if(x0.size() == 0){
        x0 = std::vector<double>(n_rows, 0.0);
    }

    std::vector<double> x = x0; 
    std::vector<double> r = subtract(B, matvecmul(A,x));
    std::vector<double> r_hat = r; 
    std::vector<double> p = r;
    
    std::vector<double> Ap = matvecmul(A,p);
    std::vector<double> s = p; 
    std::vector<double> As = matvecmul(A,s);

    int iter = 0; 
    int last_iter = 0;
    double residual = 1e6; 

    while(residual > tol && iter <= max_iter){
        p = r;
        Ap = matvecmul(A,p);
        double alpha = dot(r, r_hat)/(dot(Ap, r_hat) + VSMALL);
        s = subtract(r, scale(Ap, alpha));
        As = matvecmul(A,s);
        double omega = dot(As, s)/(dot(As,As) + VSMALL);
        x = add(add(x, scale(p, alpha)), scale(s,omega));
        double temp = dot(r, r_hat);
        r = subtract(s, scale(As, omega));
        double beta = dot(r, r_hat)/(temp + VSMALL) * (alpha/(omega + VSMALL));
        p = add(r, scale(subtract(p, scale(Ap, omega)), beta));

        residual = norm(x, x0);
        x0 = x;
        if(residual < tol){
            std::cout << GREEN << "BICGSTAB Method -> Final Residual: " << residual << " Iterations: " << iter << NC << "\n";
            std::cout << GREEN << "Solution Converged" << NC << std::endl;
            return x;
        }
        ++iter;
        if(iter >= last_iter + BATCH_SIZE){
            last_iter = iter; 
            std::cout << YELLOW << "BICGSTAB Method -> Residual: " << residual << " Iterations: " << iter << NC << "\n";
        } 
    }
    std::cout << RED << "Maximum Iterations reached or solution couldn't converge." << NC << std::endl;
    return x;
}