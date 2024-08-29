#include <iostream> // For input-output functions, cin, cout 
#include <fstream> // For writing solution 
#include <string> // for using getline
#include <cmath> // For using sqrt 
#include <chrono> // For measuring execution time 
#include <iomanip> // For using setw - to set precision 

// Get number of variables using B.txt file // 
int get_number_points(const char* file_location){
    std::ifstream file(file_location);
    std::string line; 
    int total_lines = 0; 
    while(std::getline(file, line)){
        ++total_lines; 
    }
    file.close();
    return total_lines; 
}

// Read coefficient matrix from file // 
void read_A_from_file(const char* file_location, double **A, const int& NUM_POINTS){
    std::ifstream file(file_location);
    for(int i = 0; i<NUM_POINTS; ++i){
        for(int j = 0; j<NUM_POINTS; ++j){
            file >> A[i][j];
        }
    }
    file.close();
    return; 
}

// Read source from file // 
void read_B_from_file(const char* file_location, double *B, const int& NUM_POINTS){
    std::ifstream file(file_location);
    for(int i = 0; i<NUM_POINTS; ++i){
        file >> B[i];
    }
    file.close();
    return; 
}

// Diagonal Preconditioning (Jacobi Preconditioning) of the matrix // 
void diagonal_preconditioning(double **A, double *B, const int &NUM_POINTS){
    double temp = 0.0;
    for(int i = 0; i<NUM_POINTS; ++i){
        temp = 1.0/A[i][i];
        B[i] *= temp;
        for(int j = 0; j<NUM_POINTS; ++j){
            A[i][j] *= temp;
        }
    }
    return; 
}

// L2 residual of the linear system // 
double L2residual(double **A, double *B, double *x, const int& NUM_POINTS){
    double residual = 0.0, c = 0.0, sum_B=0.0; 
    for(int i = 0; i < NUM_POINTS; ++i){
        c = B[i];
        sum_B += B[i]*B[i];
        for(int j = 0; j < NUM_POINTS; ++j){
            c -= A[i][j]*x[j];
        }
        residual += c*c;
    }
    return std::sqrt(residual)/(std::sqrt(sum_B));
}

// Performs a single jacobi iteration // 
void perform_gs_sor_iteration(double **A, double *B, double *x, double *old_x, const int& NUM_POINTS, const double& omega){
    for(int i = 0; i<NUM_POINTS; ++i){
        x[i] = omega*B[i] + (1.0-omega)*old_x[i];
        for(int j = 0; j<NUM_POINTS; ++j){
            if(j == i){
                continue; 
            }
            else if(j > i){
                x[i] -= omega*A[i][j]*old_x[j];
            }
            else{ 
                x[i] -= omega*A[i][j]*x[j];
            }
        }
    }
}

// Update the old solution with the current solution // 
void update_old_solution(double *x, double *old_x, int NUM_POINTS){
    for(int i = 0; i<NUM_POINTS; ++i){
        old_x[i] = x[i];
    }
}

// Save solution at given file location // 
void save_solution(const char* file_location, double *x, int NUM_POINTS){
    std::ofstream file(file_location);
    for(int i = 0; i<NUM_POINTS; ++i){
        file << x[i] << " ";
    }
    file << "\n"; 
    file.close();
    return; 
}

// Reset solution // 
void reset_x(double *x, double *old_x, const int &NUM_POINTS){
    for(int i = 0; i<NUM_POINTS; ++i){
        x[i] = old_x[i] = 0.0; 
    }
}

// Solve using jacobi iteration // 
void solve_gs_sor(const char* file_location, double **A, double *B, double *x, double *old_x, double tolerance, const int& NUM_POINTS, const double& omega){
    double residual = 1e6; 
    const std::chrono::time_point<std::chrono::high_resolution_clock> t_start = std::chrono::high_resolution_clock::now();
    int num_iters = 0; 
    reset_x(x, old_x, NUM_POINTS);
    while(residual >= tolerance && num_iters <= 10000){
        perform_gs_sor_iteration(A, B, x, old_x, NUM_POINTS, omega);
        update_old_solution(x, old_x, NUM_POINTS);
        residual = L2residual(A, B, x, NUM_POINTS);
        // std::cout << "Residual: " << residual; 
        // std::cout << "Num_iters: " << num_iters+1 << "\n";
        ++num_iters; 
    }
    std::cout << "GS-SOR Method took: " << num_iters << " iterations for " << NUM_POINTS << " points with omega: " << omega << "\n";
    const std::chrono::time_point<std::chrono::high_resolution_clock> t_end = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> exec_time = t_end - t_start; 
    std::cout << "GS-SOR Method took: " << std::setw(9) << exec_time.count() << " seconds for " << NUM_POINTS << " points\n"; 
    save_solution(file_location, x, NUM_POINTS);
    return; 
}

int main(int argc, char** argv){
    
    // Change this to whatever you want // 
    double tolerance = 1e-6; 
    double omega = 1.8; 

    // Custom-Made function to count number of lines // 
    int NUM_POINTS = get_number_points("./finite_element/Fvec.txt");

    // Read Coefficient Matrix A and source B // 
    double **A = new double*[NUM_POINTS];
    for(int i = 0; i<NUM_POINTS; ++i){
        A[i] = new double[NUM_POINTS];
    }
    double *B = new double[NUM_POINTS];

    read_A_from_file("./finite_element/Kmat.txt", A, NUM_POINTS);
    read_B_from_file("./finite_element/Fvec.txt", B, NUM_POINTS);
    diagonal_preconditioning(A, B, NUM_POINTS);

    double *x = new double[NUM_POINTS];
    double *old_x = new double[NUM_POINTS];
    

    solve_gs_sor("./finite_element/solution_sor.txt", A, B, x, old_x, tolerance, NUM_POINTS, omega);

    /* Free Up memory taken by finite difference matrix and source vector */
    for(int i = 0; i<NUM_POINTS; ++i){
        delete[] A[i]; 
    }
    delete[] A;
    delete[] B;
    delete[] x, old_x; 
    return 0; 
}