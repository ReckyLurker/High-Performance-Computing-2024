/* 
    Assignment-2 
    Author: Tanish Jain (recklurker)
    Roll Number: 21AE10039
    Generates Finite Difference Matrices for the Laplacian Problem
*/

#include <iostream> 
#include <armadillo> 

#define YELLOW "\033[33m" // Change Terminal color to Yellow 
#define NC "\033[0m" // Reset terminal color 

void generateDenseMatrix(size_t nx, size_t ny, double dx, double dy){
    size_t mat_size = (nx+1)*(ny+1); 
    double dx_inv = 1.0/dx;
    double dy_inv = 1.0/dy; 
    arma::Mat<double> A = arma::zeros<arma::Mat<double>>(mat_size, mat_size);

    // Boundary Condition - Bottom (y = 0)
    for(size_t i = 0; i<nx+1; ++i){
        A(i*(ny+1), i*(ny+1)) = 1.0; // j = 0 
    }

    // Boundary Condition - Left (x = 0)
    for(size_t j = 0; j<ny+1; ++j){
        A(j, j) = 1.0; // i = 0
    }

    // Boundary Condition - Right (x = Lx)
    for(size_t j = 0; j<ny+1; ++j){
        A(nx*(ny+1) + j, nx*(ny+1) + j) = 1.0; // i = nx 
    }

    // Boundary Condition - Top (y = Ly)
    for(size_t i = 0; i < nx+1 ; ++i){
        A(i*(ny+1) + ny, i*(ny+1) + ny) = 1.0; // j = ny 
    }

    // For the rest of the elements - 5 Point Stencil // 
    for(size_t i = 1; i<nx; ++i){
        for(size_t j = 1; j < ny; ++j){
            size_t k = i*(ny+1) + j; 
            A(k, k) = 2*dx_inv*dx_inv + 2*dy_inv*dy_inv; 
            A(k, k+1) = -dy_inv*dy_inv; 
            A(k, k-1) = -dy_inv*dy_inv; 
            A(k, k + ny+1) = -dx_inv*dx_inv; 
            A(k, k - (ny+1)) = -dx_inv*dx_inv; 
        }
    }
    A.save("../matrix.txt", arma::raw_ascii);
}

void generateSparseMatrix(size_t nx, size_t ny, double dx, double dy){
    size_t mat_size = (nx+1)*(ny+1);
    double dx_inv = 1.0/dx;
    double dy_inv = 1.0/dy; 

    arma::SpMat<double> A(mat_size, mat_size);

    // Boundary Condition - Bottom (y = 0)
    for(size_t i = 0; i<nx+1; ++i){
        A(i*(ny+1), i*(ny+1)) = 1.0; // j = 0 
    }

    // Boundary Condition - Left (x = 0)
    for(size_t j = 0; j<ny+1; ++j){
        A(j, j) = 1.0; // i = 0
    }

    // Boundary Condition - Right (x = Lx)
    for(size_t j = 0; j<ny+1; ++j){
        A(nx*(ny+1) + j, nx*(ny+1) + j) = 1.0; // i = nx 
    }

    // Boundary Condition - Top (y = Ly)
    for(size_t i = 0; i < nx+1 ; ++i){
        A(i*(ny+1) + ny, i*(ny+1) + ny) = 1.0; // j = ny 
    }

    // For the rest of the elements - 5 Point Stencil // 
    for(size_t i = 1; i<nx; ++i){
        for(size_t j = 1; j < ny; ++j){
            size_t k = i*(ny+1) + j; 
            A(k, k) = 2*dx_inv*dx_inv + 2*dy_inv*dy_inv; 
            A(k, k+1) = -dy_inv*dy_inv; 
            A(k, k-1) = -dy_inv*dy_inv; 
            A(k, k + ny+1) = -dx_inv*dx_inv; 
            A(k, k - (ny+1)) = -dx_inv*dx_inv; 
        }
    }
    A.save("../matrix.txt", arma::coord_ascii);
}

void generateSourceVector(size_t nx, size_t ny){
    size_t vec_size = (nx+1)*(ny+1);

    arma::Col<double> sourceVector = arma::zeros<arma::Col<double>>(vec_size);

    for(size_t i = 0; i < nx+1 ; ++i){
        sourceVector(i*(ny+1) + ny) = 1.0; // T = 1 on top  
    }

    sourceVector.save("../sourcevec.txt", arma::raw_ascii);
}


int main(void){
    // Basic IO 

    double dx, dy; // Step Sizes 
    double Lx, Ly; // Grid Sizes
    size_t nx, ny; // Number of points in Grid 

    std::cout << "Enter dx: "; 
    std::cin >> dx;
    std::cout << "Enter dy: ";
    std::cin >> dy; 
    std::cout << "Enter Lx: ";
    std::cin >> Lx;
    std::cout << "Enter Ly: ";
    std::cin >> Ly; 

    nx = (size_t)(Lx/dx);
    ny = (size_t)(Ly/dy);

    std::cout << YELLOW << "Number of Elements in X: " << nx << NC << "\n";
    std::cout << YELLOW << "Number of Elements in Y: " << ny << NC << "\n";

    if(nx*ny <= 100){
        generateDenseMatrix(nx, ny, dx, dy); // Generate A in dense format.
    }
    else{
        generateSparseMatrix(nx, ny, dx, dy); // Generate A in sparse format, avoid overheads.
    }   

    generateSourceVector(nx, ny); // Generate B 

    return 0; 
}
