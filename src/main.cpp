#include <iostream> 
#include <armadillo>
#include "matrix.hpp"


int main(void){

    // Load matrices based on fill in density //    

    arma::SpMat<double> A = loadCoordAscii("../matrix.txt");
    arma::Col<double> sourceVec;
    sourceVec.load("../sourcevec.txt", arma::raw_ascii); // Source Vector will always be in raw ascii form. 
    std::cout << A.n_rows << " " << A.n_cols << "\n";
    

}