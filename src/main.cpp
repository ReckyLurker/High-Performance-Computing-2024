#include <iostream> 
#include <armadillo>
#include "matrix.hpp"

#include <chrono> 

int main(void){

    // Load matrices based on fill in density //    

    arma::SpMat<double> A = loadCoordAscii("../matrix.txt");
    arma::Col<double> sourceVec;
    sourceVec.load("../sourcevec.txt", arma::raw_ascii); // Source Vector will always be in raw ascii form. 

    

}