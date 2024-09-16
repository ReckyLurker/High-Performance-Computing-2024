#include <iostream>
#include "densematrix.hpp"
#include "loadvector.hpp"

int main() {

    // Loading into the space // 
    DenseMatrix A;
    A.load("../../matrix.txt"); // Load matrix from file 
    std::vector<double> B;
    loadvector("../../sourcevec.txt"); // Load vector from file
 
    return 0;
}
