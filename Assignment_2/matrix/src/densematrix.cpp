#include "densematrix.hpp"
#include <iomanip>
#include <string>
#include <iostream> 
#include <fstream>

#define RED "\033[31m" // Change Terminal color to Red 
#define YELLOW "\033[33m" // Change Terminal color to Yellow 
#define NC "\033[0m" // Reset terminal color 


#define BLOCK_SIZE 32

// Constructor
DenseMatrix::DenseMatrix() : Matrix(0,0), data(std::vector<std::vector<double>>()) {}
DenseMatrix::DenseMatrix(std::size_t rows, std::size_t cols) 
    : Matrix(rows, cols), data(rows, std::vector<double>(cols, 0.0)) {}

DenseMatrix::DenseMatrix(std::size_t rows, std::size_t cols, std::vector<std::vector<double>> m_data) 
    : Matrix(rows, cols), data(m_data) {}

// Destructor
DenseMatrix::~DenseMatrix() {}

// Set value in the matrix
void DenseMatrix::set(size_t row, size_t col, double value) {
    if (row >= 0 && row < n_rows && col >= 0 && col < n_cols) {
        data[row][col] = value;
        return; 
    }
    std::string error = std::string(RED) + std::string("Index Out of Bounds") + std::string(NC) + std::string("\n");
    throw std::out_of_range(error);
}

// Get value from the matrix
double DenseMatrix::get(size_t row, size_t col) const {
    if (row >= 0 && row < n_rows && col >= 0 && col < n_cols) {
        return data[row][col];
    }
    std::string error = std::string(RED) + std::string("Index Out of Bounds") + std::string(NC) + std::string("\n");
    throw std::out_of_range(error);
}

// Print the matrix //
void DenseMatrix::print() const {
    for (std::size_t i = 0; i < n_rows; ++i) {
        for (std::size_t j = 0; j < n_cols; ++j) {
            std::cout << std::setw(10) << data[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

void DenseMatrix::print(const char* print_string) const {
    std::cout << print_string << std::endl; 
    for (std::size_t i = 0; i < n_rows; ++i) {
        for (std::size_t j = 0; j < n_cols; ++j) {
            std::cout << std::setw(10) << data[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

// Load from ASCII file // 
void DenseMatrix::load(const char* filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::string error = std::string(RED) + std::string(filename) + std::string(" couldn't be opened") + std::string(NC) + std::string("\n"); 
        throw std::runtime_error(error);
    }

    std::string line;
    std::vector<double> row; 
    double value; 

    std::size_t n_rows = 0;
    std::size_t n_cols = 0;
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
    *this = DenseMatrix(n_rows, n_cols, data); // Update the Matrix 
}

std::vector<double> DenseMatrix::diag(){
    std::vector<double> diagonal; 
    if(n_rows != n_cols){
        std::string error = std::string(RED) + std::string("diag() only defined for square matrices") + std::string(NC) + std::string("\n");
        throw std::runtime_error(error);
    }
    for(std::size_t i = 0; i<n_rows; ++i){
        diagonal.push_back((*this)(i,i));
    }
    return diagonal;
}

// Operator Overloads // 
const double DenseMatrix::operator()(size_t row, size_t col) const{ 
    if(row >= 0 && row < n_rows && col >= 0 && col < n_cols){
        return data[row][col];
    }
    std::string error = std::string(RED) + std::string("Index Out of Bounds") + std::string(NC) + std::string("\n");
    throw std::out_of_range(error);
}

double& DenseMatrix::operator()(size_t row, size_t col){
    if(row >= 0 && row < n_rows && col >= 0 && col < n_cols){
        return data[row][col];
    }
    std::string error = std::string(RED) + std::string("Index Out of Bounds") + std::string(NC) + std::string("\n");
    throw std::out_of_range(error);
}

const std::vector<double>& DenseMatrix::operator[](size_t row) const {
    if(row >= 0 && row < n_rows){
        return data[row];
    }
    std::string error = std::string(RED) + std::string("Index Out of Bounds") + std::string(NC) + std::string("\n");
    throw std::out_of_range(error);
}

std::vector<double>& DenseMatrix::operator[](size_t row){
    if(row >= 0 && row < n_rows){
        return data[row];
    }
    std::string error = std::string(RED) + std::string("Index Out of Bounds") + std::string(NC) + std::string("\n");
    throw std::out_of_range(error);
}

// Scales as O(N^2) // 
DenseMatrix DenseMatrix::operator+(const DenseMatrix& B) const {
    std::size_t n_rows = this -> getRows();
    std::size_t n_cols = this -> getCols();
    
    if (n_rows != B.getRows() || n_cols != B.getCols()) {
        std::string error = std::string(RED) + std::string("Matrix dimensions incompatible") + std::string(NC) + std::string("\n");
        throw std::invalid_argument(error);
    }

    DenseMatrix result(n_rows, n_cols);

    for(std::size_t i = 0; i < n_rows ; ++i){
        for(std::size_t j = 0; j < n_cols; ++j){
            result(i,j) = (*this)(i,j) + B(i,j);
        }
    }
    return result;
}

// Implements Block Matrix Multiplication // 
DenseMatrix DenseMatrix::operator*(const DenseMatrix& B) const {
    std::size_t m = this->getRows();
    std::size_t k = this->getCols();
    if(k != B.getRows()){
        std::string error = std::string(RED) + std::string("Incompatible Matrices") + std::string(NC) + std::string("\n");
        throw std::invalid_argument(error);
    }
    std::size_t n = B.getCols();
    
    DenseMatrix C = DenseMatrix(m,n);

    // Cache-friendly block multiplication // 
    for(std::size_t iBlock = 0; iBlock < m; iBlock += BLOCK_SIZE){
        for(std::size_t jBlock = 0; jBlock < n; jBlock += BLOCK_SIZE){
            for(std::size_t kBlock = 0; kBlock < k ; kBlock += BLOCK_SIZE){
                 // Perform Block Multiplication // 
                for(std::size_t i = iBlock; i < std::min(iBlock + BLOCK_SIZE, m); ++i){
                    for(std::size_t j = jBlock; j < std::min(jBlock + BLOCK_SIZE, n); ++j){
                        double sum = 0.0;
                        for(std::size_t p = kBlock; p < std::min(kBlock + BLOCK_SIZE, k); ++p){
                            sum += (*this)(i,p) * B(p,j);
                        }
                        C(i,j) += sum;
                    }
                }
            }
        }
    }
    return C; 
}

